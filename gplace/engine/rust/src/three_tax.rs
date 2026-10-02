//! The three-taxon tree of one query at one edge, and its branch-length optimisation
//! (pplacer's Three_tax.optimize with Minimization.brent on GSL 1.16's Brent minimiser).
//!
//! Memory-lean form, with the same values as pplacer's:
//!  - the pendant (query) edge: a query site vector is one of nu <= 64 distinct vectors and evolve()
//!    is per site, so the evolved query at site s is qp[cat[s]][qc[s]], where qp[r][u] = evolve(u[u])
//!    under P_r(pendant): identical arithmetic, computed for nu vectors instead of m sites;
//!  - the distal/proximal edges: pd[s][st] = (statd[st] * prox[s][st]) * dist[s][st], the left part
//!    of mat_log_like3's statd * x * y * z (x = prox, y = dist, z = pendant), same association;
//!  - log-likelihood: sum_s log(sum_st pd[s][st] * z[st]) in site order, + log2 * (exponent total).
//!
//! pplacer quirks reproduced: the tree keeps the LAST evaluated branch lengths (not Brent's
//! minimum); the bisection start finder evaluates f(right) before f(left); a failing
//! `assert(prox_bl >= 0.)` inside a GSL callback returns NaN (ocamlgsl's caml_callback_exn), so GSL
//! reports an error and the caller retries from the default branch lengths.

use crate::ffi::{self, log, GslFunction};
use crate::glv::{evolve_site, xp_total_twoexp};
use crate::model::xp_multi_exp;
use crate::{simd, sl, NS, SIMD20};
use std::f64::consts::LN_2;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::os::raw::{c_char, c_int, c_void};
use std::ptr;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

/// Status codes returned to Python (gplace/engine/core.py).
pub const XP_OK: c_int = 0;
pub const XP_EXCEEDED: c_int = 1;
pub const XP_GSL: c_int = 2;
pub const XP_FINDSTART: c_int = 3;
pub const XP_INVALIDSTART: c_int = 4;
pub const XP_ASSERT: c_int = 5;

/// Most distinct query site vectors (and gamma rates) per three-taxon tree.
const MAX_DISTINCT: usize = 64;

/// C's `!(x >= 0.)`, the failing case of pplacer's `assert(prox_bl >= 0.)`: x is negative or NaN.
#[inline(always)]
fn negative_or_nan(x: f64) -> bool {
    x < 0. || x.is_nan()
}

pub struct Tt {
    m: usize,               // number of (query-informative) masked sites
    cat: *const i32,        // site category per masked site (CAT); unused for gamma
    statd: *const f64,
    x: *const f64,
    l: *const f64,
    xit: *const f64,
    rates: *const f64,
    nrates: usize,
    rmask: [u8; MAX_DISTINCT], // rates occupied by these sites
    tensor_d: Vec<f64>,     // nrates x NS x NS (transposed) at the distal, proximal, pendant lengths
    tensor_p: Vec<f64>,
    tensor_q: Vec<f64>,
    dist_o: Vec<f64>,       // the edge's distal / proximal vectors at the masked sites
    prox_o: Vec<f64>,
    dist_e: Vec<i64>,
    prox_e: Vec<i64>,
    etot: f64,              // (tot(prox_e) + tot(dist_e)) + tot(pend_e = 0)
    nu: usize,              // distinct query vectors u[nu][NS], per-site index qc[m]
    u: Vec<f64>,
    qc: Vec<u8>,
    qp: Vec<f64>,           // nrates x nu x NS evolved query vectors at the pendant length
    pd: Vec<f64>,           // m x NS (gamma: m x R x NS)
    ll_core: f64,
    prox_bl: f64,
    dist_bl: f64,
    pend_bl: f64,
    n_calls: i64,
    assert_fail: bool,
    gmix: bool,             // gamma mixture: every site under all r rates
    r: usize,
}

impl Tt {
    /// The rows of query_masked (m x NS) as distinct vectors (byte-wise equality, as memcmp);
    /// None if more than MAX_DISTINCT.
    unsafe fn set_distinct_queries(&mut self, query_masked: *const f64) -> Option<()> {
        self.u = vec![0.0; MAX_DISTINCT * NS];
        self.qc = vec![0u8; self.m.max(1)];
        self.nu = 0;
        for s in 0..self.m {
            let q = sl(query_masked.add(s * NS), NS);
            let same = |u: usize| (0..NS).all(|i| self.u[u * NS + i].to_bits() == q[i].to_bits());
            let mut u = 0;
            while u < self.nu && !same(u) {
                u += 1;
            }
            if u == self.nu {
                if self.nu == MAX_DISTINCT {
                    return None;
                }
                self.u[u * NS..(u + 1) * NS].copy_from_slice(q);
                self.nu += 1;
            }
            self.qc[s] = u as u8;
        }
        self.qp = vec![0.0; self.nrates * self.nu.max(1) * NS];
        Some(())
    }

    unsafe fn multi_exp(&self, tensor: *mut f64, bl: f64) {
        xp_multi_exp(tensor, self.x, self.l, self.xit, self.rates, self.nrates as c_int, self.rmask.as_ptr(), bl);
    }

    /// The pendant length: evolve every distinct query vector.
    unsafe fn build_qp(&mut self, v: f64) {
        self.pend_bl = v;
        let p = self.tensor_q.as_mut_ptr();
        self.multi_exp(p, v);
        for r in 0..self.nrates {
            if self.rmask[r] == 0 {
                continue;
            }
            let pt = self.tensor_q.as_ptr().add(r * NS * NS);
            let q = self.qp.as_mut_ptr().add(r * self.nu * NS);
            for u in 0..self.nu {
                evolve_site(q.add(u * NS), self.u.as_ptr().add(u * NS), pt);
            }
        }
    }

    /// ll_core from pd and qp (mat_log_like3's sum).
    unsafe fn ll_pass(&mut self) {
        if self.gmix {
            self.g_ll_pass();
            return;
        }
        let mut ll_tot = 0.0;
        for s in 0..self.m {
            let pd = &self.pd[s * NS..(s + 1) * NS];
            let zo = (*self.cat.add(s) as usize * self.nu + self.qc[s] as usize) * NS;
            let zz = &self.qp[zo..zo + NS];
            let mut util = 0.0;
            for st in 0..NS {
                util += pd[st] * zz[st];
            }
            ll_tot += log(util);
        }
        self.ll_core = ll_tot;
    }

    /// pd and ll_core from the distal/proximal vectors and qp.
    unsafe fn dist_pass(&mut self) {
        if self.gmix {
            self.g_dist_pass();
            return;
        }
        let mut ll_tot = 0.0;
        let (mut dv, mut pv) = ([0.0; NS], [0.0; NS]);
        for s in 0..self.m {
            let off = *self.cat.add(s) as usize * NS * NS;
            evolve_site(dv.as_mut_ptr(), self.dist_o.as_ptr().add(s * NS), self.tensor_d.as_ptr().add(off));
            evolve_site(pv.as_mut_ptr(), self.prox_o.as_ptr().add(s * NS), self.tensor_p.as_ptr().add(off));
            let pd = self.pd.as_mut_ptr().add(s * NS);
            if SIMD20 {
                simd::pd20(pd, self.statd, pv.as_ptr(), dv.as_ptr());
            } else {
                for st in 0..NS {
                    *pd.add(st) = (*self.statd.add(st) * pv[st]) * dv[st];
                }
            }
            let zz = self.qp.as_ptr().add((*self.cat.add(s) as usize * self.nu + self.qc[s] as usize) * NS);
            let mut util = 0.0;
            for st in 0..NS {
                util += *pd.add(st) * *zz.add(st);
            }
            ll_tot += log(util);
        }
        self.ll_core = ll_tot;
    }

    /// Gamma: util[s] = sum over rates, then states, of pd * z (ten_log_like3_c), and
    /// ll = sum_s log(util[s]) - (float)m * log((float)R)   (double arithmetic; both exact ints).
    unsafe fn g_ll_pass(&mut self) {
        let (nu, r) = (self.nu, self.r);
        let mut ll_tot = 0.0;
        for s in 0..self.m {
            let mut util = 0.0;
            for k in 0..r {
                let pd = &self.pd[(s * r + k) * NS..(s * r + k + 1) * NS];
                let zo = (k * nu + self.qc[s] as usize) * NS;
                let zz = &self.qp[zo..zo + NS];
                for st in 0..NS {
                    util += pd[st] * zz[st];
                }
            }
            ll_tot += log(util);
        }
        ll_tot -= (self.m as f32 as f64) * log(r as f32 as f64);
        self.ll_core = ll_tot;
    }

    unsafe fn g_dist_pass(&mut self) {
        let r = self.r;
        let (mut dv, mut pv) = ([0.0; NS], [0.0; NS]);
        for s in 0..self.m {
            for k in 0..r {
                let o = (s * r + k) * NS;
                let off = k * NS * NS;
                evolve_site(dv.as_mut_ptr(), self.dist_o.as_ptr().add(o), self.tensor_d.as_ptr().add(off));
                evolve_site(pv.as_mut_ptr(), self.prox_o.as_ptr().add(o), self.tensor_p.as_ptr().add(off));
                for st in 0..NS {
                    self.pd[o + st] = (*self.statd.add(st) * pv[st]) * dv[st];
                }
            }
        }
        self.g_ll_pass();
    }

    fn log_like(&self) -> f64 {
        self.ll_core + (LN_2 * self.etot)
    }

    unsafe fn set_pend(&mut self, v: f64) {
        self.build_qp(v);
        self.ll_pass();
    }

    unsafe fn set_dist(&mut self, v: f64) {
        let cut = self.dist_bl + self.prox_bl;
        let prox = cut - v;
        if negative_or_nan(prox) {
            self.assert_fail = true; // alpha19: assert(prox_bl >= 0.)
        }
        self.dist_bl = v;
        let p = self.tensor_d.as_mut_ptr();
        self.multi_exp(p, v);
        self.prox_bl = prox;
        let p = self.tensor_p.as_mut_ptr();
        self.multi_exp(p, prox);
        self.dist_pass();
    }

    /// set_tt_edges: copy the masked sites of the edge's distal/proximal vectors (w doubles per
    /// site), then set the pendant, distal and proximal lengths in pplacer's order.
    unsafe fn set_edges(&mut self, w: usize, darr_a: *const f64, darr_e: *const i64, parr_a: *const f64,
                        parr_e: *const i64, reind: *const i32, cut_bl: f64, pendant: f64, distal: f64) {
        for i in 0..self.m {
            let src = *reind.add(i) as usize;
            ptr::copy_nonoverlapping(darr_a.add(src * w), self.dist_o.as_mut_ptr().add(i * w), w);
            ptr::copy_nonoverlapping(parr_a.add(src * w), self.prox_o.as_mut_ptr().add(i * w), w);
            self.dist_e[i] = *darr_e.add(src);
            self.prox_e[i] = *parr_e.add(src);
        }
        let m = self.m as c_int;
        self.etot = (xp_total_twoexp(self.prox_e.as_ptr(), m) + xp_total_twoexp(self.dist_e.as_ptr(), m)) + 0.0;
        self.build_qp(pendant); // Glv_edge.set_bl pend_edge first
        self.dist_bl = distal;
        let p = self.tensor_d.as_mut_ptr();
        self.multi_exp(p, distal);
        self.prox_bl = cut_bl - distal;
        let p = self.tensor_p.as_mut_ptr();
        self.multi_exp(p, cut_bl - distal);
        self.dist_pass();
        self.assert_fail = false;
    }
}

// --------------------------------------------------------------------- bit-level trace ---

static TRACE_ON: AtomicBool = AtomicBool::new(false);
static TRACE: Mutex<Option<BufWriter<File>>> = Mutex::new(None);

/// Debugging only: write every evaluated branch length and likelihood (C's "%a").
#[no_mangle]
pub unsafe extern "C" fn xp_trace_open(path: *const c_char) {
    let path = std::ffi::CStr::from_ptr(path).to_string_lossy().into_owned();
    if let Ok(f) = File::create(path) {
        *TRACE.lock().unwrap() = Some(BufWriter::new(f));
        TRACE_ON.store(true, Ordering::Relaxed);
    }
}

#[no_mangle]
pub extern "C" fn xp_trace_close() {
    TRACE_ON.store(false, Ordering::Relaxed);
    if let Some(mut f) = TRACE.lock().unwrap().take() {
        let _ = f.flush();
    }
}

fn trace(line: std::fmt::Arguments) {
    if TRACE_ON.load(Ordering::Relaxed) {
        if let Some(f) = TRACE.lock().unwrap().as_mut() {
            let _ = f.write_fmt(line);
        }
    }
}

/// glibc's printf("%a"): 0x1.<hex mantissa, trailing zeros dropped>p<exponent>.
fn hexfloat(v: f64) -> String {
    if v.is_nan() {
        return if v.is_sign_negative() { "-nan".into() } else { "nan".into() };
    }
    if v.is_infinite() {
        return if v < 0. { "-inf".into() } else { "inf".into() };
    }
    let sign = if v.is_sign_negative() { "-" } else { "" };
    let bits = v.to_bits();
    let (exp, man) = (((bits >> 52) & 0x7ff) as i64, bits & ((1 << 52) - 1));
    if exp == 0 && man == 0 {
        return format!("{sign}0x0p+0");
    }
    let (lead, e) = if exp == 0 { (0, -1022) } else { (1, exp - 1023) };
    let mut hex = format!("{man:013x}");
    while hex.ends_with('0') {
        hex.pop();
    }
    let dot = if hex.is_empty() { String::new() } else { format!(".{hex}") };
    format!("{sign}0x{lead}{dot}p{}{e}", if e >= 0 { "+" } else { "" })
}

// ------------------------------------------------------------------------ optimisation ---

/// The objective seen by the minimiser: which = 0 pendant, 1 distal.
struct OptCtx {
    t: *mut Tt,
    which: i32,
    in_gsl: bool,
}

extern "C" fn opt_fun(v: f64, p: *mut c_void) -> f64 {
    unsafe {
        let c = &mut *(p as *mut OptCtx);
        let t = c.t;
        (*t).n_calls += 1;
        trace(format_args!("{} {}\n", if c.which == 0 { 'P' } else { 'D' }, hexfloat(v)));
        if c.which == 1 && negative_or_nan(((*t).dist_bl + (*t).prox_bl) - v) {
            // alpha19 set_dist_bl: assert(prox_bl >= 0.) raises BEFORE any branch length changes.
            // Inside a GSL callback ocamlgsl returns NaN on exception: GSL fails with
            // GSL_EBADFUNC, and safe_ml_optimize_location retries from the default lengths.
            // Outside GSL (the bisection start finder) the assertion would abort pplacer.
            if !c.in_gsl {
                (*t).assert_fail = true;
            }
            return f64::NAN;
        }
        if c.which == 0 {
            (*t).set_pend(v);
        } else {
            (*t).set_dist(v);
        }
        trace(format_args!("L {}\n", hexfloat((*t).ll_core)));
        -(*t).log_like()
    }
}

/// Minimization.brent with bisection_start_finder (max_iters = 100), on GSL's Brent.
/// The context is only used through the raw pointer GSL also calls opt_fun with.
unsafe fn brent(c: *mut OptCtx, raw_start: f64, left: f64, right: f64, tol: f64, res: &mut f64) -> c_int {
    const MAX_ITERS: i32 = 100;
    let cp = c as *mut c_void;
    if left >= raw_start || raw_start >= right {
        return XP_INVALIDSTART;
    }
    // OCaml evaluates the tuple (f left, f right) right-to-left
    let righty = opt_fun(right, cp);
    let lefty = opt_fun(left, cp);
    let smaller = if lefty < righty { left } else { right };
    let miny = if lefty <= righty { lefty } else { righty };
    let mut prev = raw_start;
    let mut it = 1;
    let start = loop {
        let pv = opt_fun(prev, cp);
        if pv < miny {
            break prev;
        } else if (prev - smaller).abs() < tol {
            *res = if pv < miny { prev } else { smaller };
            return XP_OK;
        } else if it > MAX_ITERS {
            return XP_FINDSTART;
        }
        prev = (prev + smaller) / 2.;
        it += 1;
    };
    let mut f = GslFunction { function: opt_fun, params: cp };
    let s = ffi::gsl_min_fminimizer_alloc(ffi::gsl_min_fminimizer_brent);
    (*c).in_gsl = true;
    let finish = |status: c_int| {
        ffi::gsl_min_fminimizer_free(s);
        status
    };
    if ffi::gsl_min_fminimizer_set(s, &mut f, start, left, right) != 0 {
        (*c).in_gsl = false;
        return finish(XP_GSL);
    }
    for step in 1.. {
        if step > MAX_ITERS {
            (*c).in_gsl = false;
            return finish(XP_EXCEEDED);
        }
        if ffi::gsl_min_fminimizer_iterate(s) != 0 {
            (*c).in_gsl = false;
            return finish(XP_GSL);
        }
        let (lo, hi) = (ffi::gsl_min_fminimizer_x_lower(s), ffi::gsl_min_fminimizer_x_upper(s));
        if hi - lo < tol {
            (*c).in_gsl = false;
            *res = ffi::gsl_min_fminimizer_x_minimum(s);
            return finish(XP_OK);
        }
    }
    unreachable!()
}

/// Three_tax.optimize: alternate pendant then distal optimisation; stop when both move <= tol.
#[no_mangle]
pub unsafe extern "C" fn xp_tt_optimize(t: *mut Tt, tol: f64, max_pend: f64, max_iter: c_int) -> c_int {
    let mut cp = OptCtx { t, which: 0, in_gsl: false };
    let mut cd = OptCtx { t, which: 1, in_gsl: false };
    let (mut prev_pend, mut prev_dist) = ((*t).pend_bl, (*t).dist_bl);
    (*t).n_calls = 0;
    trace(format_args!("O\n"));
    for step in 1.. {
        if step > max_iter {
            return XP_EXCEEDED;
        }
        let mut cur_pend = 0.0;
        if 1e-8 == max_pend {
            cur_pend = (*t).pend_bl; // optimize_pend_bl: if min = max then start
        } else {
            let st = brent(&mut cp, (*t).pend_bl, 1e-8, max_pend, tol, &mut cur_pend);
            if st != XP_OK {
                return st;
            }
        }
        if (*t).assert_fail {
            return XP_ASSERT;
        }
        let cut = (*t).dist_bl + (*t).prox_bl;
        let mut cur_dist = 0.0;
        if 0. == cut {
            cur_dist = (*t).dist_bl;
        } else {
            let st = brent(&mut cd, (*t).dist_bl, 0., cut, tol, &mut cur_dist);
            if st != XP_OK {
                return st;
            }
        }
        if (*t).assert_fail {
            return XP_ASSERT;
        }
        if (prev_pend - cur_pend).abs() > tol || (prev_dist - cur_dist).abs() > tol {
            prev_pend = cur_pend;
            prev_dist = cur_dist;
        } else {
            return XP_OK;
        }
    }
    unreachable!()
}

// ---------------------------------------------------------------------- construction ---

unsafe fn new_tt(m: usize, nrates: usize, statd: *const f64, x: *const f64, l: *const f64, xit: *const f64,
                 rates: *const f64, w: usize) -> Tt {
    let site = (m * w).max(1);
    Tt {
        m, cat: ptr::null(), statd, x, l, xit, rates, nrates, rmask: [0; MAX_DISTINCT],
        tensor_d: vec![0.0; nrates * NS * NS],
        tensor_p: vec![0.0; nrates * NS * NS],
        tensor_q: vec![0.0; nrates * NS * NS],
        dist_o: vec![0.0; site], prox_o: vec![0.0; site], pd: vec![0.0; site],
        dist_e: vec![0; m.max(1)], prox_e: vec![0; m.max(1)],
        etot: 0.0, nu: 0, u: Vec::new(), qc: Vec::new(), qp: Vec::new(), ll_core: 0.0,
        prox_bl: 0.0, dist_bl: 0.0, pend_bl: 0.0, n_calls: 0, assert_fail: false, gmix: false, r: 0,
    }
}

/// A three-taxon tree for one query with m masked sites (CAT). NULL if the query has more than 64
/// distinct site vectors (impossible for protein input).
#[no_mangle]
pub unsafe extern "C" fn xp_tt_new(m: c_int, cat: *const i32, statd: *const f64, x: *const f64, l: *const f64,
                                   xit: *const f64, rates: *const f64, nrates: c_int,
                                   query_masked: *const f64) -> *mut Tt {
    let m = m as usize;
    let mut t = new_tt(m, nrates as usize, statd, x, l, xit, rates, NS);
    t.cat = cat;
    for s in 0..m {
        t.rmask[*cat.add(s) as usize] = 1;
    }
    match t.set_distinct_queries(query_masked) {
        Some(()) => Box::into_raw(Box::new(t)),
        None => ptr::null_mut(),
    }
}

/// A three-taxon tree for the gamma mixture: m masked sites, all r rates at every site.
#[no_mangle]
pub unsafe extern "C" fn xp_gtt_new(m: c_int, rr: c_int, statd: *const f64, x: *const f64, l: *const f64,
                                    xit: *const f64, rates: *const f64, query_masked: *const f64) -> *mut Tt {
    let (m, r) = (m as usize, rr as usize);
    if r > MAX_DISTINCT {
        return ptr::null_mut();
    }
    let mut t = new_tt(m, r, statd, x, l, xit, rates, r * NS);
    t.gmix = true;
    t.r = r;
    t.rmask[..r].fill(1);
    match t.set_distinct_queries(query_masked) {   // the query vector is the same at every rate
        Some(()) => Box::into_raw(Box::new(t)),
        None => ptr::null_mut(),
    }
}

#[no_mangle]
pub unsafe extern "C" fn xp_tt_free(t: *mut Tt) {
    if !t.is_null() {
        drop(Box::from_raw(t));
    }
}

/// Set the tree on an edge (CAT): darr/parr are the edge's full-site distal/proximal vectors.
#[no_mangle]
pub unsafe extern "C" fn xp_tt_set(t: *mut Tt, darr_a: *const f64, darr_e: *const i64, parr_a: *const f64,
                                   parr_e: *const i64, reind: *const i32, cut_bl: f64, pendant: f64, distal: f64) {
    (*t).set_edges(NS, darr_a, darr_e, parr_a, parr_e, reind, cut_bl, pendant, distal);
}

/// As xp_tt_set for the gamma mixture: full-site vectors are nsites x R x NS.
#[no_mangle]
pub unsafe extern "C" fn xp_gtt_set(t: *mut Tt, darr_a: *const f64, darr_e: *const i64, parr_a: *const f64,
                                    parr_e: *const i64, reind: *const i32, cut_bl: f64, pendant: f64, distal: f64) {
    let w = (*t).r * NS;
    (*t).set_edges(w, darr_a, darr_e, parr_a, parr_e, reind, cut_bl, pendant, distal);
}

/// out = (log likelihood, pendant length, distal length, number of evaluations).
#[no_mangle]
pub unsafe extern "C" fn xp_tt_results(t: *const Tt, out: *mut f64) {
    *out = (*t).log_like();
    *out.add(1) = (*t).pend_bl;
    *out.add(2) = (*t).dist_bl;
    *out.add(3) = (*t).n_calls as f64;
}
