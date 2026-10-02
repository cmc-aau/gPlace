//! The prefilter: pplacer's h_ranking score of every query against every edge.
//!
//! pplacer scores query q at edge v as
//!     sum_{s informative in q} log( sum_st qev_q[s][st] * snode[s][st] )  +  log2 * (0 + sum se[s])
//! where qev_q[s] (the query character evolved by the starting pendant length) depends only on the
//! site's rate category and the character. So for each edge we first build
//!     W[s][c] = log( sum_st qevtab[cat[s]][c][st] * snode[s][st] )     (the same operands, in order)
//! and a query's score is a sum of W lookups in site order: bit-identical to pplacer.

use crate::ffi::log;
use crate::simd::{self, EB};
use crate::{sl, sl_mut, NS};
use std::f64::consts::LN_2;
use std::os::raw::{c_int, c_long};

/// One edge against many queries; codes[q][s] is a character code, or 255 for a gap.
/// out[q * stride] receives the score. w: scratch of nsites x ncode. (Reference for the self-tests.)
#[no_mangle]
pub unsafe extern "C" fn xp_edge_scores(out: *mut f64, stride: c_long, snode: *const f64, se: *const i64,
                                        qevtab: *const f64, cat: *const i32, nsites: c_int, ncode: c_int,
                                        codes: *const u8, nq: c_int, w: *mut f64) {
    let (ns, nc, nq) = (nsites as usize, ncode as usize, nq as usize);
    let (snode, se, cat, codes) = (sl(snode, ns * NS), sl(se, ns), sl(cat, ns), sl(codes, nq * ns));
    let w = sl_mut(w, ns * nc);
    for s in 0..ns {
        let sn = &snode[s * NS..(s + 1) * NS];
        let qt = qevtab.add(cat[s] as usize * nc * NS);
        for c in 0..nc {
            let qv = sl(qt.add(c * NS), NS);
            let mut util = 0.;
            for st in 0..NS {
                util += qv[st] * sn[st];
            }
            w[s * nc + c] = log(util);
        }
    }
    for q in 0..nq {
        let cq = &codes[q * ns..(q + 1) * ns];
        let (mut ll, mut et) = (0.0f64, 0i64);
        for s in 0..ns {
            let c = cq[s];
            if c != 255 {
                ll += w[s * nc + c as usize];
                et += se[s];
            }
        }
        *out.add(q * stride as usize) = ll + (LN_2 * (0. + et as f64));
    }
}

/// The scalar form of the 8-edge prefilter: per lane and query, the W values in site order.
unsafe fn es8_scalar(out: *mut f64, ldo: usize, col0: usize, ne: usize, w: &[f64], se8: &[i64],
                     nsites: usize, ncode: usize, codes: &[u8], nq: usize) {
    let nc = ncode + 1;
    for lane in 0..ne {
        for q in 0..nq {
            let cq = &codes[q * nsites..(q + 1) * nsites];
            let (mut ll, mut et) = (0.0f64, 0i64);
            for s in 0..nsites {
                let c = cq[s] as usize;
                if c != ncode {
                    ll += w[(s * nc + c) * EB + lane];
                    et += se8[lane * nsites + s];
                }
            }
            *out.add(q * ldo + col0 + lane) = ll + (LN_2 * (0. + et as f64));
        }
    }
}

/// The prefilter for 8 edges at once (lanes = edges): for each lane and query exactly
/// xp_edge_scores' value. codes[q][s] in 0..ncode-1, or ncode for a gap ('-' or '?');
/// out[q * ldo + col0 + lane] receives the score for lane < ne. Returns 0, or -1 on bad arguments.
///
/// With SIMD, gap sites add W = +0.0 instead of being skipped: ll can never be -0.0 (log never
/// returns -0.0 and x + (-x) = +0.0 under round-to-nearest), so ll + 0.0 == ll bit for bit.
/// Exponents are summed in int32 lanes within one 128-site block and added into int64 totals
/// (per-site exponents reach ~2^18 on the full GTDB bacterial tree); integer addition is exact in
/// any order. Exponents beyond 2^23 use the scalar form.
#[no_mangle]
pub unsafe extern "C" fn xp_edge_scores8(out: *mut f64, ldo: c_long, col0: c_long, ne: c_int, sn8: *const f64,
                                         se8: *const i64, qevtab: *const f64, cat: *const i32, nsites: c_int,
                                         ncode: c_int, codes: *const u8, nq: c_int) -> c_int {
    if !(1..=EB as c_int).contains(&ne) || !(1..=254).contains(&ncode) {
        return -1;
    }
    let (ne, ns, ncode, nq) = (ne as usize, nsites as usize, ncode as usize, nq as usize);
    let nc = ncode + 1;
    let (sn8, se8, cat, codes) = (sl(sn8, EB * ns * NS), sl(se8, EB * ns), sl(cat, ns), sl(codes, nq * ns));
    let mut w = vec![0.0f64; ns * nc * EB];
    let mut se = vec![0i32; ns * 2 * EB]; // per site: 8 zeros (gap), then the 8 lanes' exponents
    let mut small = true;
    for s in 0..ns {
        let qt = qevtab.add(cat[s] as usize * ncode * NS);
        let ws = &mut w[s * nc * EB..(s + 1) * nc * EB];
        for c in 0..ncode {
            let qv = sl(qt.add(c * NS), NS);
            for lane in 0..ne {
                let sn = &sn8[(lane * ns + s) * NS..(lane * ns + s + 1) * NS];
                let mut util = 0.;
                for st in 0..NS {
                    util += qv[st] * sn[st];
                }
                ws[c * EB + lane] = log(util);
            } // lanes >= ne and the gap row stay 0.0
        }
        for lane in 0..ne {
            let v = se8[lane * ns + s];
            if !(-(1 << 23)..=(1 << 23)).contains(&v) {
                small = false;
            }
            se[s * 2 * EB + EB + lane] = v as i32;
        }
    }
    if !small || !(simd::AVX2 || simd::NEON) {
        es8_scalar(out, ldo as usize, col0 as usize, ne, &w, se8, ns, ncode, codes, nq);
        return 0;
    }
    let mut acc_d = vec![0.0f64; nq * EB]; // +0.0
    let mut acc_e = vec![0i64; nq * EB];
    simd::scores8_blocks(acc_d.as_mut_ptr(), acc_e.as_mut_ptr(), w.as_ptr(), se.as_ptr(), codes.as_ptr(),
                         ns, nc, ncode, nq);
    for q in 0..nq {
        for lane in 0..ne {
            *out.add(q * ldo as usize + col0 as usize + lane) =
                acc_d[q * EB + lane] + (LN_2 * (0. + acc_e[q * EB + lane] as f64));
        }
    }
    0
}

/// Gmix prefilter (Glv.masked_logdot of the evolved query and the supernode), 8 edges at once:
///   W[s][c] = log(util), util = sum over r ascending, st ascending of qev[r][c][st] * sn[s][r][st]
///   ll = +0.0 + sum_{s informative} W[s][code]
///   ll = ll - (float)nused[q] * log((float)R)            (ten_masked_logdot_c's subtraction)
///   score = ll + log2 * (0 + sum_{s informative} se[s])
/// sn8: 8 x nsites x R x NS; se8: 8 x nsites; qevtab: R x ncode x NS; codes as xp_edge_scores8.
#[no_mangle]
pub unsafe extern "C" fn xp_g_edge_scores8(out: *mut f64, ldo: c_long, col0: c_long, ne: c_int, sn8: *const f64,
                                           se8: *const i64, qevtab: *const f64, rr: c_int, nsites: c_int,
                                           ncode: c_int, codes: *const u8, nused: *const i32, nq: c_int) -> c_int {
    if !(1..=EB as c_int).contains(&ne) || !(1..=254).contains(&ncode) {
        return -1;
    }
    let (ne, ns, ncode, nq, r) = (ne as usize, nsites as usize, ncode as usize, nq as usize, rr as usize);
    let nc = ncode + 1;
    let (sn8, se8, qevtab) = (sl(sn8, EB * ns * r * NS), sl(se8, EB * ns), sl(qevtab, r * ncode * NS));
    let (codes, nused) = (sl(codes, nq * ns), sl(nused, nq));
    let mut w = vec![0.0f64; ns * nc * EB];
    for s in 0..ns {
        for c in 0..ncode {
            for lane in 0..ne {
                let sn = &sn8[(lane * ns + s) * r * NS..(lane * ns + s + 1) * r * NS];
                let mut util = 0.0;
                for k in 0..r {
                    let qv = &qevtab[(k * ncode + c) * NS..(k * ncode + c + 1) * NS];
                    let sv = &sn[k * NS..(k + 1) * NS];
                    for st in 0..NS {
                        util += qv[st] * sv[st];
                    }
                }
                w[(s * nc + c) * EB + lane] = log(util);
            }
        }
    }
    // log((float) R): the float argument is promoted to double; (float) nused is exact
    let log_r = log(r as f32 as f64);
    for q in 0..nq {
        let cq = &codes[q * ns..(q + 1) * ns];
        for lane in 0..ne {
            let (mut ll, mut et) = (0.0f64, 0i64);
            for s in 0..ns {
                let c = cq[s] as usize;
                if c != ncode {
                    ll += w[(s * nc + c) * EB + lane];
                    et += se8[lane * ns + s];
                }
            }
            ll -= (nused[q] as f32 as f64) * log_r;
            *out.add(q * ldo as usize + col0 as usize + lane) = ll + (LN_2 * (0. + et as f64));
        }
    }
    0
}
