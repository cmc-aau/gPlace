//! SIMD kernels: AVX2 (x86-64, 4 doubles per register) or NEON (ARM64, 2 per register).
//!
//! Each lane holds a DIFFERENT result (a different output state, or a different edge), computed
//! with exactly the operations of the scalar loop, in the same order: a multiply, then an add,
//! never fused; sums start from +0.0 and run in the scalar order. So every kernel here gives the
//! same bits as its scalar form (`glv`, `model`, `prefilter`), which the self-tests check.
//!
//! The 20-state kernels (`dediag20`, `evolve_site20`, `pd20`) are used only when `SIMD20`; the
//! prefilter kernel (`scores8_blocks`) for any state count. In a scalar build the functions exist
//! but are never called.

/// Kernels of this build.
pub const AVX2: bool = cfg!(all(target_arch = "x86_64", target_feature = "avx2", not(xp_force_scalar)));
pub const NEON: bool = cfg!(all(target_arch = "aarch64", target_feature = "neon", not(xp_force_scalar)));

/// Prefilter lanes: 8 edges per call.
pub const EB: usize = 8;
/// Sites per block in the prefilter: int32 exponent sums over one block cannot overflow when
/// every |exponent| <= 2^23 (128 * 2^23 = 2^30). Used by the AVX2 and NEON kernels only.
#[allow(dead_code)]
pub const BS: usize = 128;

pub use imp::*;

// ------------------------------------------------------------------------------- AVX2 ---
#[cfg(all(target_arch = "x86_64", target_feature = "avx2", not(xp_force_scalar)))]
mod imp {
    use super::{BS, EB};
    use std::arch::x86_64::*;

    /// dediagonalize, transposed: dst_t[j][i] = sum_k a[k][i] * uit[j][k], with
    /// a[k][i] = lambda[k] * u[i][k] (the scalar `d += lambda[k] * u[i][k] * uit[j][k]`).
    pub unsafe fn dediag20(dst_t: *mut f64, a: *const f64, uit: *const f64) {
        for j in 0..20 {
            let uj = uit.add(j * 20);
            let mut acc = [_mm256_setzero_pd(); 5];
            for k in 0..20 {
                let w = _mm256_set1_pd(*uj.add(k));
                let ak = a.add(k * 20);
                for (r, x) in acc.iter_mut().enumerate() {
                    *x = _mm256_add_pd(*x, _mm256_mul_pd(_mm256_loadu_pd(ak.add(4 * r)), w));
                }
            }
            for (r, x) in acc.iter().enumerate() {
                _mm256_storeu_pd(dst_t.add(j * 20 + 4 * r), *x);
            }
        }
    }

    /// y[i] = 0.0 + sum_j x[j] * pt[j][i]   (gslcblas dgemv; lanes run over i)
    #[inline(always)]
    pub unsafe fn evolve_site20(yv: *mut f64, xv: *const f64, pt: *const f64) {
        let z = _mm256_setzero_pd();
        let mut acc = [z; 5];
        for j in 0..20 {
            let xj = _mm256_set1_pd(*xv.add(j));
            let pj = pt.add(j * 20);
            for (r, a) in acc.iter_mut().enumerate() {
                *a = _mm256_add_pd(*a, _mm256_mul_pd(xj, _mm256_loadu_pd(pj.add(4 * r))));
            }
        }
        for (r, a) in acc.iter().enumerate() {
            _mm256_storeu_pd(yv.add(4 * r), _mm256_add_pd(z, *a));
        }
    }

    /// pd[st] = (statd[st] * pv[st]) * dv[st]
    #[inline(always)]
    pub unsafe fn pd20(pd: *mut f64, statd: *const f64, pv: *const f64, dv: *const f64) {
        for r in 0..5 {
            let o = 4 * r;
            let s = _mm256_loadu_pd(statd.add(o));
            _mm256_storeu_pd(pd.add(o), _mm256_mul_pd(_mm256_mul_pd(s, _mm256_loadu_pd(pv.add(o))),
                                                      _mm256_loadu_pd(dv.add(o))));
        }
    }

    /// Add 8 int32 lane sums into 8 int64 totals (exact).
    #[inline(always)]
    unsafe fn add64(e_tot: *mut i64, e: __m256i) {
        let lo = _mm256_cvtepi32_epi64(_mm256_castsi256_si128(e));
        let hi = _mm256_cvtepi32_epi64(_mm256_extracti128_si256::<1>(e));
        let p = e_tot as *mut __m256i;
        _mm256_storeu_si256(p, _mm256_add_epi64(_mm256_loadu_si256(p), lo));
        _mm256_storeu_si256(p.add(1), _mm256_add_epi64(_mm256_loadu_si256(p.add(1)), hi));
    }

    /// The prefilter's accumulation (see prefilter::xp_edge_scores8): for every query q and lane,
    /// acc_d[q][lane] += w[s][code][lane] over sites s in order, and acc_e[q][lane] += the site's
    /// exponent if the code is not a gap. Blocks of BS sites; 4 queries interleaved.
    pub unsafe fn scores8_blocks(acc_d: *mut f64, acc_e: *mut i64, w: *const f64, se: *const i32,
                                 codes: *const u8, nsites: usize, nc: usize, gap: usize, nq: usize) {
        let mut s0 = 0;
        while s0 < nsites {
            let s1 = (s0 + BS).min(nsites);
            let mut q = 0;
            while q + 4 <= nq {
                let a = acc_d.add(q * EB);
                let mut ad = [[_mm256_setzero_pd(); 2]; 4];
                for (k, x) in ad.iter_mut().enumerate() {
                    x[0] = _mm256_loadu_pd(a.add(k * EB));
                    x[1] = _mm256_loadu_pd(a.add(k * EB + 4));
                }
                let mut e = [_mm256_setzero_si256(); 4];           // block partial sums
                for s in s0..s1 {
                    let ws = w.add(s * nc * EB);
                    let ses = se.add(s * 2 * EB);
                    for k in 0..4 {
                        let c = *codes.add((q + k) * nsites + s) as usize;
                        let wk = ws.add(c * EB);
                        ad[k][0] = _mm256_add_pd(ad[k][0], _mm256_loadu_pd(wk));
                        ad[k][1] = _mm256_add_pd(ad[k][1], _mm256_loadu_pd(wk.add(4)));
                        let x = ses.add(((c != gap) as usize) * EB) as *const __m256i;
                        e[k] = _mm256_add_epi32(e[k], _mm256_loadu_si256(x));
                    }
                }
                for k in 0..4 {
                    _mm256_storeu_pd(a.add(k * EB), ad[k][0]);
                    _mm256_storeu_pd(a.add(k * EB + 4), ad[k][1]);
                    add64(acc_e.add((q + k) * EB), e[k]);
                }
                q += 4;
            }
            while q < nq {
                let a = acc_d.add(q * EB);
                let (mut lo, mut hi) = (_mm256_loadu_pd(a), _mm256_loadu_pd(a.add(4)));
                let mut e = _mm256_setzero_si256();
                for s in s0..s1 {
                    let c = *codes.add(q * nsites + s) as usize;
                    let wk = w.add((s * nc + c) * EB);
                    lo = _mm256_add_pd(lo, _mm256_loadu_pd(wk));
                    hi = _mm256_add_pd(hi, _mm256_loadu_pd(wk.add(4)));
                    let x = se.add(s * 2 * EB + ((c != gap) as usize) * EB) as *const __m256i;
                    e = _mm256_add_epi32(e, _mm256_loadu_si256(x));
                }
                _mm256_storeu_pd(a, lo);
                _mm256_storeu_pd(a.add(4), hi);
                add64(acc_e.add(q * EB), e);
                q += 1;
            }
            s0 = s1;
        }
    }
}

// ------------------------------------------------------------------------------- NEON ---
#[cfg(all(target_arch = "aarch64", target_feature = "neon", not(xp_force_scalar)))]
mod imp {
    use super::{BS, EB};
    use std::arch::aarch64::*;

    /// As the AVX2 dediag20, 2 lanes per register.
    pub unsafe fn dediag20(dst_t: *mut f64, a: *const f64, uit: *const f64) {
        for j in 0..20 {
            let uj = uit.add(j * 20);
            let mut acc = [vdupq_n_f64(0.0); 10];
            for k in 0..20 {
                let w = vdupq_n_f64(*uj.add(k));
                let ak = a.add(k * 20);
                for (r, x) in acc.iter_mut().enumerate() {
                    *x = vaddq_f64(*x, vmulq_f64(vld1q_f64(ak.add(2 * r)), w));
                }
            }
            for (r, x) in acc.iter().enumerate() {
                vst1q_f64(dst_t.add(j * 20 + 2 * r), *x);
            }
        }
    }

    #[inline(always)]
    pub unsafe fn evolve_site20(yv: *mut f64, xv: *const f64, pt: *const f64) {
        let z = vdupq_n_f64(0.0);
        let mut acc = [z; 10];
        for j in 0..20 {
            let xj = vdupq_n_f64(*xv.add(j));
            let pj = pt.add(j * 20);
            for (r, a) in acc.iter_mut().enumerate() {
                *a = vaddq_f64(*a, vmulq_f64(xj, vld1q_f64(pj.add(2 * r))));
            }
        }
        for (r, a) in acc.iter().enumerate() {
            vst1q_f64(yv.add(2 * r), vaddq_f64(z, *a));
        }
    }

    #[inline(always)]
    pub unsafe fn pd20(pd: *mut f64, statd: *const f64, pv: *const f64, dv: *const f64) {
        for r in 0..10 {
            let o = 2 * r;
            let s = vld1q_f64(statd.add(o));
            vst1q_f64(pd.add(o), vmulq_f64(vmulq_f64(s, vld1q_f64(pv.add(o))), vld1q_f64(dv.add(o))));
        }
    }

    /// Add 8 int32 lane sums (lo: lanes 0-3, hi: 4-7) into 8 int64 totals (exact).
    #[inline(always)]
    unsafe fn add64(e_tot: *mut i64, lo: int32x4_t, hi: int32x4_t) {
        vst1q_s64(e_tot, vaddq_s64(vld1q_s64(e_tot), vmovl_s32(vget_low_s32(lo))));
        vst1q_s64(e_tot.add(2), vaddq_s64(vld1q_s64(e_tot.add(2)), vmovl_high_s32(lo)));
        vst1q_s64(e_tot.add(4), vaddq_s64(vld1q_s64(e_tot.add(4)), vmovl_s32(vget_low_s32(hi))));
        vst1q_s64(e_tot.add(6), vaddq_s64(vld1q_s64(e_tot.add(6)), vmovl_high_s32(hi)));
    }

    /// As the AVX2 scores8_blocks: 8 lanes = 4 registers per query, 4 queries interleaved.
    pub unsafe fn scores8_blocks(acc_d: *mut f64, acc_e: *mut i64, w: *const f64, se: *const i32,
                                 codes: *const u8, nsites: usize, nc: usize, gap: usize, nq: usize) {
        let mut s0 = 0;
        while s0 < nsites {
            let s1 = (s0 + BS).min(nsites);
            let mut q = 0;
            while q + 4 <= nq {
                let a = acc_d.add(q * EB);
                let mut ad = [[vdupq_n_f64(0.0); 4]; 4];
                for (k, x) in ad.iter_mut().enumerate() {
                    for (r, v) in x.iter_mut().enumerate() {
                        *v = vld1q_f64(a.add(k * EB + 2 * r));
                    }
                }
                let mut e = [[vdupq_n_s32(0); 2]; 4];
                for s in s0..s1 {
                    let ws = w.add(s * nc * EB);
                    let ses = se.add(s * 2 * EB);
                    for k in 0..4 {
                        let c = *codes.add((q + k) * nsites + s) as usize;
                        let wk = ws.add(c * EB);
                        for (r, v) in ad[k].iter_mut().enumerate() {
                            *v = vaddq_f64(*v, vld1q_f64(wk.add(2 * r)));
                        }
                        let x = ses.add(((c != gap) as usize) * EB);
                        e[k][0] = vaddq_s32(e[k][0], vld1q_s32(x));
                        e[k][1] = vaddq_s32(e[k][1], vld1q_s32(x.add(4)));
                    }
                }
                for k in 0..4 {
                    for (r, v) in ad[k].iter().enumerate() {
                        vst1q_f64(a.add(k * EB + 2 * r), *v);
                    }
                    add64(acc_e.add((q + k) * EB), e[k][0], e[k][1]);
                }
                q += 4;
            }
            while q < nq {
                let a = acc_d.add(q * EB);
                let mut ad = [vdupq_n_f64(0.0); 4];
                for (r, v) in ad.iter_mut().enumerate() {
                    *v = vld1q_f64(a.add(2 * r));
                }
                let (mut el, mut eh) = (vdupq_n_s32(0), vdupq_n_s32(0));
                for s in s0..s1 {
                    let c = *codes.add(q * nsites + s) as usize;
                    let wk = w.add((s * nc + c) * EB);
                    for (r, v) in ad.iter_mut().enumerate() {
                        *v = vaddq_f64(*v, vld1q_f64(wk.add(2 * r)));
                    }
                    let x = se.add(s * 2 * EB + ((c != gap) as usize) * EB);
                    el = vaddq_s32(el, vld1q_s32(x));
                    eh = vaddq_s32(eh, vld1q_s32(x.add(4)));
                }
                for (r, v) in ad.iter().enumerate() {
                    vst1q_f64(a.add(2 * r), *v);
                }
                add64(acc_e.add(q * EB), el, eh);
                q += 1;
            }
            s0 = s1;
        }
    }
}

// ---------------------------------------------------------------- scalar build (stubs) ---
#[cfg(not(any(all(target_arch = "x86_64", target_feature = "avx2", not(xp_force_scalar)),
              all(target_arch = "aarch64", target_feature = "neon", not(xp_force_scalar)))))]
mod imp {
    pub unsafe fn dediag20(_: *mut f64, _: *const f64, _: *const f64) {
        unreachable!("scalar engine")
    }
    pub unsafe fn evolve_site20(_: *mut f64, _: *const f64, _: *const f64) {
        unreachable!("scalar engine")
    }
    pub unsafe fn pd20(_: *mut f64, _: *const f64, _: *const f64, _: *const f64) {
        unreachable!("scalar engine")
    }
    pub unsafe fn scores8_blocks(_: *mut f64, _: *mut i64, _: *const f64, _: *const i32, _: *const u8,
                                 _: usize, _: usize, _: usize, _: usize) {
        unreachable!("scalar engine")
    }
}
