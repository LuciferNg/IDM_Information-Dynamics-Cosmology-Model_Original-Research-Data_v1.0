"""
joint_validation_v2.py — Optimized CMB + BAO + SNe + H(z) 聯合驗證
==================================================================
Optimizations:
  1. Precompute E(z) integral on fine grid → cubic spline interpolation
  2. Coarse → fine grid strategy
  3. Vectorized χ² where possible
"""

import numpy as np
from scipy import integrate, interpolate
import time

# ============================================================
# Constants
# ============================================================
C = 299792.458
Z_REC = 1089.9
R_S = 147.09

PLANCK_R = 1.7502
PLANCK_R_ERR = 0.0046
PLANCK_LA = 301.808
PLANCK_LA_ERR = 0.090
PLANCK_R_LA_CORR = 0.47

# ============================================================
# Interpolated integrator — FAST
# ============================================================

class FastIntegrator:
    """Precomputes ∫₀ᶻ dz/E(z) on a fine grid, then interpolates."""
    def __init__(self, model, z_max=3.0, n=3000):
        self.model = model
        self.z_grid = np.linspace(0, z_max, n)
        self.z_max = z_max
        self.z_min = z_max  # will reduce as needed
        
        # Vectorized integral
        E_vals = model.E(self.z_grid)
        dz = self.z_grid[1] - self.z_grid[0]
        integrand = 1.0 / E_vals
        cum_int = np.zeros(n)
        cum_int[1:] = np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * dz)
        self.cum_int = cum_int
        
        # Build interpolation ONCE
        self._interp = interpolate.interp1d(self.z_grid, cum_int, kind='cubic')
        
    def integral(self, z):
        """∫₀ᶻ dz'/E(z') — vectorized, interpolated. z can be array."""
        z = np.asarray(z, dtype=float)
        single = z.ndim == 0
        z = np.atleast_1d(z)
        
        mask_in = z <= self.z_max
        result = np.zeros_like(z)
        
        # In-grid: use prebuilt interpolation
        if np.any(mask_in):
            result[mask_in] = self._interp(np.minimum(z[mask_in], self.z_max))
        
        # Out-grid: integrate the tail
        if np.any(~mask_in):
            f = lambda x: 1.0 / self.model.E(x)
            tail_base = float(self._interp(self.z_max))  # cache base
            for i in np.where(~mask_in)[0]:
                tail, _ = integrate.quad(f, self.z_max, z[i], limit=200, epsabs=1e-10)
                result[i] = tail_base + tail
        
        return float(result[0]) if single else result


# ============================================================
# IDM Model (same as before, with minor optimizations)
# ============================================================

class IDMModel:
    def __init__(self, H0=67.4, Om=0.315, eps=0.0, zc=0.6):
        self.H0 = H0
        self.Om = Om
        self.Ol = 1.0 - Om
        self.eps = eps
        self.zc = zc
    
    def f_idm(self, z):
        x = z / self.zc
        return 1.0 - self.eps * x * np.exp(-x)
    
    def E2(self, z):
        return self.Om * (1+z)**3 + self.Ol * self.f_idm(z)
    
    def E(self, z):
        return np.sqrt(self.E2(z))
    
    def H(self, z):
        return self.H0 * self.E(z)
    
    def comoving_distance(self, z, integrator=None):
        if integrator is not None:
            I = integrator.integral(z)
        else:
            f = lambda x: 1.0 / self.E(np.array([x]))[0]
            I, _ = integrate.quad(f, 0, z, limit=500)
        return C / self.H0 * I
    
    def lum_distance(self, z, integrator=None):
        return (1+z) * self.comoving_distance(z, integrator)
    
    def distance_modulus(self, z, integrator=None):
        dl = self.lum_distance(z, integrator)
        return 5.0 * np.log10(dl) + 25.0
    
    def D_M(self, z, integrator=None):
        return self.comoving_distance(z, integrator)
    
    def D_H(self, z):
        return C / self.H(z)
    
    def D_M_over_rs(self, z, integrator=None):
        return self.D_M(z, integrator) / R_S
    
    def D_H_over_rs(self, z):
        return self.D_H(z) / R_S
    
    def I_zrec(self, integrator=None):
        if integrator is not None:
            return integrator.integral(Z_REC)
        f = lambda x: 1.0 / self.E(np.array([x]))[0]
        I, _ = integrate.quad(f, 0, Z_REC, limit=1000)
        return I
    
    def shift_R(self, integrator=None):
        return np.sqrt(self.Om) * self.I_zrec(integrator)
    
    def l_A(self, integrator=None):
        chi = self.comoving_distance(Z_REC, integrator)
        return np.pi * chi / R_S


# ============================================================
# Data
# ============================================================

HZ_DATA = np.array([
    [0.07, 69.0, 19.5], [0.09, 69.0, 12.0], [0.12, 68.6, 26.2],
    [0.17, 83.0, 8.0], [0.179, 75.0, 4.0], [0.199, 75.0, 5.0],
    [0.20, 72.9, 29.6], [0.27, 77.0, 14.0], [0.28, 88.8, 36.6],
    [0.352, 83.0, 14.0], [0.38, 81.5, 1.9], [0.4, 82.8, 1.2],
    [0.4004, 77.0, 10.2], [0.4247, 87.1, 11.2], [0.44, 82.6, 7.8],
    [0.4448, 92.8, 12.9], [0.47, 89.0, 49.6], [0.4783, 80.9, 9.0],
    [0.48, 97.0, 62.0], [0.51, 90.4, 1.9], [0.56, 93.34, 2.3],
    [0.57, 92.4, 4.5], [0.593, 104.0, 13.0], [0.6, 87.9, 6.1],
    [0.61, 97.3, 2.1], [0.68, 92.0, 8.0], [0.73, 97.3, 7.0],
    [0.78, 97.6, 4.0], [0.781, 105.0, 12.0], [0.875, 125.0, 17.0],
    [0.88, 90.0, 40.0], [0.9, 117.0, 23.0], [1.037, 154.0, 20.0],
    [1.3, 168.0, 17.0], [1.363, 160.0, 33.6], [1.43, 177.0, 18.0],
    [1.53, 140.0, 14.0], [1.75, 202.0, 40.0], [1.965, 186.5, 50.4],
])

BAO_DATA = np.array([
    [0.38,  8.46, 0.17, 24.22, 0.72, -0.505],
    [0.51, 10.75, 0.19, 22.46, 0.60, -0.447],
    [0.61, 12.32, 0.22, 21.38, 0.60, -0.442],
])


def load_pantheon(path="/tmp/pantheon_raw.txt"):
    zs, mbs, dmbs = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            try:
                z = float(parts[1])
                mb = float(parts[4])
                dmb = float(parts[5])
                if z >= 0.01:
                    zs.append(z)
                    mbs.append(mb)
                    dmbs.append(dmb)
            except (ValueError, IndexError):
                continue
    return np.array(zs), np.array(mbs), np.array(dmbs)


# ============================================================
# χ² functions (batch-able)
# ============================================================

def chi2_hz(model):
    """Vectorized H(z) χ² — no integral needed"""
    resid = model.H(HZ_DATA[:, 0]) - HZ_DATA[:, 1]
    return np.sum((resid / HZ_DATA[:, 2]) ** 2)


def chi2_sne(model, zs, mbs, dmbs, integrator):
    """SNe χ² with M marginalization — fully vectorized."""
    # Compute ALL distance moduli in ONE vectorized call
    dl = (1 + zs) * (C / model.H0) * integrator.integral(zs)  # D_L for all SNe
    mu_model = 5.0 * np.log10(dl) + 25.0
    w = 1.0 / dmbs**2
    resid = mu_model - mbs
    M_best = -np.sum(w * resid) / np.sum(w)
    chi2 = np.sum(w * (resid + M_best)**2)
    return chi2, M_best


def chi2_bao(model, integrator):
    chi2 = 0.0
    for row in BAO_DATA:
        z, dm_d, dm_e, dh_d, dh_e, corr = row
        dm_m = model.D_M_over_rs(z, integrator)
        dh_m = model.D_H_over_rs(z)
        c11 = dm_e**2; c22 = dh_e**2; c12 = corr * dm_e * dh_e
        det = c11 * c22 - c12**2
        i11 = c22 / det; i22 = c11 / det; i12 = -c12 / det
        d_dm = dm_m - dm_d; d_dh = dh_m - dh_d
        chi2 += i11*d_dm**2 + 2*i12*d_dm*d_dh + i22*d_dh**2
    return chi2


def chi2_cmb(model, integrator):
    r_m = model.shift_R(integrator)
    la_m = model.l_A(integrator)
    cov_rr = PLANCK_R_ERR**2; cov_ll = PLANCK_LA_ERR**2
    cov_rl = PLANCK_R_LA_CORR * PLANCK_R_ERR * PLANCK_LA_ERR
    det = cov_rr * cov_ll - cov_rl**2
    i_rr = cov_ll / det; i_ll = cov_rr / det; i_rl = -cov_rl / det
    d_r = r_m - PLANCK_R; d_l = la_m - PLANCK_LA
    chi2 = i_rr*d_r**2 + 2*i_rl*d_r*d_l + i_ll*d_l**2
    return chi2, r_m, la_m


def evaluate_model(model, sne_zs, sne_mbs, sne_dmbs, integrator):
    """Single model evaluation: all 4 χ² + details."""
    t0 = time.time()
    chi2_h = chi2_hz(model)
    chi2_s, M_best = chi2_sne(model, sne_zs, sne_mbs, sne_dmbs, integrator)
    chi2_b = chi2_bao(model, integrator)
    chi2_c, r_val, la_val = chi2_cmb(model, integrator)
    total = chi2_h + chi2_s + chi2_b + chi2_c
    elapsed = time.time() - t0
    return total, {
        'h': chi2_h, 's': chi2_s, 'b': chi2_b, 'c': chi2_c,
        'M': M_best, 'R': r_val, 'lA': la_val, 't': elapsed
    }


# ============================================================
# Coarse → fine grid
# ============================================================

def grid_scan(H0_range, Om_range, eps_range, zc_fixed, sne_zs, sne_mbs, sne_dmbs):
    """Coarse 3D grid, then interpolated minima."""
    from itertools import product
    H0_vals = np.linspace(*H0_range)
    Om_vals = np.linspace(*Om_range)
    eps_vals = np.linspace(*eps_range)
    
    nH = len(H0_vals)
    nO = len(Om_vals)
    nE = len(eps_vals)
    total = nH * nO * nE
    
    print(f"  Grid: {nH}×{nO}×{nE} = {total} points")
    print(f"  Building integrator (z_max=3.0, n=3000)...", end=" ", flush=True)
    
    # Build integrator for first point, reuse for all
    # But integrator depends on model params (E(z) changes)...
    # We need a budget approach: cache integrals for similar params
    # Or: since E(z) only changes slowly with Om and eps, we can 
    # approximate by building a new integrator for each param set
    
    # Actually, the integral depends on model params. We must recompute.
    # But we can vectorize: precompute E(z) grid per model and integrate once.
    
    results = []
    t_start = time.time()
    
    for i, (H0, Om, eps) in enumerate(product(H0_vals, Om_vals, eps_vals)):
        model = IDMModel(H0=H0, Om=Om, eps=eps, zc=zc_fixed)
        integrator = FastIntegrator(model)
        total_chi2, det = evaluate_model(model, sne_zs, sne_mbs, sne_dmbs, integrator)
        results.append((total_chi2, H0, Om, eps, det))
        
        if (i+1) % 200 == 0:
            elapsed = time.time() - t_start
            rate = (i+1) / elapsed
            remaining = (total - i - 1) / rate
            print(f"\r  ... {i+1}/{total} ({rate:.1f} pts/s, ~{remaining/60:.0f}m remaining)", end="", flush=True)
    
    print(f"\r  ... {total}/{total} done! ({time.time()-t_start:.0f}s total)")
    results.sort(key=lambda x: x[0])
    return results


# ============================================================
# Main
# ============================================================

def main():
    from itertools import product
    
    print("=" * 72)
    print("  IDM 聯合驗證：CMB + BAO + SNe + H(z)")
    print("  Optimized with precomputed integral interpolation")
    print("=" * 72)
    
    # Load SNe
    print("\nLoading Pantheon...", end=" ", flush=True)
    zs, mbs, dmbs = load_pantheon()
    print(f"{len(zs)} SNe (z≥0.01)")
    print(f"  z range: [{zs.min():.3f}, {zs.max():.2f}]")
    
    # ══════════════════════════════════════════════════════
    # Phase 1: ΛCDM (ε=0) — fast 2D scan
    # ══════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  Phase 1: ΛCDM 2D scan [H₀, Ω_m]")
    print("─" * 60)
    
    H0_range_lcdm = (60, 75, 16)
    Om_range_lcdm = (0.15, 0.45, 13)
    
    H0s = np.linspace(*H0_range_lcdm)
    Oms = np.linspace(*Om_range_lcdm)
    lcdm_list = []
    t0 = time.time()
    
    for H0, Om in product(H0s, Oms):
        model = IDMModel(H0=H0, Om=Om, eps=0.0, zc=0.6)
        integrator = FastIntegrator(model)
        total, det = evaluate_model(model, zs, mbs, dmbs, integrator)
        lcdm_list.append((total, H0, Om, det))
    
    lcdm_list.sort(key=lambda x: x[0])
    bl = lcdm_list[0]
    print(f"  {len(lcdm_list)} points in {time.time()-t0:.0f}s")
    print(f"  BEST: H₀={bl[1]:.1f}, Ω_m={bl[2]:.3f}, χ²={bl[0]:.1f}")
    
    # ══════════════════════════════════════════════════════
    # Phase 2: IDM coarse grid [H₀, Ω_m, ε]
    # ══════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  Phase 2: IDM coarse 3D scan [H₀, Ω_m, ε]")
    print("─" * 60)
    
    # Coarse around plausible ranges
    idm_coarse = grid_scan(
        (62, 74, 13), (0.20, 0.40, 11), (-0.20, 0.50, 15),
        0.6, zs, mbs, dmbs
    )
    bi = idm_coarse[0]
    print(f"  COARSE BEST: H₀={bi[1]:.1f}, Ω_m={bi[2]:.3f}, ε={bi[3]:.4f}, χ²={bi[0]:.1f}")
    
    # ══════════════════════════════════════════════════════
    # Phase 3: IDM fine grid around best
    # ══════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  Phase 3: IDM fine 3D scan (refinement)")
    print("─" * 60)
    
    H0_b, Om_b, eps_b = bi[1], bi[2], bi[3]
    margin = max(4.0, min(8.0, abs(eps_b)*2))
    
    idm_fine = grid_scan(
        (H0_b - 4, H0_b + 4, 11),
        (max(0.18, Om_b - 0.06), min(0.42, Om_b + 0.06), 9),
        (eps_b - 0.15, eps_b + 0.15, 13),
        0.6, zs, mbs, dmbs
    )
    
    bi2 = idm_fine[0]
    print(f"  FINE BEST: H₀={bi2[1]:.1f}, Ω_m={bi2[2]:.3f}, ε={bi2[3]:.4f}, χ²={bi2[0]:.1f}")
    
    # ══════════════════════════════════════════════════════
    # Phase 4: z_c scan
    # ══════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  Phase 4: z_c scan")
    print("─" * 60)
    
    zc_vals = np.linspace(0.2, 1.5, 14)
    zc_results = []
    for zc in zc_vals:
        model = IDMModel(H0=bi2[1], Om=bi2[2], eps=bi2[3], zc=zc)
        integrator = FastIntegrator(model)
        total, det = evaluate_model(model, zs, mbs, dmbs, integrator)
        zc_results.append((total, zc, det))
    
    zc_results.sort(key=lambda x: x[0])
    best_zc = zc_results[0][1]
    print(f"  z_c scan results:")
    print(f"    {'z_c':^6} | {'χ²_total':^10}")
    print(f"    {'─'*6}━┿━{'─'*10}")
    for chi2, zc_val, _ in zc_results:
        star = " ← best" if abs(zc_val - best_zc) < 0.01 else ""
        print(f"    {zc_val:5.2f}  │  {chi2:8.1f} {star}")
    
    # ══════════════════════════════════════════════════════
    # Phase 5: Final refinement with best z_c
    # ══════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print(f"  Phase 5: IDM refine with z_c = {best_zc:.2f}")
    print("─" * 60)
    
    final = grid_scan(
        (bi2[1] - 3, bi2[1] + 3, 9),
        (max(0.20, bi2[2] - 0.04), min(0.40, bi2[2] + 0.04), 7),
        (bi2[3] - 0.10, bi2[3] + 0.10, 11),
        best_zc, zs, mbs, dmbs
    )
    
    bf = final[0]
    
    # ══════════════════════════════════════════════════════
    # Deep final: ΛCDM with same z_c (to be fair)
    # ══════════════════════════════════════════════════════
    print("\n" + "─" * 60)
    print("  Phase 6: ΛCDM final refit")
    print("─" * 60)
    lcdm_final = []
    for H0, Om in product(np.linspace(62, 72, 9), np.linspace(0.25, 0.38, 11)):
        model = IDMModel(H0=H0, Om=Om, eps=0.0, zc=best_zc)
        integrator = FastIntegrator(model)
        total, det = evaluate_model(model, zs, mbs, dmbs, integrator)
        lcdm_final.append((total, H0, Om, det))
    
    lcdm_final.sort(key=lambda x: x[0])
    blf = lcdm_final[0]
    
    # ══════════════════════════════════════════════════════
    # FINAL REPORT
    # ══════════════════════════════════════════════════════
    
    print("\n" + "=" * 72)
    print("  FINAL RESULTS")
    print("=" * 72)
    
    # ΛCDM
    ndof_h = len(HZ_DATA) - 1
    ndof_s = len(zs) - 1
    ndof_b = len(BAO_DATA) * 2 - 1
    ndof_c = 2
    ndof_total = ndof_h + ndof_s + ndof_b + ndof_c
    
    print(f"\n  ┌{'─'*58}┐")
    print(f"  │  ΛCDM: H₀ = {blf[1]:.1f} ± —, Ω_m = {blf[2]:.3f} ± —")
    print(f"  │  χ²_H(z)  = {blf[3]['h']:8.1f}   (ndof = {ndof_h})")
    print(f"  │  χ²_SNe   = {blf[3]['s']:8.1f}   (ndof = {ndof_s})")
    print(f"  │  χ²_BAO   = {blf[3]['b']:8.1f}   (ndof = {ndof_b})")
    print(f"  │  χ²_CMB   = {blf[3]['c']:8.1f}   (ndof = {ndof_c})")
    print(f"  │  ──────────────────────────────────────────")
    print(f"  │  χ²_total = {blf[0]:8.1f}")
    print(f"  │  χ²/ndof  = {(blf[0])/(ndof_total - 2):.2f}   (2 params)")
    print(f"  │  R = {blf[3]['R']:.4f}   l_A = {blf[3]['lA']:.1f}")
    print(f"  └{'─'*58}┘")
    
    # IDM
    print(f"\n  ┌{'─'*58}┐")
    print(f"  │  IDM: H₀ = {bf[1]:.1f} ± —, Ω_m = {bf[2]:.3f} ± —")
    print(f"  │       ε = {bf[3]:.4f} ± —, z_c = {best_zc:.2f}")
    print(f"  │  χ²_H(z)  = {bf[4]['h']:8.1f}   (ndof = {ndof_h})")
    print(f"  │  χ²_SNe   = {bf[4]['s']:8.1f}   (ndof = {ndof_s})")
    print(f"  │  χ²_BAO   = {bf[4]['b']:8.1f}   (ndof = {ndof_b})")
    print(f"  │  χ²_CMB   = {bf[4]['c']:8.1f}   (ndof = {ndof_c})")
    print(f"  │  ──────────────────────────────────────────")
    print(f"  │  χ²_total = {bf[0]:8.1f}")
    print(f"  │  χ²/ndof  = {(bf[0])/(ndof_total - 4):.2f}   (4 params)")
    print(f"  │  R = {bf[4]['R']:.4f}   l_A = {bf[4]['lA']:.1f}")
    print(f"  └{'─'*58}┘")
    
    # Comparison
    dchi2 = blf[0] - bf[0]
    aic_lcdm = blf[0] + 4   # 2 params × 2
    aic_idm  = bf[0]  + 8   # 4 params × 2
    bic_lcdm = blf[0] + 2 * np.log(ndof_total)
    bic_idm  = bf[0]  + 4 * np.log(ndof_total)
    
    print(f"\n  {'─'*60}")
    print(f"  COMPARISON")
    print(f"  {'─'*60}")
    print(f"  {'Metric':<20} {'ΛCDM':>10} {'IDM':>10} {'Δ':>10}")
    print(f"  {'─'*50}")
    print(f"  {'χ² total':<20} {blf[0]:>10.1f} {bf[0]:>10.1f} {dchi2:>+10.2f}")
    print(f"  {'AIC':<20} {aic_lcdm:>10.1f} {aic_idm:>10.1f} {aic_lcdm - aic_idm:>+10.2f}")
    print(f"  {'BIC':<20} {bic_lcdm:>10.1f} {bic_idm:>10.1f} {bic_lcdm - bic_idm:>+10.2f}")
    
    verdict = ""
    if dchi2 > 2:
        # Need to check if it's worth the extra 2 params
        if aic_idm < aic_lcdm:
            verdict = "IDM 優於 ΛCDM (AIC 改善)"
        else:
            verdict = "χ² 改善但被 AIC penalty 抵消 — 不足以拒絕 ΛCDM"
    elif dchi2 > 0:
        verdict = "IDM 輕微優於 ΛCDM (Δχ² < 2)，統計上不顯著"
    else:
        verdict = "ΛCDM 仍然優於 IDM (Δχ² < 0)"
    
    print(f"\n  VERDICT: {verdict}")
    
    # Report M_best
    print(f"\n  SNe M marginalization: M_best(IDM) = {bf[4]['M']:.3f}")
    print(f"  SNe M marginalization: M_best(ΛCDM) = {blf[3]['M']:.3f}")
    
    print(f"\n  {'='*72}")
    print(f"  Run completed. Pure Python/numpy/scipy, no LLM.")
    print(f"  {'='*72}")


if __name__ == '__main__':
    main()