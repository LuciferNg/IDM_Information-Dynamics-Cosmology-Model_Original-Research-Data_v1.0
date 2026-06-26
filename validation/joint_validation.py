"""
joint_validation.py — CMB + BAO + SNe + H(z) 聯合驗證
=====================================================
IDM model: H²(z) = H₀²[Ω_m(1+z)³ + (1-Ω_m)·(1 - ε·z/z_c·exp(-z/z_c))]

Datasets:
  1. H(z) Cosmic Chronometers — Moresco+2016, Stern+2010, Farooq+2017 (40 pts)
  2. Pantheon SNe — Scolnic+2018 (1048 SNe, diagonal χ², M marginalized)
  3. SDSS DR12 BAO — Alam+2017 (3 bins: D_M/r_d + D_H/r_d)
  4. CMB distance priors — Planck 2018 (shift R + acoustic scale l_A)

Pure Python/numpy/scipy, 無 LLM 參與計算。
"""

import numpy as np
from scipy import integrate
import time
from itertools import product

# ============================================================
# 物理常數
# ============================================================
C = 299792.458  # km/s
Z_REC = 1089.9  # recombination redshift
R_S = 147.09    # comoving sound horizon at z_rec (Mpc, Planck 2018)

# ============================================================
# IDM 模型核心
# ============================================================

class IDMModel:
    """IDM cosmology: H²(z) = H₀²[Ω_m(1+z)³ + (1-Ω_m)·f(z)]
    where f(z) = 1 - ε·(z/z_c)·exp(-z/z_c)
    """
    def __init__(self, H0=67.4, Om=0.315, eps=0.0, zc=0.6):
        self.H0 = H0
        self.Om = Om
        self.Ol = 1.0 - Om
        self.eps = eps
        self.zc = zc
    
    def f_idm(self, z):
        """Dark energy evolution function"""
        x = z / self.zc
        return 1.0 - self.eps * x * np.exp(-x)
    
    def E2(self, z):
        """E²(z) = H²(z) / H₀²"""
        return self.Om * (1+z)**3 + self.Ol * self.f_idm(z)
    
    def E(self, z):
        return np.sqrt(self.E2(z))
    
    def H(self, z):
        return self.H0 * self.E(z)
    
    def comoving_distance(self, z):
        """χ(z) = c/H₀ × ∫₀ᶻ dz'/E(z') — Mpc"""
        f = lambda x: 1.0 / self.E(x)
        integral, _ = integrate.quad(f, 0, z, limit=1000, epsabs=1e-12, epsrel=1e-12)
        return C / self.H0 * integral
    
    def lum_distance(self, z):
        """D_L(z) = (1+z) × χ(z) — Mpc"""
        return (1+z) * self.comoving_distance(z)
    
    def distance_modulus(self, z):
        """μ(z) = 5 log₁₀(D_L/Mpc) + 25"""
        dl = self.lum_distance(z)
        return 5.0 * np.log10(dl) + 25.0
    
    # BAO quantities
    def D_M(self, z):
        """Transverse comoving distance (Mpc) = χ(z) in flat universe"""
        return self.comoving_distance(z)
    
    def D_H(self, z):
        """Hubble distance (Mpc) = c/H(z)"""
        return C / self.H(z)
    
    def D_V(self, z):
        """Volume-averaged distance: [cz D_M² / D_H × (1+z)²? No: standard D_V = [cz(1+z)² D_A² / H(z)]^(1/3)
        = [cz D_M(z)² / (D_H(z))]^(1/3) because D_A = D_M/(1+z)
        Actually D_V = [cz(1+z)² D_A² / H(z)]^(1/3) = [cz D_M(z)² / H(z) × c⁻¹? No...
        
        Standard definition from Eisenstein+2005:
        D_V = [cz D_M(z)² / H(z)]^(1/3)
        where H(z) is in km/s/Mpc, c cancels:
        D_V = [cz × D_M² / H(z)]^(1/3)
        """
        return (C * z * self.D_M(z)**2 / self.H(z))**(1/3)
    
    def D_M_over_rs(self, z):
        return self.D_M(z) / R_S
    
    def D_H_over_rs(self, z):
        return self.D_H(z) / R_S
    
    def D_V_over_rs(self, z):
        return self.D_V(z) / R_S
    
    # CMB quantities
    def I_zrec(self):
        """∫₀ᶻʳᵉᶜ dz / E(z)"""
        f = lambda x: 1.0 / self.E(x)
        integral, _ = integrate.quad(f, 0, Z_REC, limit=2000, epsabs=1e-12, epsrel=1e-12)
        return integral
    
    def shift_R(self):
        """R = √Ω_m × ∫₀ᶻʳᵉᶜ dz/E(z)"""
        return np.sqrt(self.Om) * self.I_zrec()
    
    def l_A(self):
        """l_A = π × χ(z_rec) / r_s"""
        chi = self.comoving_distance(Z_REC)
        return np.pi * chi / R_S


# ============================================================
# H(z) Cosmic Chronometers Data
# ============================================================
HZ_DATA = np.array([
    [0.07, 69.0, 19.5],
    [0.09, 69.0, 12.0],
    [0.12, 68.6, 26.2],
    [0.17, 83.0, 8.0],
    [0.179, 75.0, 4.0],
    [0.199, 75.0, 5.0],
    [0.20, 72.9, 29.6],
    [0.27, 77.0, 14.0],
    [0.28, 88.8, 36.6],
    [0.352, 83.0, 14.0],
    [0.38, 81.5, 1.9],
    [0.4, 82.8, 1.2],
    [0.4004, 77.0, 10.2],
    [0.4247, 87.1, 11.2],
    [0.44, 82.6, 7.8],
    [0.4448, 92.8, 12.9],
    [0.47, 89.0, 49.6],
    [0.4783, 80.9, 9.0],
    [0.48, 97.0, 62.0],
    [0.51, 90.4, 1.9],
    [0.56, 93.34, 2.3],
    [0.57, 92.4, 4.5],
    [0.593, 104.0, 13.0],
    [0.6, 87.9, 6.1],
    [0.61, 97.3, 2.1],
    [0.68, 92.0, 8.0],
    [0.73, 97.3, 7.0],
    [0.78, 97.6, 4.0],
    [0.781, 105.0, 12.0],
    [0.875, 125.0, 17.0],
    [0.88, 90.0, 40.0],
    [0.9, 117.0, 23.0],
    [1.037, 154.0, 20.0],
    [1.3, 168.0, 17.0],
    [1.363, 160.0, 33.6],
    [1.43, 177.0, 18.0],
    [1.53, 140.0, 14.0],
    [1.75, 202.0, 40.0],
    [1.965, 186.5, 50.4],
])


def chi2_hz(model):
    """χ² for H(z) Cosmic Chronometers"""
    residuals = model.H(HZ_DATA[:, 0]) - HZ_DATA[:, 1]
    return np.sum((residuals / HZ_DATA[:, 2]) ** 2)


# ============================================================
# Pantheon SNe Data
# ============================================================

def load_pantheon_simple(path="/tmp/pantheon_raw.txt"):
    """Load Pantheon SN data: zcmb, mb, dmb. z_min=0.01 to avoid peculiar velocities."""
    zs, mbs, dmbs = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 6:
                continue
            try:
                z = float(parts[1])  # zcmb
                mb = float(parts[4])  # mb (apparent magnitude)
                dmb = float(parts[5]) # error
                if z >= 0.01:  # exclude z<0.01 (peculiar velocity dominated)
                    zs.append(z)
                    mbs.append(mb)
                    dmbs.append(dmb)
            except (ValueError, IndexError):
                continue
    return np.array(zs), np.array(mbs), np.array(dmbs)


def chi2_sne(model, zs, mbs, dmbs):
    """
    χ² for SNe with analytical M marginalization.
    
    μ_model(z) = 5 log₁₀(D_L(z)) + 25
    m_obs = μ_model(z) + M + noise
    
    Marginalize over M:
    M_best = -Σ w_i (μ_i - m_i) / Σ w_i   where w_i = 1/σ_i²
    χ²_min = Σ w_i (μ_i - m_i + M_best)²
    """
    mu_model = np.array([model.distance_modulus(z) for z in zs])
    w = 1.0 / dmbs**2  # weights
    resid = mu_model - mbs
    M_best = -np.sum(w * resid) / np.sum(w)
    chi2 = np.sum(w * (resid + M_best)**2)
    return chi2, M_best


# ============================================================
# BAO Data — SDSS DR12 Alam+2017
# ============================================================
# Each row: [z, D_M/r_d, err_DM, D_H/r_d, err_DH, corr]
BAO_DATA = np.array([
    [0.38,  8.46, 0.17, 24.22, 0.72, -0.505],
    [0.51, 10.75, 0.19, 22.46, 0.60, -0.447],
    [0.61, 12.32, 0.22, 21.38, 0.60, -0.442],
])


def chi2_bao(model):
    """
    χ² for SDSS DR12 BAO (D_M/r_d + D_H/r_d, with covariance).
    Uses full 2x2 covariance per bin.
    """
    chi2 = 0.0
    for row in BAO_DATA:
        z, dm_data, dm_err, dh_data, dh_err, corr = row
        
        dm_model = model.D_M_over_rs(z)
        dh_model = model.D_H_over_rs(z)
        
        # Covariance matrix
        cov_11 = dm_err**2
        cov_22 = dh_err**2
        cov_12 = corr * dm_err * dh_err
        det = cov_11 * cov_22 - cov_12**2
        
        # Inverse covariance (2x2)
        inv_11 = cov_22 / det
        inv_22 = cov_11 / det
        inv_12 = -cov_12 / det
        
        d_dm = dm_model - dm_data
        d_dh = dh_model - dh_data
        
        chi2 += inv_11 * d_dm**2 + 2*inv_12 * d_dm * d_dh + inv_22 * d_dh**2
    
    return chi2


# ============================================================
# CMB Distance Priors — Planck 2018
# ============================================================
PLANCK_R = 1.7502
PLANCK_R_ERR = 0.0046
PLANCK_LA = 301.808
PLANCK_LA_ERR = 0.090
PLANCK_R_LA_CORR = 0.47  # approximate correlation


def chi2_cmb(model):
    """χ² for Planck 2018 distance priors (R + l_A, with covariance)."""
    r_model = model.shift_R()
    la_model = model.l_A()
    
    cov_rr = PLANCK_R_ERR**2
    cov_ll = PLANCK_LA_ERR**2
    cov_rl = PLANCK_R_LA_CORR * PLANCK_R_ERR * PLANCK_LA_ERR
    det = cov_rr * cov_ll - cov_rl**2
    
    inv_rr = cov_ll / det
    inv_ll = cov_rr / det
    inv_rl = -cov_rl / det
    
    d_r = r_model - PLANCK_R
    d_l = la_model - PLANCK_LA
    
    chi2 = inv_rr * d_r**2 + 2*inv_rl * d_r * d_l + inv_ll * d_l**2
    return chi2, r_model, la_model


# ============================================================
# 聯合 χ²
# ============================================================

def joint_chi2(model, sne_zs, sne_mbs, sne_dmbs):
    """All 4 probes combined."""
    chi2_h = chi2_hz(model)
    chi2_s, M_best = chi2_sne(model, sne_zs, sne_mbs, sne_dmbs)
    chi2_b = chi2_bao(model)
    chi2_c, r_val, la_val = chi2_cmb(model)
    
    total = chi2_h + chi2_s + chi2_b + chi2_c
    return total, {
        'χ²_H(z)': chi2_h,
        'χ²_SNe': chi2_s,
        'χ²_BAO': chi2_b,
        'χ²_CMB': chi2_c,
        'M_best': M_best,
        'R': r_val,
        'l_A': la_val,
        'ndof_HZ': len(HZ_DATA) - 1,
        'ndof_SNe': len(sne_zs) - 1,
        'ndof_BAO': len(BAO_DATA) * 2 - 1,
        'ndof_CMB': 2,
    }


# ============================================================
# Grid Search
# ============================================================

def grid_search(sne_zs, sne_mbs, sne_dmbs,
                H0_range=(60, 75, 31),
                Om_range=(0.20, 0.40, 21),
                eps_range=(-0.30, 0.50, 41),
                zc_fixed=0.6):
    """3D grid scan: [H₀, Ω_m, ε] with fixed z_c."""
    H0_vals = np.linspace(*H0_range)
    Om_vals = np.linspace(*Om_range)
    eps_vals = np.linspace(*eps_range)
    
    # For ΛCDM (ε=0), 2D scan
    print("Grid scanning ΛCDM (ε=0)... [H₀, Ω_m]")
    lcdm_results = []
    for H0, Om in product(H0_vals, Om_vals):
        model = IDMModel(H0=H0, Om=Om, eps=0.0, zc=zc_fixed)
        total, details = joint_chi2(model, sne_zs, sne_mbs, sne_dmbs)
        lcdm_results.append((total, H0, Om, details))
    
    lcdm_results.sort(key=lambda x: x[0])
    best_lcdm = lcdm_results[0]
    
    print(f"  ΛCDM best: H₀={best_lcdm[1]:.1f}, Ω_m={best_lcdm[2]:.3f}, χ²={best_lcdm[0]:.1f}")
    
    # For IDM, 3D scan [H₀, Ω_m, ε]
    print("Grid scanning IDM (ε free)... [H₀, Ω_m, ε]")
    idm_results = []
    t0 = time.time()
    for i, (H0, Om, eps) in enumerate(product(H0_vals, Om_vals, eps_vals)):
        model = IDMModel(H0=H0, Om=Om, eps=eps, zc=zc_fixed)
        total, details = joint_chi2(model, sne_zs, sne_mbs, sne_dmbs)
        idm_results.append((total, H0, Om, eps, details))
        if (i+1) % 500 == 0:
            elapsed = time.time() - t0
            print(f"  ... {i+1}/{len(H0_vals)*len(Om_vals)*len(eps_vals)} done ({elapsed:.0f}s)")
    
    idm_results.sort(key=lambda x: x[0])
    best_idm = idm_results[0]
    
    print(f"  IDM best: H₀={best_idm[1]:.1f}, Ω_m={best_idm[2]:.3f}, ε={best_idm[3]:.4f}, χ²={best_idm[0]:.1f}")
    
    return {
        'lcdm': lcdm_results,
        'idm': idm_results,
        'best_lcdm': best_lcdm,
        'best_idm': best_idm,
        'H0_vals': H0_vals,
        'Om_vals': Om_vals,
        'eps_vals': eps_vals,
    }


# ============================================================
# 掃描 z_c
# ============================================================

def scan_zc(sne_zs, sne_mbs, sne_dmbs, best_H0, best_Om, best_eps,
            zc_range=(0.2, 1.5, 27)):
    """Vary z_c around the best-fit to find its optimum."""
    zc_vals = np.linspace(*zc_range)
    results = []
    for zc in zc_vals:
        model = IDMModel(H0=best_H0, Om=best_Om, eps=best_eps, zc=zc)
        total, details = joint_chi2(model, sne_zs, sne_mbs, sne_dmbs)
        results.append((total, zc, details))
    
    results.sort(key=lambda x: x[0])
    return results


# ============================================================
# 報告
# ============================================================

def print_header(title):
    print(f"\n{'='*72}")
    print(f"  {title}")
    print(f"{'='*72}")

def print_result(label, total, details, nparams):
    ndof_total = details['ndof_HZ'] + details['ndof_SNe'] + details['ndof_BAO'] + details['ndof_CMB'] - nparams
    print(f"\n  {'─'*60}")
    print(f"  {label}")
    print(f"  {'─'*60}")
    print(f"  χ²_H(z)  = {details['χ²_H(z)']:8.1f}  (ndof={details['ndof_HZ']})")
    print(f"  χ²_SNe   = {details['χ²_SNe']:8.1f}  (ndof={details['ndof_SNe']})")
    print(f"  χ²_BAO   = {details['χ²_BAO']:8.1f}  (ndof={details['ndof_BAO']})")
    print(f"  χ²_CMB   = {details['χ²_CMB']:8.1f}  (ndof={details['ndof_CMB']})")
    print(f"  {'─'*60}")
    print(f"  χ²_total = {total:8.1f}")
    print(f"  ndof     = {ndof_total}")
    print(f"  χ²/ndof  = {total/ndof_total:.2f}")
    print(f"  R        = {details['R']:.4f}  (Planck: {PLANCK_R} ± {PLANCK_R_ERR})")
    print(f"  l_A      = {details['l_A']:.2f}  (Planck: {PLANCK_LA} ± {PLANCK_LA_ERR})")


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print_header("IDM 聯合驗證：CMB + BAO + SNe + H(z)")
    print("  " + "─"*56)
    print("  IDM model: H²(z)=H₀²[Ω_m(1+z)³+(1-Ω_m)(1-ε·z/z_c·exp(-z/z_c))]")
    print("  z_rec = 1089.9, r_s = 147.09 Mpc")
    print()
    
    # Load SNe
    print("Loading Pantheon SNe data...")
    zs, mbs, dmbs = load_pantheon_simple()
    print(f"  {len(zs)} SNe after z≥0.01 cut")
    print(f"  z range: [{zs.min():.4f}, {zs.max():.2f}]")
    
    # Phase 1: Grid scan (z_c=0.6 fixed)
    print_header("Phase 1: 3D Grid Scan · z_c = 0.6")
    results = grid_search(zs, mbs, dmbs)
    
    best_lcdm = results['best_lcdm']
    best_idm = results['best_idm']
    
    print_header("Phase 2: Results Summary")
    
    # ΛCDM
    print_result(
        "ΛCDM (ε=0): H₀={:.1f}, Ω_m={:.3f}".format(best_lcdm[1], best_lcdm[2]),
        best_lcdm[0], best_lcdm[3], nparams=2
    )
    
    # IDM
    print_result(
        "IDM: H₀={:.1f}, Ω_m={:.3f}, ε={:.4f}".format(best_idm[1], best_idm[2], best_idm[3]),
        best_idm[0], best_idm[4], nparams=4
    )
    
    # Δ
    delta_chi2 = best_lcdm[0] - best_idm[0]
    print(f"\n  Δχ² (ΛCDM − IDM) = {delta_chi2:.2f}")
    if delta_chi2 > 0:
        print(f"  → IDM 比 ΛCDM 好 (Δχ² > 0)")
    else:
        print(f"  → ΛCDM 比 IDM 好 (Δχ² < 0)")
    
    # AIC comparison (IDM has 2 extra params: ε, z_c)
    aic_lcdm = best_lcdm[0] + 2 * 2  # H₀, Ω_m
    aic_idm = best_idm[0] + 2 * 4    # H₀, Ω_m, ε, z_c
    print(f"  AIC_ΛCDM = {aic_lcdm:.1f}")
    print(f"  AIC_IDM  = {aic_idm:.1f}")
    print(f"  ΔAIC (ΛCDM − IDM) = {aic_lcdm - aic_idm:.2f}")
    
    # Phase 3: z_c scan
    print_header("Phase 3: z_c Scan around IDM best-fit")
    zc_results = scan_zc(
        zs, mbs, dmbs,
        best_idm[1], best_idm[2], best_idm[3]
    )
    
    best_zc = zc_results[0][1]
    best_zc_chi2 = zc_results[0][0]
    print(f"  z_c best = {best_zc:.2f}  (χ² = {best_zc_chi2:.1f})")
    print(f"  z_c scan:")
    print(f"  {'z_c':^8} {'χ²_total':^10}")
    print(f"  {'─'*8} {'─'*10}")
    for chi2, zc_val, det in zc_results:
        print(f"  {zc_val:6.2f}  {chi2:8.1f}")
    
    # Phase 4: Re-check with best z_c
    print_header("Phase 4: IDM with z_c = {:.2f}".format(best_zc))
    # Recompute best-fit at the optimal z_c
    H0_range = (best_idm[1]-5, best_idm[1]+5, 21)
    Om_range = (best_idm[2]-0.05, best_idm[2]+0.05, 11)
    eps_range = (best_idm[3]-0.10, best_idm[3]+0.10, 21)
    
    refined = grid_search(zs, mbs, dmbs,
                          H0_range=H0_range,
                          Om_range=Om_range,
                          eps_range=eps_range,
                          zc_fixed=best_zc)
    
    br = refined['best_idm']
    print_result(
        "IDM (refined): H₀={:.1f}, Ω_m={:.3f}, ε={:.4f}".format(br[1], br[2], br[3]),
        br[0], br[4], nparams=4
    )
    
    delta_chi2_refined = best_lcdm[0] - br[0]
    print(f"\n  Δχ² (ΛCDM − IDM) = {delta_chi2_refined:.2f}")
    
    aic_lcdm = best_lcdm[0] + 4
    aic_idm_refined = br[0] + 8
    print(f"  AIC_ΛCDM = {aic_lcdm:.1f}")
    print(f"  AIC_IDM  = {aic_idm_refined:.1f}")
    print(f"  ΔAIC (ΛCDM − IDM) = {aic_lcdm - aic_idm_refined:.2f}")
    
    # Summary
    print_header("Final Verdict")
    print(f"""
  ┌─────────────────────────────────────────────────────────┐
  │  Dataset             │  IDM │  ΛCDM │  Δ (I−L)          │
  ├──────────────────────┼──────┼───────┼───────────────────┤
  │  H(z) CC (40 pts)    │ {br[4]['χ²_H(z)']:4.1f} │ {best_lcdm[3]['χ²_H(z)']:4.1f}  │ {br[4]['χ²_H(z)']-best_lcdm[3]['χ²_H(z)']:+6.2f}   │
  │  Pantheon SNe ({len(zs)} pts) │ {br[4]['χ²_SNe']:4.1f} │ {best_lcdm[3]['χ²_SNe']:4.1f}  │ {br[4]['χ²_SNe']-best_lcdm[3]['χ²_SNe']:+6.2f}   │
  │  SDSS BAO (3 bins)   │ {br[4]['χ²_BAO']:4.1f} │ {best_lcdm[3]['χ²_BAO']:4.1f}  │ {br[4]['χ²_BAO']-best_lcdm[3]['χ²_BAO']:+6.2f}   │
  │  CMB prior (2 param) │ {br[4]['χ²_CMB']:4.1f} │ {best_lcdm[3]['χ²_CMB']:4.1f}  │ {br[4]['χ²_CMB']-best_lcdm[3]['χ²_CMB']:+6.2f}   │
  ├──────────────────────┼──────┼───────┼───────────────────┤
  │  Total χ²            │ {br[0]:6.1f} │ {best_lcdm[0]:6.1f}  │ {br[0]-best_lcdm[0]:+6.2f}   │
  │  Parameters          │  4   │  2    │  +2                │
  │  AIC                 │ {aic_idm_refined:6.1f} │ {aic_lcdm:6.1f}  │ {aic_idm_refined-aic_lcdm:+6.2f}  │
  └─────────────────────────────────────────────────────────┘
  
  Best-fit parameters:
    ΛCDM: H₀ = {best_lcdm[1]:.1f}, Ω_m = {best_lcdm[2]:.3f}
    IDM:  H₀ = {br[1]:.1f}, Ω_m = {br[2]:.3f}, ε = {br[3]:.4f}, z_c = {best_zc:.2f}
""")