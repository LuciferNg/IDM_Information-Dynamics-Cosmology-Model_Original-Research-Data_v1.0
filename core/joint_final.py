"""
joint_final.py — AUDITED & CORRECTED validation
==============================================
Fixes applied:
  1. R_S_BAO = 147.09 (drag epoch)   → BAO only ✅
  2. R_S_CMB = 144.43 (recombination) → CMB l_A ✅
  3. Ω_r included in E(z) for z-integral → radiation matters at z_rec ✅
  4. All formulas verified against Planck 2018 definitions ✅
"""

import numpy as np
from scipy import integrate
from itertools import product
import time

C = 299792.458
Z_REC = 1089.92     # Planck recombination redshift
R_S_BAO = 147.09    # drag epoch sound horizon (BAO)
R_S_CMB = 144.43    # recombination sound horizon (CMB l_A)

# Planck CMB distance priors (TT+TE+EE+lowE+lensing)
PR = 1.7502; PRE = 0.0046
PLA = 301.808; PLAE = 0.090
PRLAC = 0.47

# Radiation density at z=0 (Planck T₀=2.725K)
# Ω_r = 4.15e-5 / h²  where h = H₀/100
# We compute it per model since it depends on H₀

class IDMModel:
    def __init__(self, H0=67.4, Om=0.315, eps=0.0, zc=0.7):
        self.H0 = H0
        self.Om = Om
        self.Ol = 1.0 - Om  # dark energy today
        self.eps = eps
        self.zc = zc
        h = H0 / 100.0
        self.Om_r = 4.15e-5 / h**2  # radiation density (T₀=2.725K)
    
    def f(self, z):
        """IDM dark energy evolution: f(z) = 1 - ε·(z/z_c)·exp(-z/z_c)"""
        x = np.asarray(z, dtype=float) / self.zc
        return 1.0 - self.eps * x * np.exp(-x)
    
    def E2_full(self, z):
        """E²(z) = H²(z)/H₀² including radiation"""
        z = np.asarray(z, dtype=float)
        matter = self.Om * (1+z)**3
        rad = self.Om_r * (1+z)**4
        de = self.Ol * self.f(z)
        return matter + rad + de
    
    def E(self, z):
        return np.sqrt(self.E2_full(z))
    
    def H(self, z):
        return self.H0 * self.E(z)
    
    def DC(self, z):
        """Comoving distance D_C(z) = c/H₀ × ∫₀ᶻ dz/E(z)"""
        f = lambda x: 1.0 / self.E(x)
        I, _ = integrate.quad(f, 0, z, limit=2000, epsabs=1e-12, epsrel=1e-12)
        return C / self.H0 * I
    
    def DM(self, z):
        """Transverse comoving distance (= D_C in flat space)"""
        return self.DC(z)
    
    def DH(self, z):
        """Hubble distance D_H(z) = c/H(z)"""
        return C / self.H(z)
    
    def DV(self, z):
        """Volume-averaged distance D_V = [cz D_M² / H(z)]^(1/3)"""
        DM = self.DM(z)
        H_z = self.H(z)
        return (C * z * DM**2 / H_z)**(1/3)
    
    # CMB quantities
    def I_zrec(self):
        """∫₀ᶻʳᵉᶜ dz/E(z) — full integral including radiation"""
        f = lambda x: 1.0 / self.E(x)
        I, _ = integrate.quad(f, 0, Z_REC, limit=2000, epsabs=1e-12, epsrel=1e-12)
        return I
    
    def shift_R(self):
        """R = √Ω_m × ∫₀ᶻʳᵉᶜ dz/E(z)"""
        return np.sqrt(self.Om) * self.I_zrec()
    
    def l_A(self):
        """l_A = π × D_M(z_rec) / r*  [using r* = recombination sound horizon]"""
        return np.pi * self.DC(Z_REC) / R_S_CMB

# ================================================================
# Data
# ================================================================

# H(z) Cosmic Chronometers (39 pts)
HZ = np.array([
    [0.07,69,19.5],[0.09,69,12],[0.12,68.6,26.2],[0.17,83,8],[0.179,75,4],
    [0.199,75,5],[0.20,72.9,29.6],[0.27,77,14],[0.28,88.8,36.6],[0.352,83,14],
    [0.38,81.5,1.9],[0.4,82.8,1.2],[0.4004,77,10.2],[0.4247,87.1,11.2],
    [0.44,82.6,7.8],[0.4448,92.8,12.9],[0.47,89,49.6],[0.4783,80.9,9],
    [0.48,97,62],[0.51,90.4,1.9],[0.56,93.34,2.3],[0.57,92.4,4.5],
    [0.593,104,13],[0.6,87.9,6.1],[0.61,97.3,2.1],[0.68,92,8],[0.73,97.3,7],
    [0.78,97.6,4],[0.781,105,12],[0.875,125,17],[0.88,90,40],[0.9,117,23],
    [1.037,154,20],[1.3,168,17],[1.363,160,33.6],[1.43,177,18],[1.53,140,14],
    [1.75,202,40],[1.965,186.5,50.4],
])

# DESI DR2 BAO — from Cobaya repo (DESI 2024 DR2)
# 12-vector: [DV@0.295, DM@0.510, DH@0.510, DM@0.706, DH@0.706,
#             DM@0.930, DH@0.930, DM@1.317, DH@1.317, DV@1.491,
#             DM@2.330, DH@2.330]
BAO_MEAN = np.array([
    7.92512927, 13.62003080, 20.98334647,
    16.84645313, 20.07872919,
    21.70841761, 17.87612922,
    27.78720817, 13.82372285,
    26.07217182,
    39.70838281, 8.52256583,
])
BAO_COV = np.array([
    [2.2723e-02,0,0,0,0,0,0,0,0,0,0,0],
    [0,6.3466e-02,-6.8534e-02,0,0,0,0,0,0,0,0,0],
    [0,-6.8534e-02,3.7297e-01,0,0,0,0,0,0,0,0,0],
    [0,0,0,1.0198e-01,-7.9940e-02,0,0,0,0,0,0,0],
    [0,0,0,-7.9940e-02,3.5445e-01,0,0,0,0,0,0,0],
    [0,0,0,0,0,7.9568e-02,-3.8011e-02,0,0,0,0,0],
    [0,0,0,0,0,-3.8011e-02,1.1994e-01,0,0,0,0,0],
    [0,0,0,0,0,0,0,4.7657e-01,-1.2941e-01,0,0,0],
    [0,0,0,0,0,0,0,-1.2941e-01,1.7827e-01,0,0,0],
    [0,0,0,0,0,0,0,0,0,4.4713e-01,0,0],
    [0,0,0,0,0,0,0,0,0,0,8.8975e-01,-7.6948e-02],
    [0,0,0,0,0,0,0,0,0,0,-7.6948e-02,2.9186e-02],
])
BAO_Z = [0.295,0.510,0.510,0.706,0.706,0.930,0.930,1.317,1.317,1.491,2.330,2.330]
BAO_T = ['DV','DM','DH','DM','DH','DM','DH','DM','DH','DV','DM','DH']
BAO_ICOV = np.linalg.inv(BAO_COV)

def load_sne(path="/tmp/pantheon_raw.txt"):
    zs, mbs, dmbs = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p = line.strip().split()
            if len(p) < 6: continue
            try:
                z = float(p[1])
                if z >= 0.01:
                    zs.append(z); mbs.append(float(p[4])); dmbs.append(float(p[5]))
            except: continue
    return np.array(zs), np.array(mbs), np.array(dmbs)


# ================================================================
# χ² functions
# ================================================================

def chi2_hz(m):
    r = m.H(HZ[:,0]) - HZ[:,1]
    return np.sum((r/HZ[:,2])**2)

def chi2_sne(m, zs, mbs, dmbs):
    mu = np.array([5*np.log10((1+z)*m.DC(z))+25 for z in zs])
    w = 1/dmbs**2
    resid = mu - mbs
    Mb = -np.sum(w*resid)/np.sum(w)
    return np.sum(w*(resid+Mb)**2), Mb

def chi2_bao(m):
    pred = np.zeros(12)
    for i in range(12):
        z = BAO_Z[i]
        if BAO_T[i] == 'DM': pred[i] = m.DM(z)/R_S_BAO
        elif BAO_T[i] == 'DH': pred[i] = m.DH(z)/R_S_BAO
        elif BAO_T[i] == 'DV': pred[i] = m.DV(z)/R_S_BAO
    d = pred - BAO_MEAN
    return d @ BAO_ICOV @ d

def chi2_cmb(m):
    r = m.shift_R()
    la = m.l_A()
    # 2D Gaussian prior with covariance
    var_r = PRE**2; var_la = PLAE**2; cov_rl = PRLAC*PRE*PLAE
    det = var_r*var_la - cov_rl**2
    i_r = var_la/det; i_la = var_r/det; i_rl = -cov_rl/det
    dr = r - PR; dla = la - PLA
    chi2 = i_r*dr**2 + 2*i_rl*dr*dla + i_la*dla**2
    return chi2, r, la

def eval_m(m, zs, mbs, dmbs):
    h = chi2_hz(m)
    s, Mb = chi2_sne(m, zs, mbs, dmbs)
    b = chi2_bao(m)
    c, r, la = chi2_cmb(m)
    return h+s+b+c, {'h':h,'s':s,'b':b,'c':c,'M':Mb,'R':r,'lA':la}


# ================================================================
# Main
# ================================================================

print("="*72)
print("  IDM JOINT VALIDATION [AUDITED]")
print("  DESI DR2 BAO + Planck CMB + Pantheon SNe + H(z) CC")
print("="*72)

zs, mbs, dmbs = load_sne()
print(f"  SNe: {len(zs)} points  |  H(z): {len(HZ)} pts  |  BAO: 12-vector (DESI DR2)")

t0 = time.time()

# === Phase 1: Planck ΛCDM baseline (verification) ===
print("\n─── Planck 2018 ΛCDM baseline ───")
m = IDMModel(H0=67.36, Om=0.3153, eps=0.0, zc=0.7)
t, d = eval_m(m, zs, mbs, dmbs)
print(f"  H₀={m.H0}, Ω_m={m.Om}, Ω_r={m.Om_r:.3e}")
print(f"  D_C(z_rec)={m.DC(Z_REC):.0f} Mpc")
print(f"  R={d['R']:.4f}  l_A={d['lA']:.1f}")
print(f"  χ²_total={t:.1f} (HZ={d['h']:.1f}, SNe={d['s']:.1f}, BAO={d['b']:.1f}, CMB={d['c']:.1f})")

# === Phase 2: ΛCDM grid ===
print("\n─── ΛCDM grid scan ───")
lcdm = []
pts = 0
for H0 in np.linspace(60, 75, 16):
    for Om in np.linspace(0.18, 0.42, 13):
        m = IDMModel(H0=H0, Om=Om, eps=0.0, zc=0.7)
        t, d = eval_m(m, zs, mbs, dmbs)
        lcdm.append((t, H0, Om, d))
        pts += 1
lcdm.sort(key=lambda x: x[0])
bl = lcdm[0]
print(f"  {pts} points | BEST: H₀={bl[1]:.1f}, Ω_m={bl[2]:.3f}, χ²={bl[0]:.1f}")

# === Phase 3: IDM coarse 3D ===
print("\n─── IDM coarse 3D ───")
idm_list = []
for H0 in np.linspace(60, 75, 11):
    for Om in np.linspace(0.18, 0.42, 9):
        for eps in np.linspace(-0.30, 0.30, 13):
            m = IDMModel(H0=H0, Om=Om, eps=eps, zc=0.7)
            t, d = eval_m(m, zs, mbs, dmbs)
            idm_list.append((t, H0, Om, eps, d))
idm_list.sort(key=lambda x: x[0])
bi = idm_list[0]
print(f"  BEST: H₀={bi[1]:.1f}, Ω_m={bi[2]:.3f}, ε={bi[3]:+.4f}, χ²={bi[0]:.1f}")

# === Phase 4: z_c scan ===
print("\n─── z_c scan ───")
zc_list = []
for zc in np.linspace(0.2, 1.5, 14):
    m = IDMModel(H0=bi[1], Om=bi[2], eps=bi[3], zc=zc)
    t, d = eval_m(m, zs, mbs, dmbs)
    zc_list.append((t, zc, d))
zc_list.sort(key=lambda x: x[0])
bz = zc_list[0]
print(f"  BEST z_c={bz[1]:.2f}, χ²={bz[0]:.1f}")

# === Phase 5: IDM refined ===
print(f"\n─── IDM refined (z_c={bz[1]:.2f}) ───")
idm_f = []
for H0 in np.linspace(bi[1]-4, bi[1]+4, 9):
    for Om in np.linspace(max(0.18,bi[2]-0.06), min(0.42,bi[2]+0.06), 7):
        for eps in np.linspace(bi[3]-0.12, bi[3]+0.12, 11):
            m = IDMModel(H0=H0, Om=Om, eps=eps, zc=bz[1])
            t, d = eval_m(m, zs, mbs, dmbs)
            idm_f.append((t, H0, Om, eps, d))
idm_f.sort(key=lambda x: x[0])
bf = idm_f[0]

# ΛCDM final with same z_c
lcdm_f = []
for H0 in np.linspace(62, 72, 9):
    for Om in np.linspace(0.24, 0.38, 10):
        m = IDMModel(H0=H0, Om=Om, eps=0.0, zc=bz[1])
        t, d = eval_m(m, zs, mbs, dmbs)
        lcdm_f.append((t, H0, Om, d))
lcdm_f.sort(key=lambda x: x[0])
blf = lcdm_f[0]

elapsed = time.time() - t0

# ================================================================
# REPORT
# ================================================================
print(f"\n{'='*72}")
print(f"  FINAL RESULTS (audited code, elapsed={elapsed:.0f}s)")
print(f"{'='*72}")

ndof_h = len(HZ)-1; ndof_s = len(zs)-1; ndof_b = 12-4; ndof_c = 2

def show(label, t, d, np_):
    ndof = ndof_h+ndof_s+ndof_b+ndof_c - np_
    print(f"\n  {label}")
    print(f"  χ²_H(z)={d['h']:7.1f}  χ²_SNe={d['s']:7.1f}")
    print(f"  χ²_BAO={d['b']:7.1f}  χ²_CMB={d['c']:7.1f}")
    print(f"  χ²_total = {t:6.1f}  (ndof={ndof}, χ²/dof={t/ndof:.2f})")
    print(f"  R={d['R']:.4f}, l_A={d['lA']:.1f}, M={d['M']:.3f}")

show(f"ΛCDM  H₀={blf[1]:.1f}  Ω_m={blf[2]:.3f}", blf[0], blf[3], 2)
show(f"IDM   H₀={bf[1]:.1f}  Ω_m={bf[2]:.3f}  ε={bf[3]:+.4f}  z_c={bz[1]:.2f}",
     bf[0], bf[4], 4)

dchi = blf[0] - bf[0]
aic_l = blf[0] + 4; aic_i = bf[0] + 8
bic_l = blf[0] + 2*np.log(ndof_h+ndof_s+ndof_b+ndof_c)
bic_i = bf[0] + 4*np.log(ndof_h+ndof_s+ndof_b+ndof_c)

print(f"\n  {'':>14} {'ΛCDM':>10} {'IDM':>12} {'Δ':>10}")
print(f"  {'─'*46}")
print(f"  {'χ² total':>12} {blf[0]:>10.1f} {bf[0]:>12.1f} {dchi:>+10.2f}")
print(f"  {'AIC':>12} {aic_l:>10.1f} {aic_i:>12.1f} {aic_l-aic_i:>+10.2f}")
print(f"  {'BIC':>12} {bic_l:>10.1f} {bic_i:>12.1f} {bic_l-bic_i:>+10.2f}")

print(f"\n  Per-dataset Δχ²:")
for lb, k in [('H(z)','h'),('SNe','s'),('BAO','b'),('CMB','c')]:
    dv = bf[4][k] - blf[3][k]
    better = "IDM" if dv<0 else "ΛCDM"
    print(f"    {lb:5s}  Δχ²={dv:+7.2f}  →  {better}")

# Check: Planck params consistency
print(f"\n  Consistency check (Planck params, Ω_r included):")
for H0, Om, label in [(67.36,0.3153,'Planck')]:
    m = IDMModel(H0=H0, Om=Om, eps=0.0, zc=0.7)
    _, d = eval_m(m, zs, mbs, dmbs)
    print(f"    {label}: R={d['R']:.4f} (Planck {PR}±{PRE}), l_A={d['lA']:.1f} (Planck {PLA}±{PLAE})")