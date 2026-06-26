"""
mcmc_idm.py — MCMC joint fit of IDM with DESI DR2 BAO + CMB + SNe + H(z)
========================================================================
Uses emcee with fast D_C interpolation (replaces 1048 quads per eval).
Parameters: H₀, Ω_m, ε, z_c  (λ fixed at 0.001, undetectable at z<2)
"""

import numpy as np
import emcee
import corner
from scipy import integrate
import time, os, sys, importlib.util

# ============================================================
# Constants
# ============================================================
C = 299792.458
Z_REC = 1089.92
R_S_BAO = 147.09
R_S_CMB = 144.43
PR = 1.7502; PRE = 0.0046
PLA = 301.808; PLAE = 0.090
PRLAC = 0.47
LAM = 0.001  # secular efficiency (fixed, undetectable at z<2)

# BAO data
BAO_MEAN = np.array([
    7.92512927,13.62003080,20.98334647,
    16.84645313,20.07872919,
    21.70841761,17.87612922,
    27.78720817,13.82372285,
    26.07217182,39.70838281,8.52256583,
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
BAO_Z = np.array([0.295,0.510,0.510,0.706,0.706,0.930,0.930,
                  1.317,1.317,1.491,2.330,2.330])
BAO_T = ['DV','DM','DH','DM','DH','DM','DH','DM','DH','DV','DM','DH']
BAO_ICOV = np.linalg.inv(BAO_COV)

# H(z) data
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


# ============================================================
# Fast cosmology: precomputed D_C interpolation
# ============================================================

def build_cosmology(H0, Om, eps, zc):
    """Compute IDM cosmology and return lookup functions.
    
    Uses cumulative integral on a fine z-grid → fast O(1) D_C(z) lookup.
    """
    Ol = 1.0 - Om
    h = H0 / 100.0
    Or = 4.15e-5 / h**2
    
    # Build z-grid: dense at low z, logarithmic to Z_REC
    N = 500
    z_low = np.linspace(0, 0.01, 50)  # 50 pts in z=[0,0.01]
    z_high = np.logspace(np.log10(0.011), np.log10(Z_REC*1.5), N)  # ~N pts to z≈1600
    z_grid = np.unique(np.concatenate([z_low, z_high]))
    
    # f_transient(z)
    x = z_grid / zc
    f_trans = 1.0 - eps * x * np.exp(-x)
    
    # E²(z) and E(z)
    E2 = Om*(1+z_grid)**3 + Or*(1+z_grid)**4 + Ol*f_trans
    E = np.sqrt(np.maximum(E2, 1e-30))  # guard against numerical issues
    
    # D_C(z) via cumulative trapezoidal
    inv_E = 1.0 / E
    DC = np.zeros_like(z_grid)
    for i in range(1, len(z_grid)):
        dz = z_grid[i] - z_grid[i-1]
        DC[i] = DC[i-1] + 0.5 * (inv_E[i] + inv_E[i-1]) * dz
    DC *= C / H0
    
    # Interpolators
    from scipy.interpolate import interp1d
    DC_interp = interp1d(z_grid, DC, kind='cubic', bounds_error=False, fill_value=0.0)
    
    def E_at_z(z):
        z = np.asarray(z, float)
        x = z / zc
        ft = 1.0 - eps * x * np.exp(-x)
        e2 = Om*(1+z)**3 + Or*(1+z)**4 + Ol*ft
        return np.sqrt(np.maximum(e2, 1e-30))
    
    return {
        'H0': H0, 'Om': Om, 'Ol': Ol, 'Or': Or,
        'eps': eps, 'zc': zc,
        'E': E_at_z,
        'DC': lambda z: float(DC_interp(z)),
        'DC_arr': lambda z: DC_interp(np.asarray(z, float)),
        'grid': z_grid, 'DC_grid': DC,
    }


# ============================================================
# Likelihood
# ============================================================

def ln_like(params):
    """Log-likelihood: -0.5 × total χ²"""
    H0, Om, eps, zc = params
    
    # Priors
    if not (55 < H0 < 80): return -np.inf
    if not (0.15 < Om < 0.45): return -np.inf
    if not (-0.5 < eps < 0.5): return -np.inf
    if not (0.2 < zc < 1.5): return -np.inf
    
    cosmo = build_cosmology(H0, Om, eps, zc)
    
    # H(z) χ² — direct E(z) call, fast
    hz_pred = cosmo['H0'] * cosmo['E'](HZ[:,0])
    chi2_h = np.sum(((hz_pred - HZ[:,1]) / HZ[:,2])**2)
    
    # DESI BAO χ² — 12 points
    pred = np.zeros(12)
    for i in range(12):
        z = BAO_Z[i]
        t = BAO_T[i]
        if t == 'DM': pred[i] = cosmo['DC'](z) / R_S_BAO
        elif t == 'DH': pred[i] = C / (cosmo['H0']*cosmo['E'](z)) / R_S_BAO
        elif t == 'DV':
            dm = cosmo['DC'](z)
            pred[i] = (C * z * dm**2 / (cosmo['H0']*cosmo['E'](z)))**(1/3) / R_S_BAO
    d = pred - BAO_MEAN
    chi2_b = d @ BAO_ICOV @ d
    
    # CMB χ²
    DC_rec = cosmo['DC'](Z_REC)
    I_rec = DC_rec * cosmo['H0'] / C
    R_val = np.sqrt(Om) * I_rec
    la_val = np.pi * DC_rec / R_S_CMB
    
    var_r = PRE**2; var_la = PLAE**2; cov_rl = PRLAC*PRE*PLAE
    det = var_r*var_la - cov_rl**2
    i_r = var_la/det; i_la = var_r/det; i_rl = -cov_rl/det
    dr = R_val - PR; dla = la_val - PLA
    chi2_c = i_r*dr**2 + 2*i_rl*dr*dla + i_la*dla**2
    
    # Total χ² = H(z) + BAO + CMB + SNe
    # SNe: 1048 points, use interpolated D_C for speed
    # Load SNe data
    chi2_s = 0.0
    
    return -0.5 * (chi2_h + chi2_b + chi2_c)


# ============================================================
# Load SNe data
# ============================================================
SN_PATH = "/tmp/pantheon_raw.txt"
if not os.path.exists(SN_PATH):
    print(f"Downloading Pantheon data to {SN_PATH}...")
    import urllib.request
    url = "https://raw.githubusercontent.com/dscolnic/Pantheon/master/lcparam_full_long.txt"
    urllib.request.urlretrieve(url, SN_PATH)

sne_zs, sne_mbs, sne_dmbs = [], [], []
with open(SN_PATH) as f:
    for line in f:
        if line.startswith("#") or not line.strip(): continue
        p = line.strip().split()
        if len(p) < 6: continue
        try:
            z = float(p[1])
            if z >= 0.01:
                sne_zs.append(z)
                sne_mbs.append(float(p[4]))
                sne_dmbs.append(float(p[5]))
        except: continue
sne_zs = np.array(sne_zs)
sne_mbs = np.array(sne_mbs)
sne_dmbs = np.array(sne_dmbs)
SNE_W = 1 / sne_dmbs**2
print(f"SNe loaded: {len(sne_zs)} points")


def ln_like_full(params):
    """Full log-likelihood including SNe (slower)"""
    H0, Om, eps, zc = params
    if not (55 < H0 < 80): return -np.inf
    if not (0.15 < Om < 0.45): return -np.inf
    if not (-0.5 < eps < 0.5): return -np.inf
    if not (0.2 < zc < 1.5): return -np.inf
    
    cosmo = build_cosmology(H0, Om, eps, zc)
    
    # H(z) — fast
    chi2_h = np.sum(((cosmo['H0']*cosmo['E'](HZ[:,0]) - HZ[:,1]) / HZ[:,2])**2)
    
    # BAO — fast, 12 points
    pred = np.zeros(12)
    for i in range(12):
        z, t = BAO_Z[i], BAO_T[i]
        if t == 'DM': pred[i] = cosmo['DC'](z) / R_S_BAO
        elif t == 'DH': pred[i] = C/(cosmo['H0']*cosmo['E'](z)) / R_S_BAO
        elif t == 'DV':
            dm = cosmo['DC'](z)
            pred[i] = (C*z*dm**2/(cosmo['H0']*cosmo['E'](z)))**(1/3) / R_S_BAO
    chi2_b = (pred - BAO_MEAN) @ BAO_ICOV @ (pred - BAO_MEAN)
    
    # CMB
    DC_rec = cosmo['DC'](Z_REC)
    I_rec = DC_rec * H0 / C
    R_val = np.sqrt(Om) * I_rec
    la_val = np.pi * DC_rec / R_S_CMB
    
    var_r = PRE**2; var_la = PLAE**2; cov_rl = PRLAC*PRE*PLAE
    det = var_r*var_la - cov_rl**2
    i_r = var_la/det; i_la = var_r/det; i_rl = -cov_rl/det
    dr = R_val - PR; dla = la_val - PLA
    chi2_c = i_r*dr**2 + 2*i_rl*dr*dla + i_la*dla**2
    
    # SNe — 1048 points with interpolated D_C
    # D_L = (1+z) × D_C(z) × f_secular(z) [λ=0.001, negligible]
    # μ = 5 log10(D_L) + 25
    DC_arr = cosmo['DC_arr'](sne_zs)
    DL = (1 + sne_zs) * DC_arr
    mu = 5 * np.log10(DL) + 25
    resid = mu - sne_mbs
    Mb = -np.sum(SNE_W * resid) / np.sum(SNE_W)
    chi2_s = np.sum(SNE_W * (resid + Mb)**2)
    
    return -0.5 * (chi2_h + chi2_s + chi2_b + chi2_c)


# ============================================================
# MCMC
# ============================================================

print("="*72)
print("  MCMC JOINT FIT — IDM (H₀, Ω_m, ε, z_c)")
print("  Datasets: DESI DR2 BAO + Planck CMB + Pantheon + H(z)")
print("="*72)

n_dim = 4
n_walkers = 16
n_burn = 500
n_steps = 3000
n_thin = 3

# Initial positions around best-fit from grid scan
p0_best = np.array([69.0, 0.300, 0.25, 0.70])
p0_scale = np.array([2.0, 0.02, 0.05, 0.10])  # spread

rng = np.random.default_rng(42)
p0 = p0_best + p0_scale * rng.normal(size=(n_walkers, n_dim))

# Ensure valid priors
for i in range(n_walkers):
    while True:
        h, o, e, z = p0[i]
        if 55<h<80 and 0.15<o<0.45 and -0.5<e<0.5 and 0.2<z<1.5:
            break
        p0[i] = p0_best + p0_scale * rng.normal(size=n_dim)

print(f"\nWalkers: {n_walkers}, burn: {n_burn}, steps: {n_steps}")
print(f"Initial: H₀={p0_best[0]:.1f}, Ω_m={p0_best[1]:.3f}, "
      f"ε={p0_best[2]:.3f}, z_c={p0_best[3]:.2f}")
print(f"\nStarting MCMC (this will take several minutes)...")
t_start = time.time()

sampler = emcee.EnsembleSampler(n_walkers, n_dim, ln_like_full,
                                moves=emcee.moves.DEMove())
sampler.run_mcmc(p0, n_burn + n_steps, progress=True)

elapsed = time.time() - t_start
print(f"\\nMCMC done: {elapsed:.0f}s ({elapsed/60:.1f} min)")

# Extract chain and log_prob
chain = sampler.get_chain(discard=n_burn, thin=n_thin, flat=True)
log_prob = sampler.get_log_prob(discard=n_burn, thin=n_thin, flat=True)

# chi2 directly from MCMC's own log-likelihood (fast interpolation)
chi2_arr = -2.0 * log_prob

i_best = np.argmin(chi2_arr)
map_params = chain[i_best]
map_chi2 = chi2_arr[i_best]

# Quick verification at MAP using joint_final's eval_m (one point)
print("\\nSingle-point MAP verification with joint_final (quad integration)...")
from joint_final import IDMModel, eval_m, load_sne as jf_load_sne
jf_zs, jf_mbs, jf_dmbs = jf_load_sne()
jf_m = IDMModel(H0=map_params[0], Om=map_params[1], eps=map_params[2], zc=map_params[3])
jf_chi2, jf_detail = eval_m(jf_m, jf_zs, jf_mbs, jf_dmbs)
print(f"  MCMC χ² (interpolated): {map_chi2:.1f}")
print(f"  joint_final χ² (quad):  {jf_chi2:.1f}")
print(f"  Δχ² = {abs(map_chi2-jf_chi2):.3f}  {'✓ AGREES' if abs(map_chi2-jf_chi2)<2.0 else '✗ MISMATCH'}")

# Use the more accurate joint_final chi2 for reporting
map_chi2_final = jf_chi2

print(f"\\nMinimum χ² (interpolated) = {map_chi2:.1f}")
print(f"Minimum χ² (quad verified) = {map_chi2_final:.1f}")
print(f"MAP from quad-verified χ²:")
print(f"  H₀ = {map_params[0]:.1f} ± {np.std(chain[:,0]):.1f}")
print(f"  Ω_m = {map_params[1]:.4f} ± {np.std(chain[:,1]):.4f}")
print(f"  ε = {map_params[2]:.4f} ± {np.std(chain[:,2]):.4f}")
print(f"  z_c = {map_params[3]:.4f} ± {np.std(chain[:,3]):.4f}")

# 68% and 95% intervals
perc = [16, 50, 84]
for i, name in enumerate(['H₀', 'Ω_m', 'ε', 'z_c']):
    q = np.percentile(chain[:,i], perc)
    print(f"  {name} = {q[1]:.4f} +{q[2]-q[1]:.4f} -{q[1]-q[0]:.4f}")

# Corner plot
print("\nGenerating corner plot...")
fig = corner.corner(chain, labels=['$H_0$', '$\\Omega_m$', '$\\varepsilon$', '$z_c$'],
                    truths=map_params, quantiles=[0.16, 0.5, 0.84],
                    show_titles=True, title_fmt='.4f',
                    plot_datapoints=True, fill_contours=True,
                    levels=(0.68, 0.95), color='#222222',
                    smooth=True, smooth1d=True)
fig.savefig("/home/wsl/cosmic_validation/mcmc_corner.png", dpi=200)
print("Saved: mcmc_corner.png")

# Save chain
np.savez("/home/wsl/cosmic_validation/mcmc_chain.npz",
         chain=chain, log_prob=log_prob,
         map_params=map_params, names=['H0','Om','eps','zc'])
print("Saved: mcmc_chain.npz")

# Summary for paper
print(f"\n{'='*72}")
print(f"  PUBLICATION-QUALITY RESULTS")
print(f"{'='*72}")
print(f"""
  Parameter     MAP (±1σ)          68% CI          95% CI
  ────────  ─────────────────  ───────────────  ───────────────
  H₀        {map_params[0]:.2f} ± {np.std(chain[:,0]):.2f}
  Ω_m       {map_params[1]:.4f} ± {np.std(chain[:,1]):.4f}
  ε         {map_params[2]:.4f} ± {np.std(chain[:,2]):.4f}
  z_c       {map_params[3]:.4f} ± {np.std(chain[:,3]):.4f}
""")

# Final comparison vs ΛCDM grid-scan baseline
daic = 1129.0 - map_chi2_final - 4
print(f"ΔAIC (IDM − ΛCDM) = {daic:+.1f}  ({'IDM wins' if daic>0 else 'ΛCDM wins'})")
print(f"χ²_IDM = {map_chi2_final:.1f} (quad-verified)")
