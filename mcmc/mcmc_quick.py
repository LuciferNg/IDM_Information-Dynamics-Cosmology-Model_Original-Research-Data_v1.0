"""
mcmc_quick.py — Fast MCMC for IDM
Uses joint_final.py's χ² computation (verified) for MCMC sampling.
Run: uv run python3 mcmc_quick.py
"""

import numpy as np
import emcee, corner, time, os, sys

# Import verified grid-scan code
sys.path.insert(0, '/home/wsl/cosmic_validation')
from joint_final import IDMModel, eval_m, load_sne

zs, mbs, dmbs = load_sne()
print(f"SNe loaded: {len(zs)} points")

def ln_posterior(params):
    H0, Om, eps, zc = params
    if not (55<H0<80): return -np.inf
    if not (0.15<Om<0.45): return -np.inf
    if not (-0.5<eps<0.5): return -np.inf
    if not (0.15<zc<1.5): return -np.inf
    m = IDMModel(H0=H0, Om=Om, eps=eps, zc=zc)
    t, _ = eval_m(m, zs, mbs, dmbs)
    return -0.5 * t

n_dim = 4; n_walkers = 24; n_burn = 1000; n_steps = 5000

p0_best = np.array([69.0, 0.300, 0.25, 0.70])
p0 = p0_best + np.array([2.0, 0.02, 0.05, 0.10]) * np.random.randn(n_walkers, n_dim)

print(f"\nWalkers: {n_walkers}, burn: {n_burn}, steps: {n_steps}")
t0 = time.time()

sampler = emcee.EnsembleSampler(n_walkers, n_dim, ln_posterior)
sampler.run_mcmc(p0, n_burn + n_steps, progress=True)
print(f"Elapsed: {time.time()-t0:.0f}s")

# Extract
chain = sampler.get_chain(discard=n_burn, thin=5, flat=True)
ln_probs = sampler.get_log_prob(discard=n_burn, thin=5, flat=True)
chi2_arr = -2 * ln_probs

i_best = np.argmin(chi2_arr)
map_p = chain[i_best]
map_chi2 = chi2_arr[i_best]

# Verify: compute chi2 via joint_final at MAP
m = IDMModel(H0=map_p[0], Om=map_p[1], eps=map_p[2], zc=map_p[3])
t_v, _ = eval_m(m, zs, mbs, dmbs)
print(f"\nMAP verification:")
print(f"  MCMC χ² = {map_chi2:.1f}  |  joint_final χ² = {t_v:.1f}  |  diff = {abs(map_chi2-t_v):.3f}")

# Parameter constraints
print(f"\n{'='*60}")
print(f"  RESULTS")
print(f"{'='*60}")
print(f"  MAP params: H₀={map_p[0]:.1f}, Ω_m={map_p[1]:.4f}, ε={map_p[2]:.4f}, z_c={map_p[3]:.4f}")
print(f"  MAP χ² = {map_chi2:.1f}")

names = ['H₀', 'Ω_m', 'ε', 'z_c']
for i, name in enumerate(names):
    vals = chain[:, i]
    p16, p50, p84 = np.percentile(vals, [16, 50, 84])
    print(f"  {name:4s} = {p50:.4f} +{p84-p50:.4f} -{p50-p16:.4f}")

# AIC vs ΛCDM
lcdm_chi2 = 1129.0  # from grid scan
daic = lcdm_chi2 - map_chi2 - 4  # IDM has 2 more params
print(f"\n  ΔAIC vs ΛCDM = {daic:.1f}")

# Corner plot
fig = corner.corner(chain, labels=['$H_0$','$\\Omega_m$','$\\varepsilon$','$z_c$'],
                    quantiles=[0.16,0.5,0.84], show_titles=True, title_fmt='.3f',
                    levels=(0.68,0.95), color='#1a5276', smooth=True, bins=50)
fig.savefig('/home/wsl/cosmic_validation/mcmc_corner_final.png', dpi=300)
print(f"\nSaved: mcmc_corner_final.png")

np.savez('/home/wsl/cosmic_validation/mcmc_final.npz', chain=chain, chi2=chi2_arr, map_params=map_p)
print(f"Saved: mcmc_final.npz")