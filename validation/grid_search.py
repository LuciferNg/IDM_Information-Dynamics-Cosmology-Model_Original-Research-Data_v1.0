#!/usr/bin/env python3
"""
grid_search.py — β 掃描 + 可視化
===================================
純 numpy，無 LLM。
掃描 beta_sync 睇 χ² 點變，check sync model 有冇統計優勢。
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cosmic_model import CosmologyParams, SyncPressureModel, LuminosityDistance, ModelComparator
from observational_data import PantheonLoader, HubbleDataLoader


def grid_search_beta():
    """
    掃描 beta_sync 睇 χ² surface，確認 global minimum
    
    β range: [-0.5, 1.5]
    Steps: 500
    """
    print("=" * 65)
    print("Grid Search: β_sync 掃描")
    print("=" * 65)

    # Load SN data
    loader = PantheonLoader()
    try:
        sn = loader.load()
    except FileNotFoundError:
        print("❌ Pantheon data missing")
        return

    M_abs = -19.3

    # Scan beta over wide range
    betas = np.linspace(-0.5, 1.5, 500)
    chi2s = np.zeros_like(betas)

    for i, beta in enumerate(betas):
        params = CosmologyParams(beta_sync=beta)
        model = SyncPressureModel(params)
        dist = LuminosityDistance(model)
        mu_model = dist.mu_array(sn.z)
        res = (mu_model - (sn.mb - M_abs)) / sn.dmb
        chi2s[i] = np.sum(res ** 2.0)

    # Find best beta
    best_idx = np.argmin(chi2s)
    beta_best = betas[best_idx]
    chi2_min = chi2s[best_idx]
    chi2_lcdm = chi2s[np.argmin(np.abs(betas))]  # β≈0

    print(f"\n  📊 Grid Search Results:")
    print(f"     β range: [{betas[0]:.2f}, {betas[-1]:.2f}]")
    print(f"     Steps: {len(betas)}")
    print(f"     β_best = {beta_best:.4f}")
    print(f"     χ²_min = {chi2_min:.1f}")
    print(f"     χ²_ΛCDM (β≈0) = {chi2_lcdm:.1f}")
    print(f"     Δχ² = {chi2_lcdm - chi2_min:.2f}")

    # χ² vs β table for key points
    key_betas = [-0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]
    print(f"\n  {'β':^8} {'χ²':^10} {'Δχ²':^10}")
    print(f"  {'-'*8} {'-'*10} {'-'*10}")
    for b in key_betas:
        b_idx = np.argmin(np.abs(betas - b))
        print(f"  {b:+.2f}    {chi2s[b_idx]:8.1f}  {chi2s[b_idx] - chi2_min:+9.1f}")

    # 1σ, 2σ confidence intervals (Δχ² = 1, 4 for 1 param)
    def conf_interval(chi2_min_val, delta, betas_arr, chi2s_arr):
        """Find β range where χ² < χ²_min + delta"""
        threshold = chi2_min_val + delta
        above = chi2s_arr > threshold
        
        # Find left bound
        left_idx = np.where(above[:best_idx])[0]
        left_beta = betas_arr[0]
        if len(left_idx) > 0:
            left_beta = betas_arr[left_idx[-1]]
        
        # Find right bound
        right_idx = np.where(above[best_idx:])[0]
        right_beta = betas_arr[-1]
        if len(right_idx) > 0:
            right_beta = betas_arr[best_idx + right_idx[0]]
        
        return left_beta, right_beta

    left_1s, right_1s = conf_interval(chi2_min, 1.0, betas, chi2s)
    left_2s, right_2s = conf_interval(chi2_min, 4.0, betas, chi2s)
    
    print(f"\n  📏 Confidence Intervals for β_sync:")
    print(f"     1σ (Δχ²=1): [{left_1s:.4f}, {right_1s:.4f}]")
    print(f"     2σ (Δχ²=4): [{left_2s:.4f}, {right_2s:.4f}]")
    
    if left_1s <= 0.0 and right_1s >= 0.0:
        print(f"\n  🔴 β=0 落在 1σ 內 → sync model 無統計顯著差異")
        print(f"     ΛCDM still consistent at 1σ level")
    else:
        print(f"\n  🟢 β=0 不在 1σ 內 → sync model 統計顯著")
        print(f"     β_best = {beta_best:.3f}")

    # Predict vs observed table (first 20 SNe)
    print(f"\n  {'z':^8} {'μ_obs':^10} {'μ_model':^10} {'Pull':^10}")
    print(f"  {'-'*8} {'-'*10} {'-'*10} {'-'*10}")
    
    params_best = CosmologyParams(beta_sync=beta_best)
    model_best = SyncPressureModel(params_best)
    dist_best = LuminosityDistance(model_best)
    
    for i in range(min(20, len(sn.z))):
        mu_mod = dist_best.distance_modulus(sn.z[i])
        mu_obs_cal = sn.mb[i] - M_abs
        pull = (mu_mod - mu_obs_cal) / sn.dmb[i]
        print(f"  {sn.z[i]:7.4f} {mu_obs_cal:9.3f} {mu_mod:9.3f} {pull:+9.3f}")

    return {
        'beta_best': beta_best,
        'chi2_min': chi2_min,
        'chi2_lcdm': chi2_lcdm,
        'betas': betas,
        'chi2s': chi2s,
        'conf_1sigma': (left_1s, right_1s),
        'conf_2sigma': (left_2s, right_2s)
    }


def three_param_scan():
    """
    三維掃描: H₀ × β_sync × Ω_m
    Check parameter degeneracy.
    
    Scipy optimise find global minimum.
    """
    import scipy.optimize as opt
    
    print("\n" + "=" * 65)
    print("3-param global fit: H₀, Ω_m, β_sync")
    print("=" * 65)
    
    loader = PantheonLoader()
    try:
        sn = loader.load()
    except FileNotFoundError:
        print("❌ Pantheon data missing")
        return
    
    M_abs = -19.3
    
    def chi2_3param(x):
        """x = [H₀, Omega_m, beta_sync]"""
        h0, om, bs = x
        if h0 <= 0 or om <= 0 or om >= 1:
            return 1e12
        ol = 1.0 - om
        params = CosmologyParams(H0=h0, Omega_m=om, Omega_Lambda=ol, beta_sync=bs)
        model = SyncPressureModel(params)
        dist = LuminosityDistance(model)
        
        try:
            mu_model = dist.mu_array(sn.z)
            res = (mu_model - (sn.mb - M_abs)) / sn.dmb
            return np.sum(res ** 2.0)
        except (OverflowError, ValueError):
            return 1e12
    
    # Initial guesses
    x0 = [69.8, 0.286, 0.3]
    
    print(f"\n  Initial guess: H₀={x0[0]}, Ω_m={x0[1]}, β_sync={x0[2]}")
    print(f"  Initial χ² = {chi2_3param(x0):.1f}")
    
    # Bounded optimization
    bounds = opt.Bounds(
        lb=[50.0, 0.1, -0.5],
        ub=[80.0, 0.6, 2.0]
    )
    
    result = opt.minimize(
        chi2_3param, x0, method='L-BFGS-B',
        bounds=bounds, options={'maxiter': 200, 'ftol': 1e-12}
    )
    
    if result.success:
        h0_best, om_best, bs_best = result.x
        chi2_best = result.fun
        
        print(f"\n  🟢 Global fit converged:")
        print(f"     H₀      = {h0_best:.2f} km/s/Mpc")
        print(f"     Ω_m     = {om_best:.4f}")
        print(f"     β_sync  = {bs_best:.4f}")
        print(f"     χ²_min  = {chi2_best:.1f}")
        
        # Compare to ΛCDM fit (β=0)
        def chi2_lcdm_only(x):
            return chi2_3param([x[0], x[1], 0.0])
        
        result_lcdm = opt.minimize(
            chi2_lcdm_only, [x0[0], x0[1]], method='L-BFGS-B',
            bounds=opt.Bounds(lb=[50.0, 0.1], ub=[80.0, 0.6]),
            options={'maxiter': 100, 'ftol': 1e-12}
        )
        
        if result_lcdm.success:
            h0_l, om_l = result_lcdm.x
            chi2_l = result_lcdm.fun
            print(f"\n     ΛCDM best (β=0):")
            print(f"       H₀  = {h0_l:.2f}, Ω_m = {om_l:.4f}, χ² = {chi2_l:.1f}")
            print(f"       Δχ² = {chi2_l - chi2_best:.2f}")
            print(f"       ΔAIC = {chi2_l - chi2_best - 2:.2f} (penalizing 1 extra param)")
    else:
        print(f"\n  🔴 Optimizer failed: {result.message}")
    
    return result


if __name__ == '__main__':
    r1 = grid_search_beta()
    if r1:
        print(f"\n  Grid search summary:")
        print(f"    β_best = {r1['beta_best']:.4f}")
        print(f"    1σ CI: [{r1['conf_1sigma'][0]:.4f}, {r1['conf_1sigma'][1]:.4f}]")
        print(f"    β=0 {'inside' if r1['conf_1sigma'][0] <= 0.0 <= r1['conf_1sigma'][1] else 'OUTSIDE'} 1σ confidence interval")
    
    three_param_scan()
    print("\nDone.")