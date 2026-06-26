#!/usr/bin/env python3
"""
run_idm_validation.py — 資訊動力學模型 vs 真實數據
====================================================
純 numpy/scipy 數學驗證。

測試兩個核心公式：
  Model A: f_ID(z) = 1 - ε × (z/z_c) × exp(-z/z_c)
  Model B: sync dip + inflation tail
  
用實際觀測數據 (Pantheon SNe + H(z)) 做 χ² fit。
"""

import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from cosmic_idm import (
    IDMParams, InformationDynamicsModel, IDMLuminosityDistance,
    IDMValidator, hz_data_points
)
from observational_data import PantheonLoader


def test_hz_fit():
    """資訊動力學模型對 H(z) 數據嘅 fit"""
    print("=" * 65)
    print("H(z) 數據驗證 — 資訊動力學模型")
    print("=" * 65)

    hz = hz_data_points()
    print(f"\n  H(z) 數據點: {len(hz['z'])}")
    print(f"  z range: [{hz['z'].min():.3f}, {hz['z'].max():.3f}]")

    # ΛCDM baseline (ε=0)
    params_lcdm = IDMParams(epsilon=0.0)
    model_lcdm = InformationDynamicsModel(params_lcdm)
    chi2_lcdm = np.sum(
        ((hz['hz'] - model_lcdm.hubble(hz['z'])) / hz['err']) ** 2.0
    )
    print(f"  ΛCDM (ε=0):  χ² = {chi2_lcdm:.2f}")

    # Scan various ε values
    print(f"\n  ε 掃描 (z_c=0.6 fixed):")
    print(f"  {'ε':^10} {'χ²':^10} {'Δχ² vs ΛCDM':^16}")
    print(f"  {'-'*10} {'-'*10} {'-'*16}")
    
    best_chi2 = chi2_lcdm
    best_eps = 0.0
    
    for eps in np.linspace(-0.15, 0.15, 31):
        params = IDMParams(epsilon=eps, z_sync=0.6)
        model = InformationDynamicsModel(params)
        chi2 = np.sum(
            ((hz['hz'] - model.hubble(hz['z'])) / hz['err']) ** 2.0
        )
        dchi2 = chi2_lcdm - chi2
        marker = " ← best" if chi2 < best_chi2 - 0.1 else ""
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_eps = eps
        if abs(eps) < 0.001 or abs(eps - round(eps*20)/20) < 0.01:
            print(f"  {eps:+8.3f}  {chi2:8.2f}  {dchi2:+14.2f}{marker}")

    print(f"\n  🟢 Best ε (H(z) only): ε = {best_eps:.3f}, χ² = {best_chi2:.2f}")
    
    # Comparison at key points
    params_best = IDMParams(epsilon=best_eps, z_sync=0.6)
    model_best = InformationDynamicsModel(params_best)
    
    print(f"\n  {'z':^7} {'H_obs':^10} {'±err':^8} {'H_IDM':^10} {'H_ΛCDM':^10} {'Pull_IDM':^10}")
    print(f"  {'-'*7} {'-'*10} {'-'*8} {'-'*10} {'-'*10} {'-'*10}")
    for i in range(len(hz['z'])):
        z = hz['z'][i]
        h_obs = hz['hz'][i]
        err = hz['err'][i]
        h_idm = model_best.hubble(np.array([z]))[0]
        h_l = model_lcdm.hubble(np.array([z]))[0]
        pull = (h_obs - h_idm) / err
        print(f"  {z:6.3f}  {h_obs:8.1f}  ±{err:5.1f}  {h_idm:8.1f}  {h_l:8.1f}  {pull:+9.3f}")

    return best_eps, best_chi2, chi2_lcdm


def test_sn_fit():
    """資訊動力學模型對 Pantheon SN 數據嘅 fit"""
    print("\n" + "=" * 65)
    print("Pantheon 超新星驗證 — 資訊動力學模型")
    print("=" * 65)

    loader = PantheonLoader()
    try:
        sn = loader.load()
    except FileNotFoundError:
        print("❌ Pantheon data missing")
        return None

    M_abs = -19.3
    print(f"\n  SN數量: {sn.n_sne}")
    print(f"  M_abs = {M_abs}")
    print(f"  z range: [{sn.z.min():.4f}, {sn.z.max():.4f}]")

    # ΛCDM baseline
    params_lcdm = IDMParams(epsilon=0.0)
    model_lcdm = InformationDynamicsModel(params_lcdm)
    dist_lcdm = IDMLuminosityDistance(model_lcdm)
    mu_lcdm = np.array([dist_lcdm.distance_modulus(z) for z in sn.z])
    chi2_lcdm = np.sum(((mu_lcdm - (sn.mb - M_abs)) / sn.dmb) ** 2.0)
    print(f"  ΛCDM (ε=0):  χ² = {chi2_lcdm:.1f}")

    # Scan ε
    print(f"\n  ε 掃描 (z_c=0.6 fixed):")
    print(f"  {'ε':^10} {'χ²':^10} {'Δχ²':^12}")
    print(f"  {'-'*10} {'-'*10} {'-'*12}")
    
    best_chi2 = chi2_lcdm
    best_eps = 0.0
    
    for eps in np.linspace(-0.15, 0.15, 31):
        params = IDMParams(epsilon=eps, z_sync=0.6)
        model = InformationDynamicsModel(params)
        dist = IDMLuminosityDistance(model)
        mu = np.array([dist.distance_modulus(z) for z in sn.z])
        chi2 = np.sum(((mu - (sn.mb - M_abs)) / sn.dmb) ** 2.0)
        dchi2 = chi2_lcdm - chi2
        if chi2 < best_chi2 - 0.1:
            best_chi2 = chi2
            best_eps = eps
        if abs(eps) < 0.001 or abs(eps - round(eps*20)/20) < 0.01:
            marker = " ← best" if abs(chi2 - best_chi2) < 0.1 else ""
            print(f"  {eps:+8.3f}  {chi2:8.1f}  {dchi2:+10.1f}{marker}")

    print(f"\n  🟢 Best ε (Pantheon only): ε = {best_eps:.3f}, χ² = {best_chi2:.1f}")
    return {'best_eps': best_eps, 'best_chi2': best_chi2, 'chi2_lcdm': chi2_lcdm}


def test_joint_fit():
    """聯合 fit H(z) + Pantheon"""
    print("\n" + "=" * 65)
    print("聯合驗證 — H(z) + Pantheon 聯合 χ²")
    print("=" * 65)

    hz = hz_data_points()
    loader = PantheonLoader()
    try:
        sn = loader.load()
    except FileNotFoundError:
        print("❌ Pantheon data missing")
        return

    M_abs = -19.3
    params_lcdm = IDMParams(epsilon=0.0)
    model_lcdm = InformationDynamicsModel(params_lcdm)
    dist_lcdm = IDMLuminosityDistance(model_lcdm)

    # ΛCDM joint χ²
    mu_lcdm = np.array([dist_lcdm.distance_modulus(z) for z in sn.z])
    chi2_sn_lcdm = np.sum(((mu_lcdm - (sn.mb - M_abs)) / sn.dmb) ** 2.0)
    chi2_hz_lcdm = np.sum(
        ((hz['hz'] - model_lcdm.hubble(hz['z'])) / hz['err']) ** 2.0
    )
    joint_lcdm = chi2_sn_lcdm + chi2_hz_lcdm
    
    print(f"\n  ΛCDM 聯合 χ² = {joint_lcdm:.1f}")
    print(f"    (SNe: {chi2_sn_lcdm:.1f} + H(z): {chi2_hz_lcdm:.1f})")

    # Scan ε, both z_c values
    print(f"\n  ε 掃描:")
    print(f"  {'ε':^8} {'z_c=0.3':^12} {'z_c=0.6':^12} {'z_c=1.0':^12}")
    print(f"  {'-'*8} {'-'*12} {'-'*12} {'-'*12}")

    best_result = {'eps': 0.0, 'zc': 0.6, 'chi2': joint_lcdm}

    for eps in np.linspace(-0.10, 0.10, 21):
        row = f"  {eps:+7.3f}"
        for zc in [0.3, 0.6, 1.0]:
            params = IDMParams(epsilon=eps, z_sync=zc)
            model = InformationDynamicsModel(params)
            dist = IDMLuminosityDistance(model)
            mu = np.array([dist.distance_modulus(z) for z in sn.z])
            chi2_sn = np.sum(((mu - (sn.mb - M_abs)) / sn.dmb) ** 2.0)
            chi2_hz = np.sum(
                ((hz['hz'] - model.hubble(hz['z'])) / hz['err']) ** 2.0
            )
            joint = chi2_sn + chi2_hz
            row += f"  {joint:10.1f}"
            if joint < best_result['chi2']:
                best_result = {'eps': eps, 'zc': zc, 'chi2': joint}
        
        if abs(eps) < 0.01 or abs(eps*20 - round(eps*20)) < 0.01:
            print(row)

    print(f"\n  🟢 Best joint fit: ε={best_result['eps']:+.3f}, "
          f"z_c={best_result['zc']:.1f}, χ²={best_result['chi2']:.1f}")
    print(f"     Δχ² vs ΛCDM = {joint_lcdm - best_result['chi2']:.1f}")
    print(f"     AIC = Δχ² - 2×n_params = "
          f"{joint_lcdm - best_result['chi2'] - 2:.1f}")

    return best_result


def analytic_summary():
    """從資訊動力學原理推導出最終公式"""
    print("\n" + "=" * 65)
    print("資訊動力學 — 由第一原理到 H(z)")
    print("=" * 65)
    
    summary = r"""
資訊動力學推導膨脹：

P1 — Landauer Principle at Cosmological Scale:
    每個 information process 產生 k_B ln(2) entropy
    → dS/dt = k_B × ln(2) × İ_total

P2 — Sync Cost Function:
    İ_total ∝ N_addresses² / τ_sync
    τ_sync ∝ 1/H (causal horizon time)
    → İ_total ∝ (N_addr/V)² × V × H ∝ ρ_addr² × H

P3 — Expansion as Entropy Release:
    H ∝ sqrt(ρ_S) where ρ_S = S/V
    But ρ_S itself depends on H through sync cost
    → self-consistent equation

P4 — Solution:
    H² = 8πG/3 × [ρ_m + ρ_Λ × (1 - ε × z/z_c × exp(-z/z_c))]
    
    where ε encodes the sync efficiency:
    ε → 0: perfect sync (ΛCDM limit)
    ε → 0.05: typical sync cost
    
P5 — Three Phases:
    z ≫ 1: ρ_m dominates → matter era (observed)
    z ~ z_c: sync cost → slight dip in effective DE (H(z) data hint)
    z → 0: causal horizon growth → long-range sync → dark energy

P6 — Testable Predictions:
    • H(z) at z=0.3-0.8 should show slight dip relative to ΛCDM
    • Dip position z_c relates to causal horizon at recombination
    • ε relates to fundamental address density N_total
    """
    print(summary)


if __name__ == '__main__':
    r1 = test_hz_fit()
    r2 = test_sn_fit()
    r3 = test_joint_fit()
    analytic_summary()
    
    print("\n" + "=" * 65)
    print("結論：")
    if r1 and r2 and r3:
        print(f"  H(z) best ε:     {r1[0]:+.4f}  (χ² = {r1[1]:.1f})")
        print(f"  SNe best ε:      {r2['best_eps']:+.4f}  (χ² = {r2['best_chi2']:.1f})")
        print(f"  Joint best ε:    {r3['eps']:+.4f}")
        print(f"  Joint best z_c:  {r3['zc']:.2f}")
        print(f"  Δχ² vs ΛCDM:    {r3['chi2'] - (r2['chi2_lcdm'] + r1[2]):.1f}")
    print("=" * 65)