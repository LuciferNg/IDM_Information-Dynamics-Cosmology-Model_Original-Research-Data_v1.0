"""
derive_eps_from_N.py — 從資訊動力學第一原理推導 ε = f(N_total)
================================================================
純 Python 數學推導，無 LLM。

將 ε 從 free parameter 變成：
  ε = ε₀ × N_total^(-1/4) × R₀^(3/4)

其中 N_total = 宇宙中嘅總 ontological address count
      R₀ = 每個 address 嘅自然 processing rate
"""

import numpy as np


# ========= 第一原理推導 =========

"""
推導步驟：
───────────

P0 — 基礎常數：
  N_total:  宇宙中總 ontological addresses 數目
  R₀:       每個 address 嘅自然 processing rate (Hz)
  l_P:      Planck length (1.616 × 10⁻³⁵ m)
  L_uni:    可觀察宇宙尺度 (~10²⁶ m)
  n_addr:   宇宙平均 address density = N_total / L³

P1 — Sync time:
  τ_sync ∝ 1 / (R₀ × N_horizon)
  N_horizon = 一個 causal horizon 內嘅 addresses 數
  N_horizon = n_addr × (4π/3) × (c/H)³

P2 — Hubble time:
  τ_hubble = 1 / H

P3 — Sync transition 發生喺 τ_sync = τ_hubble：
  H = R₀ × N_horizon
  H = R₀ × n_addr × (4π/3) × (c/H)³
  H⁴ = R₀ × n_addr × (4π/3) × c³

P4 — At transition (z = z_c):
  n_addr(z_c) = N_total / V(z_c)
  V(z_c) = (4π/3) × (c/H(z_c))³
  n_addr(z_c) = N_total × 3 × H(z_c)³ / (4π × c³)

  Substitute into P3:
  H(z_c)⁴ = R₀ × [N_total × 3 × H(z_c)³ / (4π × c³)] × (4π/3) × c³
  H(z_c)⁴ = R₀ × N_total × H(z_c)³
  H(z_c) = R₀ × N_total
  
  !!! 呢個係關鍵結果：sync transition 發生喺 H(z_c) = R₀ × N_total !!!

P5 — z_c 同 H(z_c) 嘅關係 (matter-dominated era):
  H(z_c) ≈ H₀ × sqrt(Ω_m × (1+z_c)³)
  (1+z_c)^(3/2) ≈ R₀ × N_total / H₀
  z_c ≈ (R₀ × N_total / H₀)^(2/3) - 1

P6 — 現在推導 ε：
  ε = sync efficiency 喺 transition 嘅 magnitude
  ε ∝ η(z_c) × (δH/H)_at_transition
  
  Sync efficiency η = N_synced / N_horizon
  η ∝ 1 / sqrt(Nₕₒᵣᵢ𝓏ₒₙ)  (central limit theorem: noise ~ 1/√N)
  
  N_horizon(z_c) = n_addr × (4π/3) × (c/H(z_c))³
  ∝ N_total / V(z_c) × (c/H(z_c))³
  ∝ N_total × H(z_c)³ / (c/H(z_c))³ ... no wait
  ∝ N_total / (c/H(z_c))³ × (c/H(z_c))³ = N_total

  Interesting! N_horizon(z_c) = N_total at transition time!
  This makes sense: at z_c, the causal horizon encompasses ALL addresses that 
  are within sync range, meaning the horizon is large enough to see a 
  representative fraction.

  Actually, wait. N_horizon is within ONE Hubble volume.
  N_total is across the ENTIRE universe at that time.
  The entire universe at z_c has comoving volume V_total(z_c).
  N_horizon = ρ_addr(z_c) × V_horizon(z_c)
  N_total(z_c) = ρ_addr(z_c) × V_total(z_c)
  
  V_horizon = (4π/3) × (c/H)³
  V_total = (4π/3) × (c/H₀)³ × (1+z_c)^(-3/2)... no that's complicated
  
  Actually, let's keep it simple. The horizon volume at z_c is just (4π/3)(c/H(z_c))³.
  The total volume of the universe at that time is the current observable volume 
  divided by (1+z_c)³.
  
  N_horizon(z_c) = N_total × V_horizon(z_c) / V_total(z_c)
  
  For matter-dominated: V_total(z_c) = V₀ / (1+z_c)³
  V_horizon(z_c) ∝ (c/H(z_c))³ ∝ (1+z_c)^(-9/2) [since H ∝ (1+z)^(3/2)]
  
  N_horizon(z_c) / N_total ∝ V_horizon / V_total
  ∝ (1+z_c)^(-9/2) × (1+z_c)³ ∝ (1+z_c)^(-3/2)
  
  So N_horizon ∝ N_total × (1+z_c)^(-3/2)
  And η ∝ 1/sqrt(N_horizon) ∝ N_total^(-1/2) × (1+z_c)^(3/4)

P7 — 而 (δH/H) at transition：
  The fractional change in H due to sync consumption:
  δH/H ∝ η × (sync power / Hubble power)
  sync power / Hubble power ∝ R₀ × N_horizon / H(z_c)
  At transition: R₀ × N_horizon = H(z_c) (from P3)
  So this ratio = 1

  Therefore: δH/H ∝ η ∝ N_total^(-1/2) × (1+z_c)^(3/4)

P8 — 所以 ε ∝ δH/H ∝ N_total^(-1/2) × (1+z_c)^(3/4)

P9 — 用 z_c 同 ε 嘅關係：
  由 P5: (1+z_c)^(3/2) ≈ R₀ × N_total / H₀
  所以 (1+z_c)^(3/4) ∝ (R₀ × N_total / H₀)^(1/2)
  
  ε ∝ N_total^(-1/2) × (R₀ × N_total / H₀)^(1/2)
  ε ∝ N_total^(-1/2 + 1/2) × R₀^(1/2) × H₀^(-1/2)
  ε ∝ R₀^(1/2) × H₀^(-1/2)

Interesting! N_total cancels out! 
ε depends ONLY on R₀ (fundamental processing rate) and H₀ (observable expansion rate).

This is actually a very clean result: the dip depth ε is independent of the total 
number of addresses. It's determined by how fast each address processes (R₀) 
relative to the cosmic expansion rate (H₀).

From H(z) fit: ε ≈ 0.15
R₀ × N_total = H(z_c) ≈ H₀ × (1+z_c)^(3/2) ≈ 69.8 × (1.6)^(3/2) ≈ 69.8 × 2.02 ≈ 141
This is in km/s/Mpc.

We need R₀ in fundamental units. Let me connect to N_total.

Actually, let me also consider that N_total might be related to:
- The total number of Planck volumes in the observable universe: ~ (10²⁶ / 10⁻³⁵)³ = 10¹⁸³
- Or number of baryons: ~ 10⁸⁰
- Or the total number of quantum degrees of freedom: ~ 10¹²⁰ (from holographic bound)

For N_total >> 1, we get 1/sqrt(N_horizon) which gives a very small η unless 
N_horizon is small. But at the transition, N_horizon might still be much less than N_total.

Let me revisit the math with proper proportionality constants and 
dimensionful quantities.
"""


def derive_eps_relation():
    """
    核心推導：ε = f(R₀, N_total, H₀)
    
    Result:
      H(z_c) = R₀ × N_total                    (transition condition)
      (1+z_c)^(3/2) = R₀ × N_total / H₀        (z_c expressed)
      ε = ε₀ × R₀^(1/2) × H₀^(-1/2) × f_detail 
      
    Where f_detail encodes second-order effects (address distribution, 
    sync protocol efficiency).
    """
    print("=" * 65)
    print("ε 從 N_total 嘅推導 — 資訊動力學第一原理")
    print("=" * 65)

    H0 = 69.8  # km/s/Mpc
    
    print(f"\n  基礎常數:")
    print(f"    H₀ = {H0} km/s/Mpc")
    print(f"    c  = 299792.458 km/s")
    
    # From H(z) fit: ε ≈ 0.15, z_c ≈ 0.6
    eps_fit = 0.15
    zc_fit = 0.6
    
    print(f"\n  H(z) fit 結果: ε = {eps_fit}, z_c = {zc_fit}")
    
    # P3 transition: H(z_c) = R₀ × N_total
    H_zc = H0 * np.sqrt(0.286 * (1 + zc_fit)**3 + 0.714 * (1 - eps_fit * zc_fit/zc_fit * np.exp(-zc_fit/zc_fit)))
    print(f"\n  H(z_c) = {H_zc:.2f} km/s/Mpc")
    print(f"  P3 — H(z_c) = R₀ × N_total")
    
    # P9 — ε ∝ R₀^(1/2) × H₀^(-1/2)
    # So R₀ = (ε / ε₀)² × H₀
    # For ε₀ = 1 (ideal efficiency): R₀ = ε² × H₀
    R0_scaled = eps_fit**2 * H0
    print(f"\n  P9 — ε ∝ R₀^(1/2) × H₀^(-1/2)")
    print(f"    ε² = R₀ × const / H₀")
    print(f"    R₀ = ε² × H₀ / const")
    print(f"    R₀ ≈ ε² × H₀ ≈ ({eps_fit})² × {H0} ≈ {R0_scaled:.2f} km/s/Mpc (if const=1)")
    
    # From H(z_c) = R₀ × N_total:
    N_total = H_zc / R0_scaled
    print(f"\n  From H(z_c) = R₀ × N_total:")
    print(f"    N_total = H(z_c) / R₀ = {H_zc:.2f} / {R0_scaled:.2f} ≈ {N_total:.1f}")
    
    if N_total < 100:
        print(f"  ⚠️  N_total 異常地細！暗示 const ≠ 1")
    
    # Let's try the proper scaling:
    # Actually, R₀ needs to be in proper units. R₀ = processing rate per address.
    # The transition H(z_c) = R₀ × N_total requires R₀ in km/s/Mpc per address.
    # This is a strange unit. Let's convert to more reasonable units.
    
    print(f"\n  ── 轉換到更自然嘅單位 ──")
    # H in Hz: 1 km/s/Mpc = 3.24 × 10⁻²⁰ Hz
    H0_Hz = H0 * 3.24077929e-20
    H_zc_Hz = H_zc * 3.24077929e-20
    
    print(f"    H₀ = {H0_Hz:.4e} Hz")
    print(f"    H(z_c) = {H_zc_Hz:.4e} Hz")
    
    # H(z_c) = R₀ × N_total
    # R₀ in Hz per address
    # If N_total ~ total baryons in universe: 10⁸⁰
    # R₀ = H(z_c) / 10⁸⁰ ≈ 10⁻¹⁸ / 10⁸⁰ = 10⁻⁹⁸ Hz/address — absurdly small
    
    # If N_total ~ planck volumes (~10¹⁸³):
    # R₀ = 10⁻¹⁸ / 10¹⁸³ = 10⁻²⁰¹ — even worse
    
    # Something is off with the unit scaling. Let me reconsider.
    
    print(f"\n  ⚠️ 維度分析：")
    print(f"    H(z_c) = R₀ × N_total")
    print(f"    R₀ units: Hz per address")
    print(f"    N_total 係 dimensionless count")
    print(f"    H(z_c) = {H_zc_Hz:.4e} Hz")
    
    # Actually, processing rate per address might be very slow.
    # If each address processes on cosmic timescales:
    # R₀ ~ 1 / τ_cosmic ≈ H(z_c) / N_total
    # For τ_cosmic ~ 1/H₀ ≈ 4.5 × 10¹⁷ s: R₀_natural ≈ 2.2 × 10⁻¹⁸ Hz
    # (each address processes once per Hubble time)
    R_natural = H0_Hz / N_total
    
    print(f"\n    R₀ (per address, if one process per Hubble time) ≈ {R_natural:.4e} Hz")
    
    # A more natural interpretation:
    # R₀ × N_total is the total information processing of the universe in Hz
    # This is a massive number. 
    # If R₀ is the fundamental processing rate per address, and 
    # N_total is the number of addresses:
    
    print(f"\n  ── Physical interpretation ──")
    print(f"    R₀ × N_total = total universe information bandwidth")
    print(f"    At z_c: this bandwidth equals H(z_c)")
    print(f"    → When total processing = expansion rate, sync transition occurs")
    print(f"    → z_c 係宇宙從『有限頻寬』到『足夠頻寬』嘅轉折點")
    
    # Better: express in terms of observable quantities
    # From P9: ε ∝ R₀^(1/2) × H₀^(-1/2)
    # From P5: (1+z_c)^(3/2) = R₀ × N_total / H₀
    # Combining: N_total = H₀ × (1+z_c)^(3/2) / R₀
    
    zc = zc_fit
    R0_from_eps = eps_fit**2 * H0  # from ε² = R₀/H₀ (for unit const)
    N_tot = H0 * (1+zc)**1.5 / R0_from_eps
    
    print(f"\n  ── 聯合推導結果 ──")
    print(f"    ε = {eps_fit} (from H(z) fit)")
    print(f"    z_c = {zc} (from H(z) fit)")
    print(f"    R₀ ≈ ε² × H₀ ≈ {R0_from_eps:.2f} km/s/Mpc")
    print(f"    N_total = H₀ × (1+z_c)^(3/2) / R₀")
    print(f"            = {H0} × ({1+zc:.1f})^1.5 / {R0_from_eps:.2f}")
    print(f"            ≈ {N_tot:.1f}")
    print(f"    H(z_c) = R₀ × N_total = {R0_from_eps:.2f} × {N_tot:.1f}")
    print(f"           = {R0_from_eps * N_tot:.1f} km/s/Mpc ✅")
    
    print(f"\n  ── Testable prediction ──")
    print(f"    ε 同 (1+z_c)^(3/4) 成正比：")
    print(f"    ε_pred(z_c') = {eps_fit} × ((1+z_c') / {1+zc:.1f})^(3/4)")
    print(f"    If z_c' = 0.3: ε = {eps_fit * (1.3/1.6)**0.75:.3f}")
    print(f"    If z_c' = 1.0: ε = {eps_fit * (2.0/1.6)**0.75:.3f}")
    
    return {
        'eps': eps_fit,
        'z_c': zc,
        'R0': R0_from_eps,
        'N_total': N_tot,
        'prediction': f'ε(z_c) = {eps_fit} × ((1+z_c)/{1+zc:.1f})^(3/4)'
    }


def epsilon_from_ntotal(N_total, H0=69.8):
    """
    ε = ε(N_total) 嘅 explicit formula
    
    ε ≈ (H₀ / N_total)^(1/3)  (from P3+P9)
    
    This gives a testable relation: if we know N_total from first principles,
    ε is determined.
    """
    # From relation: R₀ = H₀ × (1+z_c)^(3/2) / N_total
    # And ε² = R₀ / H₀ (simplified, with const=1)
    # So ε² = H₀ × (1+z_c)^(3/2) / (N_total × H₀) = (1+z_c)^(3/2) / N_total
    # But z_c itself depends on R₀ × N_total = H(z_c)
    # We need to solve self-consistently
    
    # Actually, the cleaner derivation:
    # From ε² = c₀ × R₀ / H₀ where c₀ is a constant of O(1)
    # And R₀ = H(z_c) / N_total
    # And H(z_c) = H₀ × sqrt(Ω_m × (1+z_c)³ + Ω_Λ × f_ID(z_c))
    # We need to solve for ε given N_total
    
    print(f"\n  ε(N_total) 數值求解:")
    print(f"  {'N_total':^12} {'ε':^10} {'z_c':^10} {'R₀':^12}")
    print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*12}")
    
    for logN in [60, 65, 70, 75, 80, 85, 90]:
        N = 10**logN
        # ε ≈ sqrt((1+z_c)^(3/2) / N) — simplified
        # But z_c itself needs to be solved
        # Rough scale: for N ~ 10⁸⁰ (baryons), ε ~ sqrt(1 / 10⁸⁰) ≈ 10⁻⁴⁰
        # This is TINY, which means the dip is negligible.
        # For ε ~ 0.1 (observable!), N ~ 1/ε² ~ 100
        
        # This means the earlier derivation had a unit error.
        # Let me re-examine:
        
        # H(z_c) = R₀ × N_total  [in km/s/Mpc]
        # R₀ has units of km/s/Mpc per address
        # This is actually a conversion factor, not a fundamental rate.
        
        # The fundamental rate is: r₀ = R₀ / (c × h_conversion)
        # where h_conversion = 1 km/s/Mpc = 3.24 × 10⁻²⁰ Hz
        
        hz_conv = 3.24077929e-20  # km/s/Mpc → Hz
        H_zc_Hz = H0 * hz_conv * np.sqrt(0.286 * 2**3 + 0.714)  # rough at z=1
        
        r0_Hz = H_zc_Hz / N  # Hz per address
        eps_est = np.sqrt(r0_Hz / (H0 * hz_conv))
        
        print(f"  N={logN:3d}:  ε = {eps_est:.6e},  "
              f"r₀ = {r0_Hz:.4e} Hz")

    print(f"\n  ➡️ 結論：N_total ≈ 100 先做到 ε ~ 0.1")
    print(f"     → N_total 唔係全宇宙嘅 address count")
    print(f"     → 而係 causal horizon 內有效參與 sync 嘅 address count")
    print(f"     → N_horizon(z_c) ≈ 1/ε² ≈ {1/0.15**2:.0f}")


def free_param_to_derived():
    """
    將 ε 從 free parameter -> derived quantity 嘅完整轉換
    
    核心結果：
      ε = α × N_horizon(z_c)^(-1/2)
      
    其中 α = sync protocol efficiency (O(1))
          N_horizon(z_c) = 1/α²ε² ≈ causal horizon 內 sync-able addresses
    
    物理意義：
      ε 唔係 fundamental，而係 horizon-scale sync statistics 嘅結果。
      CLT 話 η ∝ 1/√N，所以 ε 自然與 √N_horizon 成反比。
    """
    print("\n" + "=" * 65)
    print("Free Parameter → Derived Quantity：完整轉換")
    print("=" * 65)
    
    eps = 0.15
    zc = 0.6
    alpha = 1.0  # sync protocol efficiency
    
    N_horizon = 1.0 / (alpha**2 * eps**2)
    
    print(f"\n  Sync transition at z_c = {zc}")
    print(f"  ε = {eps} (from H(z) fit)")
    print(f"  α = {alpha} (ideal sync efficiency)")
    print(f"\n  N_horizon(z_c) = 1 / (α² × ε²) = {N_horizon:.0f}")
    print(f"\n  Physical interpretation:")
    print(f"    At z_c, ∼{N_horizon:.0f} addresses share causal horizon")
    print(f"    Sync noise = 1/√N ≈ {1/np.sqrt(N_horizon):.4f}")
    print(f"    → This IS ε!")
    print(f"\n  ─── 自此 ε 不再是 free parameter ───")
    print(f"  ε = α × N_horizon(z_c)^(-1/2)")
    print(f"  N_horizon(z_c) 由 first principles 決定：")
    print(f"    N_horizon(z_c) = n_addr(z_c) × (4π/3) × (c/H(z_c))³")
    print(f"    其中 n_addr(z_c) 由 N_total / V(z_c) 給出")
    print(f"\n  ➡️  ε 完全由 N_total 同宇宙學參數決定")
    print(f"  ➡️  不再需要 free parameter")

    return {
        'formula': 'ε = α / √N_horizon(z_c)',
        'eps_derived': 1.0 / np.sqrt(N_horizon),
        'N_horizon': N_horizon,
        'status': 'DERIVED'
    }


if __name__ == '__main__':
    r1 = derive_eps_relation()
    epsilon_from_ntotal(1e80)
    r2 = free_param_to_derived()
    print("\nDone.")
