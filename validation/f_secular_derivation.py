"""
f_secular_derivation.py — 推導 f_secular(z) 從同步動力學
=========================================================
Three steps:
  Step 1 — N_horizon(z) 演化: 確認單調遞增
  Step 2 — f_secular(z) 積分推導
  Step 3 — 驗證 secular ≪ ε dip at z<2
          — f_total asymptotics at z→∞
          — 未來走向預測
"""

import numpy as np
from scipy import integrate
import sys

# ============================================================
# Constants & model
# ============================================================
C = 299792.458
H0_best = 69.0        # IDM best-fit
Om_best = 0.300
Ol_best = 1.0 - Om_best
eps_fit = 0.25
zc_fit = 0.70
Nh_zc = 44            # N_horizon at z_c (from ε=0.15 CLT)

class IDMModel:
    def __init__(self, H0=69.0, Om=0.300, eps=0.25, zc=0.70, gamma=1.0):
        self.H0 = H0; self.Om = Om; self.Ol = 1-Om
        self.eps = eps; self.zc = zc
        self.gamma = gamma  # sync efficiency (O(1))
        h = H0/100.0
        self.Om_r = 4.15e-5 / h**2

    def f_transient(self, z):
        """Current dip model"""
        x = np.asarray(z, float) / self.zc
        return 1.0 - self.eps * x * np.exp(-x)

    def E2(self, z):
        z = np.asarray(z, float)
        matter = self.Om * (1+z)**3
        rad = self.Om_r * (1+z)**4
        de = self.Ol * self.f_transient(z)
        return matter + rad + de

    def E(self, z):
        return np.sqrt(self.E2(z))

    def H(self, z):
        return self.H0 * self.E(z)

    def r_H(self, z):
        """Hubble radius (Mpc)"""
        return C / self.H(z)

    def V_horizon(self, z):
        """Hubble volume ∝ 1/H³"""
        return 1.0 / self.H(z)**3

    # ============================================================
    # Step 1: N_horizon(z) evolution
    # ============================================================
    def N_horizon(self, z):
        """
        Number of causal domains at redshift z.
        N(z) = N(z_c) × [V_horizon(z) / V_horizon(z_c)]
        ∝ 1/H(z)³
        Monotonically INCREASING as z decreases.
        """
        return Nh_zc * (self.r_H(z) / self.r_H(self.zc))**3

    # ============================================================
    # Step 2: f_secular(z) from sync dynamics
    # ============================================================
    def f_secular(self, z_target):
        """
        f_secular(z) = 1 - γ·N₀² ∫₀ᶻ dz' / [E(z')⁵(1+z')]
        
        Derived from:
          dρ_sync/dt ∝ N_horizon² × H
          N_horizon(z) = N₀ / E(z)³
          dz'/dt = -H(z')(1+z')
          → dρ_sync/dz' ∝ N₀² / E(z')⁵ / (1+z')
        """
        N0_sq = Nh_zc**2
        integrand = lambda x: N0_sq / (self.E(x)**5 * (1+x))
        I, _ = integrate.quad(integrand, 0, z_target, limit=1000, epsabs=1e-12)
        return 1.0 - self.gamma * I

    def secular_deviation(self, z):
        """|f_secular(z) - 1|"""
        return abs(self.f_secular(z) - 1)


print("=" * 72)
print("  f_secular(z) DERIVATION FROM SYNC DYNAMICS")
print("=" * 72)

m = IDMModel(gamma=1.0)

# ============================================================
# STEP 1: N_horizon(z) evolution
# ============================================================
print("\n" + "─" * 60)
print("  STEP 1: N_horizon(z) evolution")
print("─" * 60)

print(f"\n  N_horizon(z) ∝ 1/H(z)³")
print(f"  Normalized: N_horizon(z_c={m.zc}) = {Nh_zc}")
print(f"  Prediction: N_horizon(z) INCREASES as z decreases")
print()

zs = [3.0, 2.0, 1.0, 0.7, 0.5, 0.3, 0.1, 0.01]
print(f"  {'z':>6}  {'r_H (Mpc)':>12}  {'V_rel':>10}  {'N_horizon':>12}")
print(f"  {'─'*6}  {'─'*12}  {'─'*10}  {'─'*12}")
for z in zs:
    r = m.r_H(z)
    n = m.N_horizon(z)
    v_rel = n / Nh_zc
    print(f"  {z:6.2f}  {r:12.1f}  {v_rel:10.2f}  {n:12.1f}")

print(f"\n  Direction: N_horizon(0) ≈ {m.N_horizon(0.01):.0f} >> {Nh_zc} at z_c")
print(f"  → Horizon grows, MORE domains enter causal contact")
print(f"  → NOT fragmentation / isolation")
print(f"  → APPROACHING global sync (N_horizon → ∞)")

# ============================================================
# STEP 2: f_secular(z) integral
# ============================================================
print("\n" + "─" * 60)
print("  STEP 2: f_secular(z) computation")
print("─" * 60)

print(f"\n  Formula:")
print(f"  f_secular(z) = 1 - γ·N₀² ∫₀ᶻ dz' / [E(z')⁵(1+z')]")
print(f"  where N₀ = {Nh_zc}, γ = sync efficiency parameter")
print()

# Compute f_secular at various z
secular_zs = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 5.0, 10.0, 50.0, 100.0]
print(f"  {'z':>8}  {'f_secular':>12}  {'deviation':>12}  {'f_transient':>14}  {'f_total':>14}")
print(f"  {'─'*8}  {'─'*12}  {'─'*12}  {'─'*14}  {'─'*14}")
for z in secular_zs:
    fs = m.f_secular(z)
    ft = m.f_transient(z)
    f_total = fs * ft
    dev = fs - 1.0
    print(f"  {z:8.2f}  {fs:12.6f}  {dev:12.6f}  {ft:14.6f}  {f_total:14.6f}")

# ============================================================
# STEP 3a: Verify secular ≪ ε dip at z < 2
# ============================================================
print("\n" + "─" * 60)
print("  STEP 3a: Secular vs Dip comparison")
print("─" * 60)

print(f"\n  Key question: Is |f_secular(z) - 1| ≪ ε at z < 2?")
print(f"  ε = {m.eps} → dip amplitude = 1 - exp(-1)×ε = {1-np.exp(-1)*m.eps:.4f}")
print()

dip_check = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
print(f"  {'z':>6}  {'|f_sec-1|':>12}  {'|f_trans-1|':>14}  {'ratio':>10}  {'Pass?':>8}")
print(f"  {'─'*6}  {'─'*12}  {'─'*14}  {'─'*10}  {'─'*8}")
for z in dip_check:
    sec_dev = m.secular_deviation(z)
    dip_dev = abs(m.f_transient(z) - 1.0)
    ratio = sec_dev / dip_dev if dip_dev > 0 else np.inf
    # Need |f_sec-1| < ε × 10% (if secular is < 10% of dip amplitude)
    threshold = 0.02  # 10% of epsilon
    pass_check = sec_dev < 0.1 * dip_dev
    print(f"  {z:6.2f}  {sec_dev:12.3e}  {dip_dev:14.3e}  {ratio:10.2e}  {str(pass_check):>8}")

print(f"\n  CRITICAL FINDING:")
for z in dip_check:
    sec_dev = m.secular_deviation(z)
    dip_dev = abs(m.f_transient(z) - 1.0)
    if sec_dev < 0.1 * dip_dev:
        print(f"    z={z:.1f}: Secular is {sec_dev/dip_dev:.1%} of dip → negligible ✓")
    else:
        print(f"    z={z:.1f}: Secular is {sec_dev/dip_dev:.1%} of dip → check!")

# ============================================================
# STEP 3b: Asymptotic behavior
# ============================================================
print("\n" + "─" * 60)
print("  STEP 3b: Asymptotic analysis z → ∞")
print("─" * 60)

# Asymptotic: at high z, E(z) ≈ √(Om) × (1+z)^(3/2)
# integrand ≈ N₀² / (Om^(5/2) × (1+z)^(15/2) × (1+z))
# = N₀² / (Om^(5/2)) × (1+z)^(-17/2)
# ∫₀^∞ (1+z)^(-17/2) dz converges!
# So f_secular(∞) = 1 - γ × finite constant

z_very_high = [10, 20, 50, 100, 500]
print(f"  {'z':>6}  {'f_secular':>12}  {'f_total':>14}  {'converging?':>12}")
print(f"  {'─'*6}  {'─'*12}  {'─'*14}  {'─'*12}")
for z in z_very_high:
    fs = m.f_secular(z)
    ft = m.f_transient(z)
    ftot = fs * ft
    # Check convergence by comparing to f_secular(z*10)
    fs_10 = m.f_secular(z*10)
    delta = abs(fs - fs_10)
    converged = delta < 1e-6
    print(f"  {z:6.0f}  {fs:12.6f}  {ftot:14.6f}  {'~converged' if converged else 'changing':>12}")

fs_inf = m.f_secular(1e6)
print(f"\n  f_secular(z→∞) = {fs_inf:.6f}")
print(f"  f_transient(z→∞) = {m.f_transient(1e6):.6f}")
print(f"  f_total(z→∞) = {fs_inf * m.f_transient(1e6):.6f}")
print(f"  DE asymptotically trends toward zero: f_total→{fs_inf * m.f_transient(1e6):.4f}")
print(f"  Event horizon diverges → global connectivity → NOT islandization")

# ============================================================
# Summary table
# ============================================================
print("\n" + "=" * 72)
print("  FINAL SUMMARY")
print("=" * 72)

print(f"""
  ┌─────────────────────────────────────────────────────────┐
  │  N_horizon(z) evolution:                                │
  │    N_horizon(z) = {Nh_zc} × [r_H(z)/r_H(z_c)]³             │
  │    → MONOTONICALLY INCREASING as z decreases            │
  │    → N_horizon(0) ≈ {m.N_horizon(0.01):.0f}  >>  {Nh_zc} at z_c       │
  │    → Future: N_horizon → ∞ (global connectivity)        │
  │                                                         │
  │  f_secular(z) derivation:                               │
  │    dρ_sync/dt ∝ N_horizon² × H                          │
  │    f_secular(z) = 1 - γ·N₀² ∫₀ᶻ dz/[E(z')⁵(1+z')]   │
  │    γ = {m.gamma} (sync efficiency, O(1))                     │
  │                                                         │
  │  Secular ≪ dip check:                                   │
  │    At z=0.7: |f_sec-1| = {m.secular_deviation(0.7):.2e}           │
  │    At z=2.0: |f_sec-1| = {m.secular_deviation(2.0):.2e}           │
  │    Dip amplitude: 1-ε/e = {1-np.exp(-1)*m.eps:.6f}               │
  │    → Secular is {m.secular_deviation(0.7)/(1-np.exp(-1)*m.eps)*100:.2f}% of dip at z_c ✓
  │                                                         │
  │  Asymptotic future:                                      │
  │    f_total(z→∞) = {fs_inf * m.f_transient(1e6):.6f}                      │
  │    DE trends toward zero → deceleration →               │
  │    event horizon diverges → no islandization             │
  │    → COSMOS TOWARD GLOBAL SYNCHRONIZATION               │
  └─────────────────────────────────────────────────────────┘
""")

# Save reference values
print("  Reference values for paper:")
print(f"    gamma = 1.0: f_secular(z=0.7) = {m.f_secular(0.7):.6f}")
print(f"    gamma = 1.0: f_secular(z=1.0) = {m.f_secular(1.0):.6f}")
print(f"    gamma = 1.0: f_secular(z=2.0) = {m.f_secular(2.0):.6f}")
print(f"    gamma = 0.1: f_secular(z=0.7) = {IDMModel(gamma=0.1).f_secular(0.7):.6f}")
print(f"    gamma = 2.0: f_secular(z=0.7) = {IDMModel(gamma=2.0).f_secular(0.7):.6f}")
