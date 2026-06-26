"""
f_secular_v3.py — Corrected: 1D particle horizon N_horizon + √N scaling
=======================================================================
Key revisions:
  1. N_horizon(z) = 44 × D_C(z)/D_C(z_c)  [1D: linear, not cubic]
  2. At z→0: N→0, at z→∞: N→~200 (finite)
  3. dρ_sync/dz ∝ √N_horizon(z)/(1+z)  [total sync power]
"""

import numpy as np
from scipy import integrate

C = 299792.458
H0_best = 69.0; Om = 0.300; Ol = 1-Om; eps_f = 0.25; zc_f = 0.70
Nh_zc = 44; z_c = 0.70

class Model:
    def __init__(self, H0=69.0, Om=0.300, eps=0.25, zc=0.70, lam=1.0):
        self.H0=H0; self.Om=Om; self.Ol=1-Om; self.eps=eps; self.zc=zc; self.lam=lam
        h=H0/100; self.Om_r=4.15e-5/h**2
        # Precompute D_C(z_c) for normalization
        self.DC_zc = None
    
    def f_trans(self, z):
        x=np.asarray(z,float)/self.zc
        return 1-self.eps*x*np.exp(-x)
    
    def E(self, z):
        z=np.asarray(z,float)
        return np.sqrt(self.Om*(1+z)**3+self.Om_r*(1+z)**4+self.Ol*self.f_trans(z))
    
    def D_C(self, z):
        """Comoving particle horizon at epoch z (Mpc)"""
        f=lambda x:1/self.E(x)
        I,_=integrate.quad(f,0,z,limit=2000,epsabs=1e-12)
        return C/self.H0*I
    
    def N1D(self, z, dc_zc=None):
        """1D causal domain count at redshift z"""
        dc = self.D_C(z)
        if dc_zc is None:
            dc_zc = self.D_C(self.zc)
        return 44 * dc / dc_zc if dc_zc > 0 else 0.0
    
    def f_secular(self, z_target):
        """
        df/dz = -λ × √N_1D(z) / (1+z)
        Sync consumption ∝ total sync rate = pairs × 1/√N = N × 1/√N = √N
        """
        dc_zc = self.D_C(self.zc)
        integrand = lambda x: np.sqrt(max(44*self.D_C(x)/dc_zc,1e-30))/(1+x)
        # Integrate from 0 to z_target (past)
        I,_=integrate.quad(integrand,0,z_target,limit=1000,epsabs=1e-12)
        return 1.0 - self.lam * I
    
    def f_total(self, z):
        return self.f_secular(z) * self.f_trans(z)


m = Model(lam=1.0)
dc_zc = m.D_C(m.zc)

print("=" * 72)
print("  f_SECULAR(z) v3 — 1D particle horizon + √N scaling")
print("=" * 72)

# Step 1: N_horizon evolution
print(f"\n─── STEP 1: N_horizon(z) = 44 × D_C(z)/D_C(z_c={m.zc}) ───")
print(f"  D_C(z_c) = {dc_zc:.0f} Mpc")
print(f"  D_C(∞) ≈ {m.D_C(1000):.0f} Mpc")
print(f"  N1D(∞) ≈ {44*m.D_C(1000)/dc_zc:.0f}")
print()

zs = [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 10.0, 100.0]
print(f"  {'z':>6}  {'D_C(Mpc)':>10}  {'N1D':>6}  {'√N1D':>7}  {'df/dz':>10}")
print(f"  {'─'*6}  {'─'*10}  {'─'*6}  {'─'*7}  {'─'*10}")
for z in zs:
    dc = m.D_C(z)
    n = 44*dc/dc_zc
    sqrt_n = np.sqrt(max(n,1e-30))
    integrand = sqrt_n/(1+z)
    print(f"  {z:6.2f}  {dc:10.0f}  {n:6.1f}  {sqrt_n:7.3f}  {integrand:10.4f}")

# Step 2: f_secular computation
print(f"\n─── STEP 2: f_secular(z) = 1 - λ∫ √N/(1+z) dz ───")
for lam in [0.001, 0.01, 0.05, 0.1, 0.5, 1.0]:
    ml = Model(lam=lam)
    print(f"\n  λ = {lam}:")
    for z in [0.01, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 10.0, 100.0]:
        fs = ml.f_secular(z)
        print(f"    z={z:6.2f}  f_sec={fs:9.6f}  dev={fs-1:+.3e}")

# Step 3: Secular ≪ dip at z<2?
print(f"\n─── STEP 3: Secular vs dip ───")
for lam in [0.001, 0.01, 0.05, 0.1]:
    ml = Model(lam=lam)
    all_ok = True
    for z in [0.3, 0.5, 0.7, 1.0, 2.0]:
        sec = abs(ml.f_secular(z)-1)
        dip = abs(ml.f_trans(z)-1)
        ratio = sec/dip*100 if dip>0 else 0
        ok = ratio < 10
        if not ok: all_ok = False
        print(f"  λ={lam:.3f} z={z:.1f}: sec={sec:.2e}, dip={dip:.4f}, ratio={ratio:.2f}% {'✓' if ok else '✗'}")
    best_lam = lam if all_ok else None
    if all_ok:
        print(f"  → λ = {lam} satisfies secular ≪ dip ✓")
        best_ok_lam = lam

# Step 4: Find max λ such that secular < 10% of dip at z<2
print(f"\n─── STEP 4: Constraining λ ───")
for lam_test in [0.001, 0.003, 0.005, 0.007, 0.01, 0.02, 0.05]:
    ml = Model(lam=lam_test)
    max_ratio = 0
    for z in [0.3, 0.5, 0.7, 1.0, 2.0]:
        sec = abs(ml.f_secular(z)-1)
        dip = abs(ml.f_trans(z)-1)
        ratio = sec/dip*100 if dip>0 else 0
        max_ratio = max(max_ratio, ratio)
    ok = max_ratio < 10
    print(f"  λ={lam_test:.4f}: max_sec/dip_ratio={max_ratio:.2f}%  {'✓' if ok else '✗'}")

# Step 5: Future asymptotic
print(f"\n─── STEP 5: Asymptotic future ───")
I_total,_ = integrate.quad(lambda x: np.sqrt(max(44*m.D_C(x)/dc_zc,1e-30))/(1+x), 0, 1e6, limit=5000)
print(f"  I_total(∫₀^∞ √N/(1+z) dz) = {I_total:.4f}")
for lam_test in [0.001, 0.005, 0.01]:
    print(f"  λ={lam_test:.4f}: f_secular(∞) = {1-lam_test*I_total:.6f}")
print(f"\n  f_secular(∞) converges to 1 - λ·I∞")
print(f"  I∞ ≈ {I_total:.4f}")
print(f"  → NEED λ ≈ 1/I∞ ≈ {1/I_total:.4f} for f_secular→0")
print(f"  → But λ ≈ {0.001:.4f} needed for secular ≪ dip at z<2")
print(f"  → At best λ, f_secular(z→∞) ≈ {1-0.001*I_total:.4f}")

# Step 6: Gemini vs sync resolution
print(f"\n─── STEP 6: Gemini's islandization vs IDM's global sync ───")
print(f"""
  For λ ≈ 0.001 (secular negligible):
    f_secular(∞) ≈ {1-0.001*I_total:.4f}
    → DE stays at nearly full Λ₀ forever
    → De Sitter expected (Gemini's islandization partially correct)
    → But residual sync continues slowly
  
  For λ ≈ 1 (sync dominant):
    f_secular(∞) ≈ {1-1.0*I_total:.4f}
    → Unphysical (f < 0)
    → λ must be constrained by data
  
  λ ≈ 1/I∞ ≈ {1/I_total:.5f}:
    f_secular(∞) = 0 exactly
    → sync consumes ALL DE asymptotically
    → deceleration, no De Sitter, no islandization
    → global sync (IDM's ultimate target)

  The transition between these regimes determines the ultimate fate.
  Current data constrains λ (need joint fit).
""")