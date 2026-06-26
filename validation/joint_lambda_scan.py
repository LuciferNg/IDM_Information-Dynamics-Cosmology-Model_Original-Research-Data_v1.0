"""
joint_lambda_scan.py — Joint fit including λ (secular sync efficiency)
======================================================================
Strategy:
  1. Pre-compute I(z) = ∫₀ᶻ √N(z')/(1+z') dz'  [independent of λ]
  2. Scan λ, compute f_total = (1-λ·I) × f_transient
  3. Find best-fit λ and its constraint on DE consumption
"""

import numpy as np
from scipy import integrate
from itertools import product
import time

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

Nh_zc = 44

# ============================================================
# Model
# ============================================================
class Model:
    def __init__(self, H0=69.0, Om=0.300, eps=0.25, zc=0.70, lam=0.001):
        self.H0=H0; self.Om=Om; self.Ol=1-Om; self.eps=eps; self.zc=zc; self.lam=lam
        h=H0/100; self.Om_r=4.15e-5/h**2
        self._I_cache = {}  # cache for secular integral

    def f_trans(self, z):
        x=np.asarray(z,float)/self.zc
        return 1-self.eps*x*np.exp(-x)

    def E(self, z):
        z=np.asarray(z,float)
        return np.sqrt(self.Om*(1+z)**3+self.Om_r*(1+z)**4+self.Ol*self.f_trans(z))

    def H(self, z):
        return self.H0*self.E(z)

    def D_C(self, z):
        f=lambda x:1/self.E(x)
        I,_=integrate.quad(f,0,z,limit=2000,epsabs=1e-12,epsrel=1e-12)
        return C/self.H0*I

    def DC_zc(self):
        return self.D_C(self.zc)

    def N1D(self, z, dc_zc):
        """1D domain count at z"""
        dc = self.D_C(z)
        return 44*dc/dc_zc

    def I_secular(self, z, dc_zc):
        """I(z) = ∫₀ᶻ √N(z')/(1+z') dz' — independent of λ"""
        if z in self._I_cache:
            return self._I_cache[z]
        def integrand(x):
            dc = self.D_C(x)
            n = 44*dc/dc_zc
            return np.sqrt(max(n,1e-30))/(1+x)
        I,_ = integrate.quad(integrand,0,z,limit=2000,epsabs=1e-12)
        self._I_cache[z] = I
        return I

    def f_secular(self, z, dc_zc):
        return 1 - self.lam * self.I_secular(z, dc_zc)

    def DM(self, z):
        return self.D_C(z)

    def DH(self, z):
        return C/self.H(z)

    def DV(self, z):
        DM=self.DM(z); H_z=self.H(z)
        return (C*z*DM**2/H_z)**(1/3)

    # CMB
    def I_zrec(self):
        f=lambda x:1/self.E(x)
        I,_=integrate.quad(f,0,Z_REC,limit=2000,epsabs=1e-12,epsrel=1e-12)
        return I

    def shift_R(self):
        return np.sqrt(self.Om)*self.I_zrec()

    def l_A(self):
        return np.pi*self.D_C(Z_REC)/R_S_CMB


# ============================================================
# Data
# ============================================================

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

BAO_MEAN = np.array([
    7.92512927,13.62003080,20.98334647,
    16.84645313,20.07872919,
    21.70841761,17.87612922,
    27.78720817,13.82372285,
    26.07217182,
    39.70838281,8.52256583,
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
    zs,mbs,dmbs=[],[],[]
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip(): continue
            p=line.strip().split()
            if len(p)<6: continue
            try:
                z=float(p[1])
                if z>=0.01:
                    zs.append(z);mbs.append(float(p[4]));dmbs.append(float(p[5]))
            except: continue
    return np.array(zs),np.array(mbs),np.array(dmbs)


# ============================================================
# χ² functions
# ============================================================

def chi2_hz(m):
    r=m.H(HZ[:,0])-HZ[:,1]
    # f_total at H(z) redshifts
    return np.sum((r/HZ[:,2])**2)

def chi2_sne(m,zs,mbs,dmbs):
    dc_zc=m.DC_zc()
    mu=np.array([5*np.log10(
        (1+z)*m.DM(z)*m.f_secular(z,dc_zc)  # D_L w/ secular corr
    )+25 for z in zs])
    w=1/dmbs**2
    resid=mu-mbs
    Mb=-np.sum(w*resid)/np.sum(w)
    return np.sum(w*(resid+Mb)**2),Mb

def chi2_bao(m):
    pred=np.zeros(12)
    for i in range(12):
        z=BAO_Z[i]; t=BAO_T[i]
        if t=='DM': pred[i]=m.DM(z)/R_S_BAO
        elif t=='DH': pred[i]=m.DH(z)/R_S_BAO
        elif t=='DV': pred[i]=m.DV(z)/R_S_BAO
    d=pred-BAO_MEAN
    return d@BAO_ICOV@d

def chi2_cmb(m):
    r=m.shift_R(); la=m.l_A()
    vr=PRE**2; vla=PLAE**2; crl=PRLAC*PRE*PLAE
    det=vr*vla-crl**2
    i_r=vla/det; i_la=vr/det; i_rl=-crl/det
    dr=r-PR; dla=la-PLA
    return i_r*dr**2+2*i_rl*dr*dla+i_la*dla**2,r,la

def eval_m(m,zs,mbs,dmbs):
    h=chi2_hz(m)
    s,Mb=chi2_sne(m,zs,mbs,dmbs)
    b=chi2_bao(m)
    c,r,la=chi2_cmb(m)
    return h+s+b+c,{'h':h,'s':s,'b':b,'c':c,'M':Mb,'R':r,'lA':la}


# ============================================================
# Main
# ============================================================
print("="*72)
print("  JOINT FIT WITH λ (SECULAR SYNC EFFICIENCY)")
print("="*72)

zs,mbs,dmbs=load_sne()
print(f"  SNe: {len(zs)}, H(z): {len(HZ)}, BAO: 12-vector, CMB: R+l_A")

t0=time.time()

print("\n─── Step 1: Baseline at best IDM (no λ) ───")
m0=Model(H0=69.0,Om=0.300,eps=0.25,zc=0.70,lam=0.0)
t0v,d0v=eval_m(m0,zs,mbs,dmbs)
print(f"  λ=0: χ²={t0v:.1f}")

# Step 2: Precompute reference integral values for the λ scan
print("\n─── Step 2: Precomputing I(z) for λ scan ───")
m_ref=Model(H0=69.0,Om=0.300,eps=0.25,zc=0.70,lam=0.001)
dc_zc=m_ref.DC_zc()

# We need I(z) at each BAO/HZ/SNe redshift + z_rec
# For efficiency, compute I(z) at key redshifts
print("  Computing I(z) at BAO redshifts...")
IBAO=[m_ref.I_secular(z,dc_zc) for z in BAO_Z]

print("  Computing I(z) at H(z) redshifts...")
IHZ=[m_ref.I_secular(z,dc_zc) for z in HZ[:,0]]

print("  Computing I(z_rec) for CMB...")
I_zrec=m_ref.I_secular(Z_REC,dc_zc)

# For SNe, we use approximation: D_L(z) = (1+z)*D_C(z)*f_secular(z)
# This is ≈ (1+z)*D_C(z)*(1-λ·I(z))
# For the SN χ², we need to marginalize over M, so we need the μ values
# We'll compute I(z) at the SNe redshifts for a subset
# But that's 1048 points... let's just sample

print("  Computing I(z) at H(z)+SNe+BAO reference grid...")
all_zs = sorted(set([z for z in HZ[:,0]] + BAO_Z + [0.01, 0.1, 0.3, 0.5, 0.7, 1.0, 2.0, 3.0, 10.0, Z_REC]))
I_grid = {z: m_ref.I_secular(z, dc_zc) for z in all_zs}

elapsed_p1 = time.time() - t0
print(f"  Precomputation done: {elapsed_p1:.0f}s")

# Step 3: Scan λ
print("\n─── Step 3: λ scan ───")
t2=time.time()

results=[]
for loglam in np.linspace(-5, 0, 26):  # λ = 10^(-5) to 10^0
    lam = 10**loglam
    # Check f_secular stays physical (>0) at all redshifts
    fsec_inf = 1 - lam * I_grid[100.0]
    if fsec_inf < 0:
        print(f"  λ={lam:.6f}: SKIP (f_sec<0 at high z)")
        continue
        
    m=Model(H0=69.0,Om=0.300,eps=0.25,zc=0.70,lam=lam, _I_cache={})
    t,d=eval_m(m,zs,mbs,dmbs)
    results.append((lam,t,d))
    
elapsed_l = time.time() - t2

results.sort(key=lambda x:x[1])

print(f"  Scan done: {elapsed_l:.0f}s")
print(f"  Evaluated {len(results)} λ values")

# Step 4: Report
print("\n─── RESULTS ───")
print(f"\n  λ (best) = {results[0][0]:.6f}")
print(f"  χ²_total = {results[0][1]:.1f}")
print(f"  vs λ=0: χ²={t0v:.1f}, Δχ²={t0v-results[0][1]:+.2f}")

print(f"\n  {'λ':>12}  {'χ²':>8}  {'χ²_H':>7}  {'χ²_S':>7}  {'χ²_B':>7}  {'χ²_C':>7}  {'f_sec(∞)':>10}")
print(f"  {'─'*12}  {'─'*8}  {'─'*7}  {'─'*7}  {'─'*7}  {'─'*7}  {'─'*10}")

# Show top 5 + λ=0 + λ where f_sec(∞)=0.5, 0.8
show_lams = set([r[0] for r in results[:5]] + [0.0])
for lam_val in [0.0, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3]:
    match=[r for r in results if abs(r[0]-lam_val)<1e-10]
    if match:
        r=match[0]
        fsec_inf=1-r[0]*I_grid[100.0]
        print(f"  {r[0]:12.6f}  {r[1]:8.1f}  {r[2]['h']:7.1f}  {r[2]['s']:7.1f}  {r[2]['b']:7.1f}  {r[2]['c']:7.1f}  {fsec_inf:10.4f}")

# Also show top 5 from sorted results
print("\n  Top 5 λ:")
for i,r in enumerate(results[:5]):
    fsec_inf=1-r[0]*I_grid[100.0]
    print(f"  {i+1}. λ={r[0]:.6f}  χ²={r[1]:.1f}  f_sec(∞)={fsec_inf:.4f}")

# Best λ details
best_lam = results[0][0]
best_d = results[0][2]
fsec_inf_best = 1 - best_lam * I_grid[100.0]

print(f"\n  Best-fit summary:")
print(f"    λ = {best_lam:.6f}")
print(f"    f_secular(∞) = {fsec_inf_best:.4f}")
print(f"    → {100*fsec_inf_best:.1f}% of DE remains fundamental")
print(f"    → {(1-fsec_inf_best)*100:.1f}% consumed by sync")
print()

# Per-dataset comparison
for k,lbl in [('h','H(z)'),('s','SNe'),('b','BAO'),('c','CMB')]:
    dv = results[0][2][k] - d0v[k]
    better = "λ>0" if dv<0 else "λ=0"
    print(f"    {lbl:5s}: Δχ²={dv:+7.2f}  →  {better}")

elapsed_total = time.time() - t0
print(f"\n  Total elapsed: {elapsed_total:.0f}s")
