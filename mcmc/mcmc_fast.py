"""mcmc_fast.py — Fast interpolation-based MCMC for IDM"""
import numpy as np, emcee, corner, time, os, sys
from scipy import integrate, interpolate as ip

C=299792.458; ZR=1089.92; RB=147.09; RC=144.43; PR=1.7502; PE=0.0046
PL=301.808; PLE=0.090; PRC=0.47

# ─── Data ───
HZ=np.array([[0.07,69,19.5],[0.09,69,12],[0.12,68.6,26.2],[0.17,83,8],[0.179,75,4],
  [0.199,75,5],[0.20,72.9,29.6],[0.27,77,14],[0.28,88.8,36.6],[0.352,83,14],
  [0.38,81.5,1.9],[0.4,82.8,1.2],[0.4004,77,10.2],[0.4247,87.1,11.2],
  [0.44,82.6,7.8],[0.4448,92.8,12.9],[0.47,89,49.6],[0.4783,80.9,9],[0.48,97,62],
  [0.51,90.4,1.9],[0.56,93.34,2.3],[0.57,92.4,4.5],[0.593,104,13],[0.6,87.9,6.1],
  [0.61,97.3,2.1],[0.68,92,8],[0.73,97.3,7],[0.78,97.6,4],[0.781,105,12],
  [0.875,125,17],[0.88,90,40],[0.9,117,23],[1.037,154,20],[1.3,168,17],
  [1.363,160,33.6],[1.43,177,18],[1.53,140,14],[1.75,202,40],[1.965,186.5,50.4]])
BM=np.array([7.92512927,13.62003080,20.98334647,16.84645313,20.07872919,21.70841761,
  17.87612922,27.78720817,13.82372285,26.07217182,39.70838281,8.52256583])
BC=np.array([[2.2723e-02,0,0,0,0,0,0,0,0,0,0,0],[0,6.3466e-02,-6.8534e-02,0,0,0,0,0,0,0,0,0],
  [0,-6.8534e-02,3.7297e-01,0,0,0,0,0,0,0,0,0],[0,0,0,1.0198e-01,-7.9940e-02,0,0,0,0,0,0,0],
  [0,0,0,-7.9940e-02,3.5445e-01,0,0,0,0,0,0,0],[0,0,0,0,0,7.9568e-02,-3.8011e-02,0,0,0,0,0],
  [0,0,0,0,0,-3.8011e-02,1.1994e-01,0,0,0,0,0],[0,0,0,0,0,0,0,4.7657e-01,-1.2941e-01,0,0,0],
  [0,0,0,0,0,0,0,-1.2941e-01,1.7827e-01,0,0,0],[0,0,0,0,0,0,0,0,0,4.4713e-01,0,0],
  [0,0,0,0,0,0,0,0,0,0,8.8975e-01,-7.6948e-02],
  [0,0,0,0,0,0,0,0,0,0,-7.6948e-02,2.9186e-02]])
BZ=np.array([0.295,0.510,0.510,0.706,0.706,0.930,0.930,1.317,1.317,1.491,2.330,2.330])
BT=['DV','DM','DH','DM','DH','DM','DH','DM','DH','DV','DM','DH']; BIC=np.linalg.inv(BC)

zs,ms,dms=[],[],[]
with open('/tmp/pantheon_raw.txt') as f:
    for l in f:
        if l.startswith("#")or not l.strip():continue
        p=l.strip().split()
        if len(p)<6:continue
        z=float(p[1])
        if z>=0.01:zs.append(z);ms.append(float(p[4]));dms.append(float(p[5]))
zs=np.array(zs);ms=np.array(ms);dms=np.array(dms);w=1/dms**2

import sys; print = lambda *a,**kw: __builtins__.print(*a,**kw,flush=True)
np.random.seed(42)
zg=np.unique(np.concatenate([np.linspace(0,0.005,25),
    np.logspace(np.log10(0.006),np.log10(2000),500)]))
CHI2_CACHE={}

def chi2_total(params):
    key=tuple(params)
    if key in CHI2_CACHE: return CHI2_CACHE[key]
    H0,Om,eps,zc=params
    if not all([55<H0<80,0.15<Om<0.45,-0.5<eps<0.5,0.15<zc<1.5]):return 1e30
    Ol=1-Om;h=H0/100;Or=4.15e-5/h**2
    
    # D_C via cumulative integral on grid
    x=np.asarray(zg)/zc; ft=1-eps*x*np.exp(-x)
    E2=Om*(1+zg)**3+Or*(1+zg)**4+Ol*ft
    E=np.sqrt(np.maximum(E2,1e-30))
    DC=np.zeros_like(zg)
    for i in range(1,len(zg)):
        dz=zg[i]-zg[i-1]
        DC[i]=DC[i-1]+0.5*(1/E[i]+1/E[i-1])*dz
    DC*=C/H0
    DCI=ip.interp1d(zg,DC,kind='linear',bounds_error=False,fill_value=0.0)
    Ef=lambda z:np.sqrt(Om*(1+z)**3+Or*(1+z)**4+Ol*(1-eps*(np.asarray(z)/zc)*np.exp(-np.asarray(z)/zc)))
    
    # H(z)
    hz_pred=H0*Ef(HZ[:,0])
    ch=np.sum(((hz_pred-HZ[:,1])/HZ[:,2])**2)
    
    # BAO
    pr=np.zeros(12)
    for i in range(12):
        z=BZ[i];t=BT[i]
        dc=float(DCI(z));dh=C/(H0*Ef(z))
        if t=='DM':pr[i]=dc/RB
        elif t=='DH':pr[i]=dh/RB
        elif t=='DV':pr[i]=(C*z*dc**2/(H0*Ef(z)))**(1/3)/RB
    cb=(pr-BM)@BIC@(pr-BM)
    
    # CMB
    DCR=float(DCI(ZR))
    Rv=np.sqrt(Om)*DCR*H0/C
    la=np.pi*DCR/RC
    vr=PE**2;vla=PLE**2;crl=PRC*PE*PLE
    det=vr*vla-crl**2; dr=Rv-PR; dla=la-PL
    cc=(vla/det)*dr**2+2*(-crl/det)*dr*dla+(vr/det)*dla**2
    
    # SNe
    DCA=DCI(zs)
    DL=(1+zs)*DCA; mu=5*np.log10(DL)+25
    resid=mu-ms; Mb=-np.sum(w*resid)/np.sum(w)
    cs=np.sum(w*(resid+Mb)**2)
    
    total=ch+cs+cb+cc
    CHI2_CACHE[key]=total
    return total

def ln_like(params):
    c=chi2_total(params)
    if c>1e20:return -np.inf
    return -0.5*c

# Fast inline accuracy check
print("Verifying χ² at grid best-fit (H₀=69, Ω_m=0.30, ε=0.25, z_c=0.70)...")
c_test=chi2_total([69.0,0.300,0.25,0.70])
print(f"  χ²_fast = {c_test:.1f}")

# Also verify at MCMC MAP-like point
c_test2=chi2_total([68.3,0.306,0.18,0.40])
print(f"  χ²_fast(zc=0.4) = {c_test2:.1f}")

# MCMC
print(f"\nMCMC: 16 walkers × 6000 steps (1000 burn + 5000 prod)")
nw,bn,ns=24,500,3000
p0=np.array([69.0,0.300,0.25,0.70])+np.array([2,0.02,0.05,0.10])*np.random.randn(nw,4)
t0=time.time()
sampler=emcee.EnsembleSampler(nw,4,ln_like)
sampler.run_mcmc(p0,bn+ns,progress=True)
print(f"Elapsed: {time.time()-t0:.0f}s")

chain=sampler.get_chain(discard=bn,thin=5,flat=True)
lnp=sampler.get_log_prob(discard=bn,thin=5,flat=True)
chi2_chain=np.array([chi2_total(list(p)) for p in chain])

ib=np.argmin(chi2_chain)
mp=chain[ib]; mc2=chi2_chain[ib]

print(f"\n{'='*60}")
print(f"  RESULTS")
print(f"{'='*60}")
print(f"  MAP: H₀={mp[0]:.1f}, Ω_m={mp[1]:.4f}, ε={mp[2]:.4f}, z_c={mp[3]:.4f}")
print(f"  χ²_min = {mc2:.1f}")

# Results (already saved in chain)
print(f"\n  68% CI:")
names=['H₀','Ω_m','ε','z_c']
for i,nm in enumerate(names):
    p=np.percentile(chain[:,i],[16,50,84])
    print(f"  {nm:4s} = {p[1]:.4f} +{p[2]-p[1]:.4f} -{p[1]-p[0]:.4f}")

# Corner plot
import corner
fig=corner.corner(chain,labels=['$H_0$','$\\Omega_m$','$\\varepsilon$','$z_c$'],
    quantiles=[0.16,0.5,0.84],show_titles=True,title_fmt='.3f',
    levels=(0.68,0.95),color='#1a5276',smooth=True,bins=50)
fig.savefig('/home/wsl/cosmic_validation/mcmc_final.png',dpi=300)
print(f"\nSaved: mcmc_final.png")

np.savez('/home/wsl/cosmic_validation/mcmc_final.npz',chain=chain,chi2=chi2_chain,map=mp)
print(f"Saved: mcmc_final.npz")
print(f"\nDone!")