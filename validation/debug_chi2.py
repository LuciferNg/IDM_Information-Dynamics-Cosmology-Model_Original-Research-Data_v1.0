"""debug_chi2.py — Debug χ² discrepancy in fast interpolation"""
import numpy as np; from scipy import integrate, interpolate as ip

C=299792.458; ZR=1089.92; RB=147.09; RC=144.43
PR=1.7502; PE=0.0046; PL=301.808; PLE=0.090; PRC=0.47

HZ=np.array([[0.38,81.5,1.9],[0.51,90.4,1.9]])
BM=np.array([7.925,13.620,20.983,16.846,20.079,21.708,17.876,27.787,13.824,26.072,39.708,8.523])
BZ=np.array([0.295,0.510,0.510,0.706,0.706,0.930,0.930,1.317,1.317,1.491,2.330,2.330])
BT=['DV','DM','DH','DM','DH','DM','DH','DM','DH','DV','DM','DH']

H0,Om,eps,zc=69.0,0.300,0.25,0.70
Ol=1-Om;h=H0/100;Or=4.15e-5/h**2;Ef=lambda z:np.sqrt(Om*(1+z)**3+Or*(1+z)**4+Ol*(1-eps*(z/zc)*np.exp(-z/zc)))

# D_C via quad (reference)
def DC_quad(z):
    f=lambda x:1/Ef(x)
    I,_=integrate.quad(f,0,z,limit=2000,epsabs=1e-12,epsrel=1e-12)
    return C/H0*I

# D_C via interpolation (MCMC method)
zg=np.unique(np.concatenate([np.linspace(0,0.005,25),np.logspace(np.log10(0.006),np.log10(2000),500)]))
x=zg/zc;ft=1-eps*x*np.exp(-x)
E2=Om*(1+zg)**3+Or*(1+zg)**4+Ol*ft
E=np.sqrt(np.maximum(E2,1e-30))
DC=np.zeros_like(zg)
for i in range(1,len(zg)):
    dz=zg[i]-zg[i-1]
    DC[i]=DC[i-1]+0.5*(1/E[i]+1/E[i-1])*dz
DC*=C/H0
DCI=ip.interp1d(zg,DC,kind='linear',bounds_error=False,fill_value=0.0)

# Compare components
print("Component comparison (quad vs interp):")

# H(z)
hz_q=H0*Ef(HZ[:,0]); hz_i=hz_q  # same
print(f"H(z): same")

# BAO
print(f"\nBAO:")
for i,(z,t) in enumerate(zip(BZ,BT)):
    dc_q=DC_quad(z); dc_i=float(DCI(z))
    h_q=C/(H0*Ef(z)); h_i=h_q
    if t=='DM': v_q=dc_q/RB; v_i=dc_i/RB
    elif t=='DH': v_q=h_q/RB; v_i=h_i/RB
    elif t=='DV': v_q=(C*z*dc_q**2/h_q)**(1/3)/RB; v_i=(C*z*dc_i**2/h_i)**(1/3)/RB
    print(f"  z={z:.3f} {t}: quad={v_q:.4f} interp={v_i:.4f} diff={v_q-v_i:.6f} data={BM[i]:.4f}")

# Find where the chi2 divergence comes from
print(f"\nD_C comparison at SN redshifts:")
zs,ms,dms=[],[],[]
with open('/tmp/pantheon_raw.txt') as f:
    for l in f:
        if l.startswith("#")or not l.strip():continue
        p=l.strip().split()
        if len(p)<6:continue
        z=float(p[1])
        if z>=0.01:zs.append(z);ms.append(float(p[4]));dms.append(float(p[5]))
zs=np.array(zs);ms=np.array(ms);dms=np.array(dms)

# Check first 10 SNe
for i in range(10):
    z=zs[i]
    dc_q=DC_quad(z); dc_i=float(DCI(z))
    mu_q=5*np.log10((1+z)*dc_q)+25
    mu_i=5*np.log10((1+z)*dc_i)+25
    print(f"  z={z:.4f}: DC_q={dc_q:.1f} DC_i={dc_i:.1f} μ_q={mu_q:.4f} μ_i={mu_i:.4f} data={ms[i]:.4f}")

# Now compute full chi2 from interpolation code vs direct
print(f"\nFull chi2 with BAO_ICOV...")
BC=np.array([[2.2723e-02,0,0,0,0,0,0,0,0,0,0,0],
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
  [0,0,0,0,0,0,0,0,0,0,-7.6948e-02,2.9186e-02]])
BIC=np.linalg.inv(BC)
pr=np.zeros(12)
for i in range(12):
    z=BZ[i];t=BT[i]
    dc=float(DCI(z));dh=C/(H0*Ef(z))
    if t=='DM':pr[i]=dc/RB
    elif t=='DH':pr[i]=dh/RB
    elif t=='DV':pr[i]=(C*z*dc**2/dh)**(1/3)/RB
cb=(pr-BM)@BIC@(pr-BM)
print(f"χ²_BAO (interp) = {cb:.1f}")

# Quad-based BAO chi2
pr2=np.zeros(12)
for i in range(12):
    z=BZ[i];t=BT[i]
    dc=DC_quad(z);dh=C/(H0*Ef(z))
    if t=='DM':pr2[i]=dc/RB
    elif t=='DH':pr2[i]=dh/RB
    elif t=='DV':pr2[i]=(C*z*dc**2/dh)**(1/3)/RB
cb2=(pr2-BM)@BIC@(pr2-BM)
print(f"χ²_BAO (quad)  = {cb2:.1f}")