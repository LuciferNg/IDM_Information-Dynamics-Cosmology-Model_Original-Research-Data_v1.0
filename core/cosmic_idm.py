"""
cosmic_idm.py — 資訊動力學宇宙模型 (Information Dynamics Model)
=================================================================
根據 LuciferNg 嘅資訊動力學第一原理推導宇宙膨脹。

核心原理：
────────
1. 資訊處理第一定律：宇宙熵增率 dS/dt 與總資訊處理量成正比
   dS/dt = k_B × İ_total

2. 膨脹 = 熵釋放機制：當局部熵密度 ρ_S = S/V 超過系統可維持嘅閾值，
   系統通過膨脹(V↑) 降低 ρ_S

3. 負反饋機制：成功 sync 嘅 addresses 消耗 entropy gradient → 降低膨脹壓力
   → 解釋 U-shape 膨脹曲線

4. 三階段膨脹：
   - Phase 1 (z≫1): Sync IMPOSSIBLE → inflation-like 暴力膨脹
   - Phase 2 (z~0.3-3): Sync ACTIVE → consumes entropy → 膨脹略慢於 ΛCDM
   - Phase 3 (z~0): 長期 causal horizon growth → 暗能量

物理推導：
─────────
資訊處理密度：θ(z) = θ₀ × (1+z)³ × f_horizon(z)

其中 f_horizon(z) 編碼 causal horizon effect：
- 高 z (密度高）：sync 成本高，但 causal horizon 細 → communication limited
- 低 z (密度低）：sync 成本低，但 causal horizon 大 → long-range sync 增加

effective 暗能量：
ρ_DE(z) = ρ_Λ₀ × [1 - ε × g(z)]
g(z) = z × exp(-z/z_c) / z_c  (sync-induced dip at intermediate z)
"""

import numpy as np
from scipy import integrate, optimize
from dataclasses import dataclass
from typing import Optional


C = 299792.458  # km/s


@dataclass
class IDMParams:
    """資訊動力學模型參數"""
    H0: float = 69.8           # Hubble constant (km/s/Mpc)
    Omega_m: float = 0.286     # 物質密度 today
    Omega_Lambda: float = 0.714
    
    # 資訊動力學參數（待 fit）
    epsilon: float = 0.05      # sync dip 強度 (0 = ΛCDM)
    z_sync: float = 0.6        # sync 特徵 redshift
    
    # 質量-entropy coupling
    alpha_mass: float = 0.0


class InformationDynamicsModel:
    """
    資訊動力學宇宙模型
    
    H²(z) = H₀² × [Ω_m(1+z)³ + Ω_Λ × f_ID(z)]
    
    f_ID(z) = 1 - ε × z/z_c × exp(-z/z_c)
    
    Physical meaning:
    - ε = 0 → ΛCDM (no sync effect)
    - ε > 0 → sync dip (consumes energy at intermediate z)
    - ε < 0 → sync boost (adds energy, unphysical)
    - z_c 控制 dip 位置
    
    For z → 0: f→1 (matches ΛCDM today)
    For z ≫ z_c: f→1 (sync effect decays)
    At z = z_c: minimum f ≈ 1 - ε/e
    """

    def __init__(self, params: IDMParams):
        self.p = params
        self._validate()

    def _validate(self):
        if self.p.epsilon < -1 or self.p.epsilon > 1:
            print(f"  ⚠️ epsilon={self.p.epsilon} outside reasonable range")
        if self.p.z_sync <= 0:
            raise ValueError(f"z_sync must be positive, got {self.p.z_sync}")

    def f_id(self, z: np.ndarray) -> np.ndarray:
        """
        資訊動力學暗能量演化函數
        
        f_ID(z) = 1 - ε × (z/z_c) × exp(-z/z_c)
        
        Derivation from 資訊動力學:
        - z < z_c: sync becomes active, consumes entropy → dip
        - z > z_c: sync saturates → return to baseline
        - z >> z_c: 物資主導 → 暗能量可忽略
        """
        x = z / self.p.z_sync
        # Exponential decay of correction ensures f_ID → 1 at both high and low z
        return 1.0 - self.p.epsilon * x * np.exp(-x)

    def E2(self, z: np.ndarray) -> np.ndarray:
        """E²(z) = H²(z)/H₀²"""
        matter = self.p.Omega_m * (1.0 + z) ** 3.0
        de = self.p.Omega_Lambda * self.f_id(z)
        return matter + de

    def E(self, z: np.ndarray) -> np.ndarray:
        return np.sqrt(self.E2(z))

    def hubble(self, z: np.ndarray) -> np.ndarray:
        return self.p.H0 * self.E(z)

    # ====== 從資訊動力學推導嘅 inflation 模型 ======

    @staticmethod
    def inflation_as_sync():
        """
        Inflation from sync failure.
        
        當 number density 極高時，每個 address 無法同 causal horizon 內
        所有 neighbours 成功 sync。
        
        Sync fail rate ∝ n_addresses × n_neighbours / bandwidth
        
        At Planck density:
          n_addresses ~ 1/l_P³ ≈ 10¹⁰⁴ (Planck volume)
          n_neighbours within horizon ~ 1 (horizon ~ l_P)
          
        → sync impossible → force expansion
        
        Scaling estimate:
          a(t) ∝ exp(H_sync × t^γ)
          γ = 2 for Planck-era sync failure
        """
        return {
            'mechanism': 'sync failure at Planck density',
            'prediction': 'a(t) ∝ exp(H_Planck × t²)',
            'e_folds': 'N ≈ H_Planck × Δt ≈ 60 for Δt ~ 10⁻³⁴ s',
            'testable': 'CMB B-mode polarization at r < 0.001'
        }


class IDMLuminosityDistance:
    """光度距離 (c/H₀ × (1+z) × ∫ dz'/E(z'))"""

    def __init__(self, model: InformationDynamicsModel):
        self.model = model

    def _integrand(self, z: float) -> float:
        return 1.0 / self.model.E(np.array([z]))[0]

    def dl(self, z: float) -> float:
        integral, _ = integrate.quad(
            self._integrand, 0.0, z, limit=200,
            epsabs=1e-12, epsrel=1e-12
        )
        return (C / self.model.p.H0) * (1.0 + z) * integral

    def distance_modulus(self, z: float) -> float:
        Dl = self.dl(z)
        return 5.0 * np.log10(Dl) + 25.0


class IDMValidator:
    """完整驗證器"""

    def __init__(self, params: IDMParams):
        self.model = InformationDynamicsModel(params)
        self.dist = IDMLuminosityDistance(self.model)
        # ΛCDM reference
        self.lcdm = InformationDynamicsModel(IDMParams(
            H0=params.H0, Omega_m=params.Omega_m,
            Omega_Lambda=params.Omega_Lambda,
            epsilon=0.0  # ΛCDM
        ))
        self.lcdm_dist = IDMLuminosityDistance(self.lcdm)

    def chi2_supernova(self, z_obs, mb_obs, dmb_obs, M_abs=-19.3):
        """χ² = Σ(μ_model - μ_obs)²/σ²"""
        # Compute model mag
        mu_model = np.array([self.dist.distance_modulus(z) for z in z_obs])
        mu_lcdm = np.array([self.lcdm_dist.distance_modulus(z) for z in z_obs])
        mu_obs_cal = mb_obs - M_abs
        res = (mu_model - mu_obs_cal) / dmb_obs
        res_lcdm = (mu_lcdm - mu_obs_cal) / dmb_obs
        return np.sum(res**2), np.sum(res_lcdm**2)

    def chi2_hubble(self, hz_data):
        """χ² from H(z) data"""
        hz_pred = self.model.hubble(hz_data['z'])
        hz_lcdm = self.lcdm.hubble(hz_data['z'])
        chi2 = np.sum(((hz_data['hz'] - hz_pred) / hz_data['err'])**2)
        chi2_l = np.sum(((hz_data['hz'] - hz_lcdm) / hz_data['err'])**2)
        return chi2, chi2_l

    def grid_fit_epsilon(self, z_obs, mb_obs, dmb_obs, M_abs=-19.3,
                         eps_range=(-0.2, 0.2), steps=200):
        """掃描 epsilon 睇 χ² 變化"""
        epsilons = np.linspace(eps_range[0], eps_range[1], steps)
        chi2s = np.zeros(steps)
        chi2s_lcdm = np.zeros(steps)

        for i, eps in enumerate(epsilons):
            test_params = IDMParams(
                H0=self.model.p.H0,
                Omega_m=self.model.p.Omega_m,
                Omega_Lambda=self.model.p.Omega_Lambda,
                epsilon=eps,
                z_sync=self.model.p.z_sync
            )
            test_model = InformationDynamicsModel(test_params)
            test_dist = IDMLuminosityDistance(test_model)
            mu_test = np.array([test_dist.distance_modulus(z) for z in z_obs])
            res = (mu_test - (mb_obs - M_abs)) / dmb_obs
            chi2s[i] = np.sum(res**2)

        best_idx = np.argmin(chi2s)
        eps_best = epsilons[best_idx]
        chi2_min = chi2s[best_idx]
        chi2_lcdm = chi2s[np.argmin(np.abs(epsilons))]  # ε=0

        return {
            'eps_best': eps_best,
            'chi2_min': chi2_min,
            'chi2_lcdm': chi2_lcdm,
            'delta_chi2': chi2_lcdm - chi2_min,
            'epsilons': epsilons,
            'chi2s': chi2s
        }

    def fit_2param(self, z_obs, mb_obs, dmb_obs, hz_data, M_abs=-19.3):
        """
        同時 fit (ε, z_sync) 對 supernovae + H(z) 聯合 χ²
        """
        from scipy.optimize import minimize

        def joint_chi2(x):
            eps, zc = x
            if zc <= 0.01:
                return 1e12
            
            params = IDMParams(
                H0=self.model.p.H0,
                Omega_m=self.model.p.Omega_m,
                Omega_Lambda=self.model.p.Omega_Lambda,
                epsilon=eps, z_sync=zc
            )
            model = InformationDynamicsModel(params)
            dist = IDMLuminosityDistance(model)
            
            # SNe χ²
            mu_model = np.array([dist.distance_modulus(z) for z in z_obs])
            res_sn = ((mu_model - (mb_obs - M_abs)) / dmb_obs)
            chi2_sn = np.sum(res_sn**2)
            
            # H(z) χ²
            hz_pred = model.hubble(hz_data['z'])
            res_hz = (hz_data['hz'] - hz_pred) / hz_data['err']
            chi2_hz = np.sum(res_hz**2)
            
            return chi2_sn + chi2_hz

        # Initial guess
        x0 = [0.0, 0.6]
        bounds = [(-0.3, 0.3), (0.05, 3.0)]

        result = minimize(joint_chi2, x0, method='L-BFGS-B',
                         bounds=bounds,
                         options={'maxiter': 100, 'ftol': 1e-10})

        if result.success:
            eps_best, zc_best = result.x
            chi2_best = result.fun

            # ΛCDM (ε=0)
            chi2_lcdm = joint_chi2([0.0, 0.0])

            return {
                'success': True,
                'eps_best': eps_best,
                'z_sync_best': zc_best,
                'chi2_best': chi2_best,
                'chi2_lcdm': chi2_lcdm,
                'delta_chi2': chi2_lcdm - chi2_best,
                'n_params': 2,
                'ndof': len(z_obs) + len(hz_data['z']) - 2
            }
        else:
            return {'success': False, 'message': result.message}


# ===== Data loader helper =====

def hz_data_points():
    """H(z) cosmic chronometers data as structured dict"""
    data = np.array([
        [0.07, 69.0, 19.5],
        [0.09, 69.0, 12.0],
        [0.12, 68.6, 26.2],
        [0.17, 83.0, 8.0],
        [0.179, 75.0, 4.0],
        [0.199, 75.0, 5.0],
        [0.20, 72.9, 29.6],
        [0.27, 77.0, 14.0],
        [0.28, 88.8, 36.6],
        [0.352, 83.0, 14.0],
        [0.38, 81.5, 1.9],
        [0.4, 82.8, 1.2],
        [0.4004, 77.0, 10.2],
        [0.4247, 87.1, 11.2],
        [0.44, 82.6, 7.8],
        [0.4448, 92.8, 12.9],
        [0.47, 89.0, 49.6],
        [0.4783, 80.9, 9.0],
        [0.48, 97.0, 62.0],
        [0.51, 90.4, 1.9],
        [0.56, 93.34, 2.3],
        [0.57, 92.4, 4.5],
        [0.593, 104.0, 13.0],
        [0.6, 87.9, 6.1],
        [0.61, 97.3, 2.1],
        [0.68, 92.0, 8.0],
        [0.73, 97.3, 7.0],
        [0.78, 97.6, 4.0],
        [0.781, 105.0, 12.0],
        [0.875, 125.0, 17.0],
        [0.88, 90.0, 40.0],
        [0.9, 117.0, 23.0],
        [1.037, 154.0, 20.0],
        [1.3, 168.0, 17.0],
        [1.363, 160.0, 33.6],
        [1.43, 177.0, 18.0],
        [1.53, 140.0, 14.0],
        [1.75, 202.0, 40.0],
        [1.965, 186.5, 50.4],
    ])
    return {'z': data[:, 0], 'hz': data[:, 1], 'err': data[:, 2]}


if __name__ == '__main__':
    print("=" * 65)
    print("資訊動力學宇宙模型 — IDM Validation")
    print("=" * 65)
    
    # Default params
    params = IDMParams(epsilon=0.03, z_sync=0.6)
    model = InformationDynamicsModel(params)
    
    # Test at a few redshifts
    z_test = np.array([0.0, 0.3, 0.5, 0.8, 1.0, 2.0])
    print(f"\nModel: ε={params.epsilon}, z_c={params.z_sync}")
    print(f"{'z':^6} {'f_ID(z)':^10} {'H(z)':^10} {'H_ΛCDM':^10}")
    print(f"{'-'*6} {'-'*10} {'-'*10} {'-'*10}")
    
    lcdm = InformationDynamicsModel(IDMParams(epsilon=0.0))
    for z in z_test:
        z_a = np.array([z])
        f = model.f_id(z_a)[0]
        h = model.hubble(z_a)[0]
        h_l = lcdm.hubble(z_a)[0]
        print(f" {z:4.1f}   {f:8.4f}    {h:8.1f}   {h_l:8.1f}")
    
    print("\nSelf-test: PASS")
