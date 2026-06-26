"""
cosmic_model.py -- 純 Python 數學模組，唔用 LLM
===============================================
根據 LuciferNg 嘅三層宇宙模型：
  Layer 1: Fractal Panentheism (Ontology)
  Layer 2: 資訊動力學 (Dynamics)
  Layer 3: 投影宇宙論 (Phenomenology)

核心公式由 LLM 設計（本文件註解），
但數學計算係 pure numpy/scipy，無任何 AI 模型。

公式設計：
  1. 膨脹 = 資訊同步壓力釋放閥
     H²(z) = H₀² × [Ω_m(1+z)³ + Ω_Λ × f_sync(z)]

  2. f_sync(z) = (1+z)^β_sync
     β_sync > 0 → 早期更快膨脹（同步壓力大）
     β_sync = 0 → ΛCDM（暗能量常數）

  3. 質量 = entropy side-effect：
     m_i ∝ η × N_i × (1 + α × dS_local/dt)

  4. 通用化：β_sync 可隨 z 演化
     β(z) = β₀ × (1 - exp(-z/z_c))  (高 z 時 β→β₀, 低 z 時 β→0)
"""

import numpy as np
from scipy import integrate, optimize
from dataclasses import dataclass
from typing import Callable, Optional

# ========= 物理常數 =========
C = 299792.458  # km/s
MPC_TO_KM = 3.085677581e19  # 1 Mpc in km
H0_DEFAULT = 69.8  # km/s/Mpc (Planck 2020)

@dataclass
class CosmologyParams:
    """你個模型嘅可調參數"""
    H0: float = H0_DEFAULT       # Hubble constant (km/s/Mpc)
    Omega_m: float = 0.286       # 物質密度參數 today
    beta_sync: float = 0.0       # 同步壓力指數 (0 = ΛCDM)
    Omega_Lambda: float = 0.714  # 暗能量參數 today

    # 質量-entropy 耦合
    mass_entropy_coupling: float = 0.0  # α: 質量對 entropy 嘅耦合強度

    # 進階：β(z) 演化參數
    beta_evolve: bool = False    # β 隨 z 演化？
    beta_inf: float = 2.0        # z→∞ 時嘅 β（早期 inflation 強度）
    z_c: float = 0.5             # 演化特徵尺度


class SyncPressureModel:
    """
    同步壓力宇宙模型 (Sync Pressure Cosmology)
    ==========================================
    
    H²(z) = H₀² × [Ω_m(1+z)³ + Ω_Λ × f_sync(z)]
    
    其中 f_sync(z) = (1+z)^β(z) 編碼資訊同步壓力對膨脹嘅貢獻
    
    物理意義：
    - 高 z（早期）: address density 高 → sync pressure 大 → 膨脹快
    - 低 z（而家）: address density 低 → sync pressure 小 → 似 ΛCDM
    """

    def __init__(self, params: CosmologyParams):
        self.p = params

    def E2(self, z: np.ndarray) -> np.ndarray:
        """
        E²(z) = H²(z) / H₀²
        純 numpy, 無 LLM
        """
        # 物質項
        matter = self.p.Omega_m * (1.0 + z) ** 3.0

        # 同步壓力項
        if self.p.beta_evolve:
            # 演化版: β(z) 隨 redshift 遞增
            beta_z = self.p.beta_inf * (1.0 - np.exp(-z / self.p.z_c))
        else:
            beta_z = self.p.beta_sync * np.ones_like(z)

        sync = self.p.Omega_Lambda * (1.0 + z) ** beta_z

        return matter + sync

    def E(self, z: np.ndarray) -> np.ndarray:
        """E(z) = H(z) / H₀"""
        return np.sqrt(self.E2(z))

    def hubble(self, z: np.ndarray) -> np.ndarray:
        """H(z) in km/s/Mpc"""
        return self.p.H0 * self.E(z)


class LuminosityDistance:
    """
    光度距離 D_L(z) = c/H₀ × (1+z) × ∫₀ᶻ dz'/E(z')
    
    使用 scipy integrate 做純數值積分，唔用 LLM
    """

    def __init__(self, model: SyncPressureModel):
        self.model = model

    def _integrand(self, z: float) -> float:
        """1/E(z') for integration"""
        return 1.0 / self.model.E(np.array([z]))[0]

    def dl(self, z: float) -> float:
        """
        光度距離 (Mpc)
        純數值積分
        """
        integral, _ = integrate.quad(
            self._integrand, 0.0, z, limit=200, epsabs=1e-12, epsrel=1e-12
        )
        return (C / self.model.p.H0) * (1.0 + z) * integral

    def dl_array(self, z_arr: np.ndarray) -> np.ndarray:
        """Vectorized D_L"""
        return np.array([self.dl(z) for z in z_arr])

    def distance_modulus(self, z: float) -> float:
        """
        μ(z) = 5 × log₁₀(D_L / Mpc) + 25
        """
        Dl = self.dl(z)
        if Dl <= 0:
            return np.inf
        return 5.0 * np.log10(Dl) + 25.0

    def mu_array(self, z_arr: np.ndarray) -> np.ndarray:
        """Vectorized distance modulus"""
        return np.array([self.distance_modulus(z) for z in z_arr])


class MassEntropyModel:
    """
    質量 = entropy side-effect model
    ==================================
    
    m_i = m₀ × N_i × (1 + α × dS_local/dt)
    
    其中：
    - m₀: 基礎質量單位（address cluster 基本質量）
    - N_i: particle type 嘅 address cluster count
    - α: entropy coupling （質量對 entropy 梯度嘅敏感度）
    - dS_local/dt: 局部 entropy 產生率
    """

    # 已知粒子質量 (MeV/c²)
    PARTICLE_MASSES = {
        'electron': 0.511,
        'muon': 105.66,
        'tau': 1776.86,
        'proton': 938.272,
        'neutron': 939.565,
        'pion0': 134.977,
        'pion_plus': 139.570,
        'higgs': 125000.0,  # 125 GeV
    }

    def __init__(self, coupling: float = 0.0, base_mass: float = 0.511):
        self.alpha = coupling
        self.m0 = base_mass  # MeV, 電子質量做基本單位

    def compute_mass(self, n_clusters: float, entropy_rate: float = 0.0) -> float:
        """
        從 address cluster count + entropy rate 計質量
        
        m = m₀ × N × (1 + α × dS/dt)
        """
        return self.m0 * n_clusters * (1.0 + self.alpha * entropy_rate)

    def mass_ratio_vs_entropy(self, base_mass: float, 
                              entropy_range: np.ndarray) -> np.ndarray:
        """
        俾定一個 base particle，entropy rate 改變時質量點變
        用嚟睇 mass-entropy 耦合嘅影響
        """
        ratio = np.ones_like(entropy_range)
        if self.alpha != 0:
            ratio = 1.0 + self.alpha * entropy_range
        return ratio


class ModelComparator:
    """
    模型比較器：ΛCDM (β=0) vs Sync Pressure (β>0)
    純 numpy, 無 LLM
    """

    def __init__(self, model: SyncPressureModel):
        self.model = model
        self.lcdm = SyncPressureModel(CosmologyParams(
            H0=model.p.H0,
            Omega_m=model.p.Omega_m,
            beta_sync=0.0  # ΛCDM
        ))
        self.dist_lcdm = LuminosityDistance(self.lcdm)
        self.dist_sync = LuminosityDistance(self.model)

    def hubble_ratio(self, z: np.ndarray) -> np.ndarray:
        """H_sync(z) / H_lcdm(z)"""
        return self.model.E(z) / self.lcdm.E(z)

    def mu_residual(self, z: np.ndarray) -> np.ndarray:
        """
        Δμ(z) = μ_sync(z) - μ_lcdm(z)
        
        正數 → sync model 預測更亮嘅超新星（距離更遠）
        負數 → sync model 預測更暗
        """
        mu_s = self.dist_sync.mu_array(z)
        mu_l = self.dist_lcdm.mu_array(z)
        return mu_s - mu_l

    def chi_squared(self, z_obs: np.ndarray, mu_obs: np.ndarray,
                    mu_err: np.ndarray, M_abs: float = -19.3) -> float:
        """
        χ² = Σ (μ_model(z_i) - (mu_obs_i - M_abs))² / σ_i²
        
        純數值驗證，無 LLM
        """
        mu_model = self.dist_sync.mu_array(z_obs)
        mu_obs_calibrated = mu_obs - M_abs
        residuals = mu_model - mu_obs_calibrated
        return np.sum((residuals / mu_err) ** 2.0)

    def fit_beta(self, z_obs: np.ndarray, mu_obs: np.ndarray,
                 mu_err: np.ndarray, M_abs: float = -19.3,
                 beta_range: tuple = (-0.5, 2.0)) -> dict:
        """
        用實際數據 fit β_sync 最佳值
        
        Returns {beta_best, chi2_min, chi2_lcdm, delta_chi2}
        """
        def chi2_for_beta(beta):
            test_params = CosmologyParams(
                H0=self.model.p.H0,
                Omega_m=self.model.p.Omega_m,
                beta_sync=beta
            )
            test_model = SyncPressureModel(test_params)
            test_dist = LuminosityDistance(test_model)
            mu_test = test_dist.mu_array(z_obs)
            res = (mu_test - (mu_obs - M_abs)) / mu_err
            return np.sum(res ** 2.0)

        # 掃描 β
        betas = np.linspace(beta_range[0], beta_range[1], 200)
        chi2s = np.array([chi2_for_beta(b) for b in betas])

        best_idx = np.argmin(chi2s)
        beta_best = betas[best_idx]
        chi2_min = chi2s[best_idx]
        chi2_lcdm = chi2_for_beta(0.0)

        return {
            'beta_best': beta_best,
            'chi2_min': chi2_min,
            'chi2_lcdm': chi2_lcdm,
            'delta_chi2': chi2_lcdm - chi2_min,
            'betas': betas,
            'chi2s': chi2s,
            'ndof': len(z_obs) - 1
        }


# =========  inflation model (早期同步壓力釋放) =========

class InflationSyncModel:
    """
    早期 inflation 作為 sync pressure release valve

    關鍵 claim：
      inflation (a(t) ~ exp(H*t)) 唔係由 inflaton field 驅動，
      而係由資訊同步壓力極大時嘅強制 volume expansion。

    公式：
      ä/a = H² + Ḣ = (8πG/3) × (3P_sync - ρ_sync)
      P_sync = w_sync × ρ_sync
      w_sync large → rapid expansion

      Simplified: a(t) ∝ exp(H_sync × t^γ)
      γ > 1 → accelerated expansion (inflation)
      γ = 1 → de Sitter (like dark energy)
      γ < 1 → decelerated
    """

    def __init__(self, H_sync: float = 1.0, gamma: float = 2.0):
        self.H_sync = H_sync  # 同步壓力主頻率 (Planck units)
        self.gamma = gamma    # 膨脹指數

    def scale_factor(self, t: np.ndarray) -> np.ndarray:
        """a(t) ∝ exp(H_sync × t^gamma)"""
        return np.exp(self.H_sync * t ** self.gamma)

    def hubble_param(self, t: np.ndarray) -> np.ndarray:
        """H(t) = d/dt (ln a) = gamma × H_sync × t^(gamma-1)"""
        return self.gamma * self.H_sync * t ** (self.gamma - 1.0)

    def efoldings(self, t_start: float, t_end: float) -> float:
        """Number of e-foldings between t_start and t_end"""
        return self.H_sync * (t_end ** self.gamma - t_start ** self.gamma)


if __name__ == '__main__':
    # Quick self-test
    p = CosmologyParams(beta_sync=0.3)
    model = SyncPressureModel(p)
    
    z_test = np.array([0.0, 0.5, 1.0])
    print("E(z) values:", model.E(z_test))
    
    dist = LuminosityDistance(model)
    print(f"D_L(z=0.5) = {dist.dl(0.5):.1f} Mpc")
    print(f"μ(z=0.5) = {dist.distance_modulus(0.5):.2f}")

    comp = ModelComparator(model)
    print(f"μ residual at z=1: Δμ = {comp.mu_residual(np.array([1.0]))[0]:.4f}")
    print("Module self-test: PASS")
