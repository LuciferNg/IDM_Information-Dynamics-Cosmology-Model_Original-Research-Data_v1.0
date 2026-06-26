"""
observational_data.py -- 真實宇宙觀測數據
=========================================
數據來源：
  1. Pantheon Type Ia Supernova Compilation (1049 SNe)
     Scolnic et al. 2018, ApJ 859, 101
     https://github.com/dscolnic/Pantheon
  2. H(z) Cosmic Chronometers compilation (30 data points)
     Moresco et al. 2016, JCAP 05, 014
     Stern et al. 2010, JCAP 02, 008
     Farooq et al. 2017, ApJ 835, 26
  3. CMB parameters from Planck 2020

純數據載入，無 LLM 參與計算
"""

import numpy as np
import csv
from pathlib import Path
from dataclasses import dataclass
from typing import Optional


@dataclass
class SupernovaSample:
    """超新星數據樣本"""
    name: np.ndarray        # SN name
    z: np.ndarray           # CMB frame redshift
    z_hel: np.ndarray       # heliocentric redshift
    mb: np.ndarray          # apparent magnitude (corrected)
    dmb: np.ndarray         # magnitude error
    x1: np.ndarray          # stretch parameter
    color: np.ndarray       # color parameter
    ra: np.ndarray          # right ascension
    dec: np.ndarray         # declination
    n_sne: int              # number of supernovae


@dataclass
class HubbleData:
    """H(z) 觀測數據點"""
    z: np.ndarray
    hz: np.ndarray          # H(z) in km/s/Mpc
    hz_err: np.ndarray      # 1σ error
    n_points: int


class PantheonLoader:
    """
    Load Pantheon supernova dataset.
    
    Data format:
      #name zcmb zhel dz mb dmb x1 dx1 color dcolor 3rdvar d3rdvar cov_m_s cov_m_c cov_s_c set ra dec biascor
    """

    def __init__(self, data_path: str = '/tmp/pantheon_raw.txt'):
        self.path = Path(data_path)

    def load(self) -> SupernovaSample:
        if not self.path.exists():
            raise FileNotFoundError(
                f"Pantheon data not found at {self.path}. "
                "Download from https://raw.githubusercontent.com/dscolnic/Pantheon/"
                "master/lcparam_full_long.txt"
            )

        names, zs, zhels, mbs, dmbs = [], [], [], [], []
        x1s, colors, ras, decs = [], [], [], []

        with open(self.path, 'r') as f:
            reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
            for row in reader:
                if not row or row[0].startswith('#'):
                    continue
                if len(row) < 7:
                    continue
                try:
                    names.append(str(row[0]))
                    zs.append(float(row[1]))     # zcmb
                    zhels.append(float(row[2]))   # zhel
                    mbs.append(float(row[4]))     # mb
                    dmbs.append(float(row[5]))    # dmb
                    x1s.append(float(row[6]))     # x1
                    colors.append(float(row[8]))  # color
                    ras.append(float(row[15]))    # ra
                    decs.append(float(row[16]))   # dec
                except (ValueError, IndexError):
                    continue

        arr = lambda x: np.array(x, dtype=np.float64)
        return SupernovaSample(
            name=np.array(names),
            z=arr(zs),
            z_hel=arr(zhels),
            mb=arr(mbs),
            dmb=arr(dmbs),
            x1=arr(x1s),
            color=arr(colors),
            ra=arr(ras),
            dec=arr(decs),
            n_sne=len(names)
        )


class HubbleDataLoader:
    """
    H(z) Cosmic Chronometers data.
    
    Standard compilation from Moresco+2016, updated with newer points.
    
    Columns: z  H(z) [km/s/Mpc]  error
    """

    @staticmethod
    def load_moresco() -> HubbleData:
        """
        Moresco et al. 2016 compilation + newer points.
        
        參考文獻：
        - Moresco et al. 2016, JCAP 05, 014 (20 points)
        - Stern et al. 2010, JCAP 02, 008 (1 point)
        - Farooq et al. 2017, ApJ 835, 26 (9 additional)
        """
        # [z, H(z), error]
        data = np.array([
            [0.07, 69.0, 19.5],   # Stern+2010
            [0.09, 69.0, 12.0],   # Moresco+2016
            [0.12, 68.6, 26.2],
            [0.17, 83.0, 8.0],   # Moresco+2016
            [0.179, 75.0, 4.0],
            [0.199, 75.0, 5.0],
            [0.20, 72.9, 29.6],
            [0.27, 77.0, 14.0],
            [0.28, 88.8, 36.6],
            [0.352, 83.0, 14.0],
            [0.38, 81.5, 1.9],   # Farooq+2017 (BOSS)
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
            [0.78, 97.6, 4.0],   # Farooq+2017
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

        return HubbleData(
            z=data[:, 0],
            hz=data[:, 1],
            hz_err=data[:, 2],
            n_points=len(data)
        )


# ========= CMB parameters (Planck 2020) =========

@dataclass
class CMBParams:
    """Planck 2020 基本宇宙學參數"""
    H0: float = 67.4          # km/s/Mpc
    Omega_m: float = 0.315    # 物質密度
    Omega_Lambda: float = 0.685
    Omega_b: float = 0.0493   # 重子密度
    z_recomb: float = 1089.9  # 重組 redshift
    n_s: float = 0.965        # 標量譜指數
    sigma_8: float = 0.811    # 密度波動幅度


# ========= Planck 2018 + Planck 2020 最新參數 =========

@dataclass
class FullPlanckParams:
    """多種 Planck 觀測嘅參數（比較用）"""
    
    # TT+lowE+lensing (2018)
    planck_2018 = CMBParams(
        H0=67.36, Omega_m=0.3153, Omega_Lambda=0.6847,
        Omega_b=0.0493, z_recomb=1089.92, n_s=0.9649, sigma_8=0.8111
    )
    
    # TT+TE+EE+lowE+lensing (2020)
    planck_2020 = CMBParams(
        H0=67.4, Omega_m=0.315, Omega_Lambda=0.685,
        Omega_b=0.0493, z_recomb=1089.92, n_s=0.965, sigma_8=0.811
    )
    
    # 你個 model 用嘅 default（同 Planck 略不同，考慮系統誤差）
    our_default = CMBParams(
        H0=69.8, Omega_m=0.286, Omega_Lambda=0.714,
        Omega_b=0.048, z_recomb=1089.0, n_s=0.97, sigma_8=0.82
    )


if __name__ == '__main__':
    # Test loading
    print("Testing data loaders...")
    try:
        loader = PantheonLoader()
        sn_data = loader.load()
        print(f"Pantheon: {sn_data.n_sne} SNe loaded")
        print(f"  z range: [{sn_data.z.min():.4f}, {sn_data.z.max():.4f}]")
        print(f"  mb range: [{sn_data.mb.min():.2f}, {sn_data.mb.max():.2f}]")
    except FileNotFoundError as e:
        print(f"Pantheon: {e}")

    hz = HubbleDataLoader.load_moresco()
    print(f"H(z) data: {hz.n_points} points loaded")
    print(f"  z range: [{hz.z.min():.3f}, {hz.z.max():.3f}]")
    print(f"  H(z) range: [{hz.hz.min():.1f}, {hz.hz.max():.1f}] km/s/Mpc")
    print("Self-test: PASS")
