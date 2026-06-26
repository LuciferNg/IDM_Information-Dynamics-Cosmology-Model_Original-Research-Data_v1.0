# 資訊動力學宇宙論（IDM）：探索過程與延伸想法

## — 從第一次 code 到 CLASS Boltzmann solver 嘅演變

---

## 五、驗證旅程：從 H(z) 到六個獨立測試

### 5.1 第一版：cosmic_idm.py（6月24日）

第一個 Python 檔案 `cosmic_idm.py` 係最原始嘅骨架。當時我嘅 parameters 係 H₀=69.8, Ω_m=0.286, ε=0.05, z_sync=0.6——呢啲數值全部係粗略估計。

第一個測試係單一 dataset：39 個宇宙計時器（Cosmic Chronometers）H(z) 數據點。結果 IDM χ²=24.30 vs ΛCDM χ²=30.14——改善 Δχ²=+5.84。

**呢個 5.84 嘅改善係第一個 signal：dip 存在。**

但呢個時候我仲唔敢信。39 點數據，改善 5.84 點 χ²，可能係 noise。我需要更多獨立測試。

### 5.2 BUG 大災難：聲音交錯的混亂期

呢個 project 最痛苦嘅時期係 6月24-25日——LLM API 頻繁故障，context window 不斷 compaction，導致：

1. **Sound horizon swap bug** — CMB l_A 要用 r*=144.43 Mpc（重組期聲學視界），BAO 要用 r_d=147.09 Mpc（drag epoch）。但我嘅 code 初期用錯咗 r_d 去計 l_A，結果 χ²_CMB 大咗 ~2000 點，ε 出現負值（−0.26）——完全違反 CLT 預測嘅 dip（ε>0）。呢個 bug 浪費咗 2+ hours chasing phantom discrepancy。

2. **缺少 radiation term Ω_r(1+z)⁴** — 積分到 recombination (z≈1090) 時，radiation 佔 ~38% 嘅 total energy density。冇加 radiation term 令 D_C(z_rec) 偏差 ~0.5%，l_A 偏差 ~1.5（≈17σ in Planck units）。

3. **D_V formula inverted** — D_V = [c·z·D_M²/H(z)]^(1/3) 俾我寫錯做 [c·z·D_M²/D_H]^(1/3) 其中 D_H=c/H(z)——即係 H(z) 喺 numerator 而唔係 denominator。呢個令 D_V 細咗 ~3.6 倍，從兩個受影響嘅 BAO points 貢獻 χ²~1400。

呢啲 bug 嘅教訓係：

> **純靠 LLM 生成嘅 code 必須做 systematic audit。一個 sign error、一個 constant swap，可以殺死成個 project。**

### 5.3 真正的 breakthrough：DESI DR2 BAO 數據

之前我係用 SDSS DR12 BAO 數據（10點，BOSS + eBOSS）。但呢啲數據有個 convention trap——SDSS 嘅 D_M_over_rs 欄位其實係 D_M in Mpc（如 1512 Mpc at z=0.38），唔係 ratio。

**轉用 DESI DR2 BAO 之後，成個分析改變咗。**

DESI DR2 提供 12-vector consensus data（6個 redshift bins，D_M/r_d, D_H/r_d, D_V/r_d），附帶完整嘅 block-diagonal 12×12 協方差矩陣。呢個係 2024 年最強嘅 BAO 數據。

結果：DESI DR2 BAO + Planck CMB + Pantheon SNe + H(z) CC 四組數據聯合 fit，IDM χ²=1078.2 vs ΛCDM χ²=1129.0（Δχ²=-50.8, ΔAIC=-46.8）。

### 5.4 MCMC：從 grid 到置信區間

Initial grid scan 發現 ε=+0.25（z_c=0.70）。但呢個係 coarse-grid artifact——13 個 ε steps 太粗，resolve 唔到 ε-z_c 嘅 degeneracy ridge。

之後升級到 MCMC（emcee），24 walkers × 3500 steps，38 seconds run time。關鍵 speed trick：用 cumulative trapezoidal integral + interpolation 代替 1048 個 scipy quad integration。

**MCMC 結果完全改變咗 picture：**

| Method | ε | α = ε√44 | z_c |
|:-------|:--:|:--------:|:---:|
| Grid scan | +0.25 | 1.66 | 0.70 |
| **MCMC MAP** | **0.1545** | **1.025** | **0.24** |
| MCMC 68% CI | [0.109, 0.297] | [0.72, 1.97] | [0.45, 1.35] |

**MCMC MAP 嘅 ε=0.1545 偏離 CLT 預測（1/√44=0.151）僅 2.5%**。呢個係 non-trivial 嘅——model 嘅核心公式 ε=α/√N 由 CLT 推導，完全冇 tuning，居然 match 到數據最佳值。

Grid scan 嘅 ε=+0.25（α=1.66）係 coarse-grid artifact。

### 5.5 三碼 MCMC cross-validation

為咗確保結果唔係 code bug 嘅 artifact，我用三套獨立實作 cross-validate：

| Implementation | Engine | D_C(z) method | D_V bug? |
|:--------------|:------:|:--------------:|:--------:|
| mcmc_quick.py | emcee | Quadrature (verified) | ❌ No |
| mcmc_fast.py | emcee | 500-pt interpolation | ❌ No |
| mcmc_idm.py | emcee | 500-pt interpolation | ✅ Fixed |

最大 inter-code difference：ΔH₀<0.03, Δε<0.002, Δz_c<0.05——三套 code 嘅結果 statistically indistinguishable。

---

## 六、延伸想法：ε-z_c degeneracy

### 6.1 發現 degeneracy

MCMC posterior reveal 咗一個 characteristic degeneracy：ε 同 z_c 有雙峰結構。

| Family | ε | z_c | α | χ²_min |
|:-------|:----:|:---:|:---:|:-------:|
| **Low-z_c (CLT branch)** | 0.154 | 0.61 | 1.02 | 1074.6 |
| **High-z_c (deep dip)** | 0.209 | 1.26 | 1.39 | 1074.9 |

Δχ²<0.3——current data 完全無法分辨。

### 6.2 點解會有呢個 degeneracy？

所有距離數據（BAO/SNe/H(z)/CMB compressed）都係積分 ∫dz/H(z)。一個 dip 嘅 integrated area ε·z_c 係固定嘅（≈0.09），但位置可以 swap：

- CLT branch：shallow dip 喺早期（z_c≈0.6, ε≈0.15）
- Deep branch：deep dip 喺後期（z_c≈1.3, ε≈0.21）

兩者嘅 integrated dip area 一樣。

### 6.3 延伸思考：點樣打破 degeneracy？

我諗咗好幾個方向：

1. **fσ₈ growth rate** — RSD 數據來自結構形成嘅動力學。直覺上，dip 嘅位置應該影響 growth rate。但 full ODE 解出嚟，兩個分支嘅 f(z) 差異 <0.4%（數據誤差 ±4-12%）。**fσ₈ 無法分辨**——growth equation 都係積分。

2. **First-principles N_horizon** — 用純幾何方法計 N_horizon（完全唔用 ε）。三種計數方法：
   - 1D radial（D_C(z_rec)/D_H(z_c)）= 4.38
   - 2D surface = 19.2
   - 3D volume = 84.1
   - **√(N₂×N₃) ≈ 40.2** 同 CLT-derived N=41.9 差 4%

   呢個 4% 吻合只喺 CLT branch (z_c=0.6) 成立。Deep branch (z_c=1.3) 嘅 ratio 係 4.96——差好遠。

3. **CLASS CMB Boltzmann solver** — 直接 patch CLASS 源碼，加入 IDM dark energy。結果：CLT branch 達 χ²_CMB=0.7（近乎完美 fit），Deep branch χ²=5.2（disfavored）。

**三層 degeneracy breaking 完成：距離數據 → 幾何 N_horizon → CMB 功率譜。**

---

## 七、延伸想法：Secular Sync

### 7.1 核心問題

IDM 嘅瞬態 dip 函數 f(z)=1-ε·(z/z_c)·exp(-z/z_c) 有個問題：佢預測 f(0)=f(∞)=1——sync 完成後完全恢復。但 IDM 嘅物理框架話 sync 係 ongoing process，唔應該有「恢復」。

呢個矛盾牽涉到一個更大嘅問題：**宇宙嘅最終 fate 係乜？** Gemini 嘅 islandization 分析話全部暗能量會 remain（De Sitter 隔離），但 IDM 嘅 global sync 目標話 sync 會持續消耗暗能量。邊個啱？

### 7.2 推導 f_secular(z)

從同步動力學推導第二層修正：

$$N_{\rm horiz}(z) = 44 \times \frac{D_C(z)}{D_C(z_c=0.7)}$$

$$f_{\rm secular}(z) = 1 - \lambda \int_0^z \frac{\sqrt{N_{\rm 1D}(z')}}{1+z'} dz'$$

I_total = ∫₀^∞ √N/(1+z) dz ≈ 55.72（收斂到有限值）

### 7.3 結論：~80/20 分割

| λ | f_sec(∞) | DE consumed | feasibility |
|:-:|:---------:|:-----------:|:-----------:|
| 0 | 100% | 0% | baseline |
| 0.001 | 94% | 6% | preferred |
| 0.003 | 83% | 17% | upper bound |
| 0.018 | 0% | 100% | physical limit |

**~83-94% 嘅暗能量係 fundamental（量子真空）**，只有 ~6-17% 可以被 sync 消耗。

呢個 reconcile 咗 Gemini 同 IDM：
- **Gemini 啱**：大部分 Λ 係 fundamental → De Sitter islandization 基本正確
- **IDM 都啱**：sync 消耗一部分 Λ → global sync 有有限目標
- **結論**：兩個 model 講緊同一個 universe 嘅兩個層面

### 7.4 宇宙擴張嘅三個 Regime

| Regime | Redshift | 物理過程 | 觀測信號 |
|:-------|:---------|:---------|:---------|
| **Transient dip** | z~0.7 | 主要 sync event，CLT 1/√44 | Δχ²=-54.4 vs ΛCDM |
| **Secular consumption** | z<2 | 持續 sync，λ<0.003 | <1% correction（undetectable） |
| **Fundamental Λ** | all z | 量子真空 | ≥83% of observed DE |

---

**續 Part 3：最終結論與深層意義**