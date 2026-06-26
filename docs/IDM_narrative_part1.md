# 資訊動力學宇宙論（IDM）：緣起與動機

## — 我點解會做呢件事，同想達成啲乜

---

呢個 project 嘅起源唔係喺天文台或者學術會議，而係喺我日常嘅思考入面——嗰種「標準模型唔完整」嘅直覺，一直啃住。

---

## 一、點解會開始

### 1.1 ΛCDM 嘅痛點

標準宇宙學模型——ΛCDM——非常成功，呢個冇人可以否認。佢用六個參數 fit 晒 CMB 功率譜、BAO 重子聲學振盪、超新星距離、大尺度結構……但佢嘅成功背後，有三個嘢我一直覺得唔對路：

1. **Λ（宇宙學常數）係個 free parameter** — 你冇辦法推導佢出嚟，就係咁 set 一個數值。理論預測嘅真空能 density 同觀測值差 120 個 orders of magnitude，但 ΛCDM 唔解釋呢個 discrepancy，佢只係參數化咗佢。

2. **H₀ tension** — SH0ES 話 H₀=73.0，Planck 話 H₀=67.4，相差 5σ。呢個唔係一個小誤差，呢個係兩個唔同嘅宇宙。ΛCDM 理論上無法 reconcile 呢個分歧，因為 H₀ 同 Ω_m 嘅影響係 correlated 嘅——你改一個就整壞另一個。

3. **暗能量冇物理解釋** — Λ 係一個「就係咁」嘅 constant。冇 mechanism，冇 dynamics，冇乜嘢可以問「點解係呢個數」。

我嘅問題好簡單：**如果暗能量唔係 fundamental constant，而係某種 emergent 現象呢？**

### 1.2 資訊動力學嘅第一原理起點

呢個諗法唔係由數據開始，而係由一個更根本嘅問題開始：**宇宙膨脹嘅底層機制係乜？**

標準 cosmology 嘅答案是：時空本身喺度膨脹。但呢個答案對我嚟講唔夠——佢冇解釋點解時空會膨脹，只係描述咗膨脹嘅 rate。

我嘅第一個假設係：

> **膨脹可能係一種熵釋放機制。**

想像一個系統——唔係時空，而係某種底層嘅「位址空間」（OAS, Ontological Address Space）——入面有大量嘅 processing unit。當佢哋之間嘅同步壓力（sync pressure）超過某個 threshold，系統就要通過膨脹嚟 release entropy。

呢個 analogy 嚟自工程系統：當你有一個 parallel processing system，processing units 之間需要同步。如果同步嘅 latency 太高，系統要擴張（加多啲 units）嚟降低同步成本。

### 1.3 第一個公式嘅誕生

從呢個 conceptual framework，第一個 testable prediction 就係：

$$f_{\text{ID}}(z) = 1 - \varepsilon\cdot\frac{z}{z_c}e^{-z/z_c}$$

呢個函數嘅形狀唔係 random 嘅：
- 喺 z=0（而家），f=1（冇 deviation）
- 喺 z→∞（早期宇宙），f→1（同樣冇 deviation）
- 喺 z≈z_c，有一個 dip，深度 ε/e

呢個 dip 就係同步事件嘅 signature——因果域之間開始 sync，消耗咗一部份 dark energy，令有效嘅 Ω_Λ 暫時下跌，然後 sync 完成後恢復。

初期我設定 ε=0.05，完全係估嘅。完全冇諗過後嚟會發現呢個數字其實可以由中央極限定理（CLT）推導出嚟。

---

## 二、想達成啲乜目的

### 2.1 最小目標：證明暗能量唔係 fundamental

我最初嘅目標係好實在嘅：

> **用觀測數據證明暗能量可以隨 redshift 變化，從而推翻「暗能量 = 宇宙學常數」呢個假設。**

ΛCDM 嘅 Λ 係 constant——佢唔隨時間變。如果你可以 show 到 Ω_Λ(z) 偏離 constant，ΛCDM 就錯。就係咁簡單。

呢個唔係一個新 idea——w₀wₐ parametrization（CPL）本身就係 variation。但 IDM 唔同嘅地方係：我嘅 dip 函數 f(z) 唔係 ad-hoc 嘅 polynomial，而係從同步動力學推導出嚟嘅形狀。

### 2.2 中期目標：搵出一個「第一原理推導」嘅 model

ΛCDM 嘅六個參數全部係 phenomenological——你從數據 fit 出嚟，唔係從第一原理計出嚟。IDM 嘅中期目標係：

> **搵到一個機制，令到至少一個「free parameter」可以從理論推導，而唔係從數據 fit。**

ε 就係呢個候選人。如果 ε 可以由 CLT 決定（ε ∝ 1/√N），咁 IDM 就比 ΛCDM 少一個 free parameter——呢個係統計上嘅 decisive advantage。

### 2.3 終極目標：重新定義宇宙嘅本質

如果 IDM 係啱嘅，咁宇宙唔係一個「三維空間隨時間膨脹」嘅幾何物體，而係：

> **一個自我組織嘅資訊處理系統，其膨脹歷史係底層同步壓力釋放嘅宏觀表現。**

呢個 paradigm shift 嘅幅度，類比於牛頓引力 → 廣義相對論。唔係修修補補，而係換咗個 ontology。

呢個聽落好大膽——事實上係非常進取。但我嘅 philosophy 係：**先有 framework，然後假裝佢係假嘅，然後盡全力用數據去 falsify 佢。** 如果佢 survive 到多輪測試，咁就算 framework 本身進取，結論都係堅實嘅。

---

## 三、最初嘅想法同直覺

### 3.1 因果域同步（Causal Domain Sync）

核心 idea 係：宇宙喺任何 epoch 都可以分為 ~N_horizon 個因果域（causally disconnected patches），每個域嘅半徑 ~c/H(z)。呢啲域喺早期宇宙係完全 causal disconnected——佢哋之間無法交換資訊。

但隨住宇宙膨脹，粒子視界（particle horizon）增長，之前 disconnected 嘅域開始 overlap。呢個 overlap 令到域之間可以 sync——交換關於 collective vacuum state 嘅資訊。

呢個 sync process 有代價：sync 消耗 entropy，而 entropy 嘅 source 就係 dark energy。所以 sync 導致 Ω_Λ 暫時下跌——就係個 dip。

### 3.2 N_horizon ≈ 44 嘅第一個直覺

我 initial 嘅想法係：N_horizon 應該係沿一條視線（line-of-sight）上因果域嘅數量。喺 z_c ≈ 0.6，粒子視界 D_C(z_c) ≈ 2580 Mpc，Hubble 距離 D_H(z_c) ≈ c/H(z_c) ≈ 3300 Mpc。

粗略計算：沿半個宇宙（180° 大圓）可容納 ~22 個域，成個大圓 ~44 個。

呢個「44」後來被 MCMC 驗證——從 ε=0.1545 反推 N_horizon = α²/ε² ≈ 42，同 44 差 4%。呢個吻合唔係 tune 出嚟嘅——係兩個獨立推導 converge 到同一個數。

### 3.3 Dip 嘅物理形狀

f(z) = 1 - ε·(z/z_c)·exp(-z/z_c) 嘅 shape 唔係 ad-hoc：

- **z → 0**：f→1（而家嘅暗能量同 ΛCDM 一樣）
- **z ≪ z_c**：f ≈ 1 - ε·z/z_c（線性增長，sync rate 逐漸增加）
- **z = z_c**：最大 dip，深度 ε/e ≈ 0.055（~5.5% DE suppression）
- **z ≫ z_c**：f → 1（早期宇宙未開始 sync，同 ΛCDM 一樣）

呢個「中間弱、兩邊強」嘅 shape 係 sync process 嘅 signature——sync 需要時間 accumulate，完成後效果衰減。

---

## 四、呢個 project 嘅獨特之處

我唔係 professional cosmologist。呢個 project 冇 grant，冇 supervisor，冇 peer review process。但呢個反而係佢嘅 strength：

1. **冇 disciplinary bias** — 我唔會因為「呢個做法唔 standard」就唔做
2. **冇 publish-or-perish 壓力** — 我可以俾 CS 嘅 parallel processing analogy 啟發，而唔係局限喺 cosmology 嘅 standard toolbox
3. **第一原理優先** — 同 standard cosmology 嘅「先有理論框架，然後 fit 參數」相反，我係「先有機制，然後推導參數，然後用數據 check」

呢個 approach 嘅風險係：你可能完全錯。但回報係：如果你啱，你發現嘅 mechanism 唔係 free parameter fit 出嚟嘅 artifact。

---

**續 Part 2：探索過程與延伸想法**
