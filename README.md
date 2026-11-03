# IDM v1.0 — Information-Dynamics-Cosmology

**Information Dynamics Cosmology Model（資訊動力學宇宙論模型）— 1st Generation**

> `Information-Dynamics-Cosmology` 係版本代號，唔係理論名稱。理論名稱由 v1.0 到 v11.0 一直係 Information Dynamics Cosmology Model，從未改變。
>
> v1.0 的職責是**立框架**：把宇宙膨脹由「時空拉伸」重述為「資訊底層的同步壓力釋放」，並令原本自由擬合的凹陷深度 ε 第一次由中央極限定理閉合。
>
> 本 repo 同時係早期檔案庫：`Entropy-Gradient-Validation/` 與 `Lagrangian-Reverse-Engineering/` 兩個子專案以 embedded git 形式並存，屬本版產出的延伸工作。

## 緣起

起點唔係數據，而係三個對 ΛCDM 的不滿：

1. **Λ 係自由參數** — 真空能理論值與觀測值差約 120 個數量級，ΛCDM 只係把差異參數化，唔解釋
2. **H₀ tension** — SH0ES 73.0 ± 1.0 與 Planck 67.4 ± 0.5 相差 5σ，而 H₀ 與 Ω_m 相關，改一個就整壞另一個
3. **暗能量無機制** — Λ 係「就係咁」的常數，無 dynamics 可問「點解係呢個數」

第一個假設係：**膨脹係熵釋放機制**，唔係幾何屬性。類比來自平行運算系統——processing unit 之間的同步 latency 太高，系統就擴張以降低同步成本。由這個框架落出的第一條可測預言是凹陷函數：

```
f_ID(z) = 1 − ε · (z/z_c) · exp(−z/z_c)
```

z = 0 與 z → ∞ 皆回 1，只在 z ≈ z_c 有一個深度 ε/e 的凹陷。此凹陷就是同步事件的 signature——域開始同步、消耗一部分暗能量、同步完成後恢復。

## 版本演進：緣起 → v1.0

IDM 係**逆向工程**：由已觀測的常數與結構反推生成機制，唔預設公理、唔設「核心假設」。因此版本之間的分野唔係假設改變，而係**逆推對象改變**。

### 前版

v1.0 之前無論文、無 repo、無 register。當時唯一可對照的「前版」係 ΛCDM 本身：六個參數 fit 通 CMB 功率譜、BAO、超新星距離、大尺度結構。它的成功無可否認，問題在於成功的性質。

### 根因：參數化 ≠ 解釋

ΛCDM 每一項成功都係曲線貼合，而唔係機制。更麻煩的是，本版初稿犯同一個病：凹陷深度最初設 ε = 0.05，純屬估計，深度無來源、無約束、可任意調——與 Λ 同一個結構。

換言之，**如果 ε 繼續係自由參數，IDM 只係換咗名字的 ΛCDM**，只係多一個可調旋鈕。判別方法只有一個：ε 能否由框架自己算出，而唔係由數據 fit 出。

### 本版

v1.0 因此只做三件事：

1. **立三層** — OAS（本體）→ 資訊動力（機制）→ 投影宇宙論（可測），三層職責分離，唔互相借用結論
2. **閉合 ε** — 用中央極限定理把 ε 寫成 `α/√N_horizon`：因果域之間同步的統計雜訊，α 為協議效率（理想同步 α = 1），N_horizon 為 z_c 時的因果域數
3. **對觀測** — 六項獨立檢驗 + 四資料集聯合擬合 + 三套獨立 MCMC + CLASS Boltzmann 實作

驗收標準：凹陷深度須在無 tuning 下落於觀測值附近。MCMC MAP 給出 ε = 0.1545，對 CLT 預言 1/√44 = 0.151 差 2.5%，α = 1.025。這一步令 IDM 由「換名 ΛCDM」變成零自由參數的候選模型。

## 研究路線

六層，每層只做一件事，上層不能借用下層的結論：

**L1 本體層** — 基底係咩

- 基底唔係時空，而係可數的 N_total 個本體位址（OAS, Ontological Address Space）
- 此層無空間延展；「距離」唔係基本量

**L2 動力層** — 位址之間發生咩

- Landauer 宇宙論原理：每次資訊處理事件產生 k_B ln 2 熵，`dS/dt = η k_B İ_total`（η ~ O(1)）
- 同步成本函數：`İ_total ∝ ρ_addr² · V_horizon / τ_sync`，配合 `τ_sync ∝ 1/H` 得 `İ_total ∝ ρ_addr² H⁻²`
- 熵密度超臨界則系統以擴張釋放熵：`H² ∝ (1/V)(dS/dt)`

**L3 投影層** — 觀測到咩

- 3+1 維時空、Friedmann 方程、Hubble flow 皆為 OAS 動力學的平均場投影
- `H²(z) = H₀²[Ω_m(1+z)³ + Ω_r(1+z)⁴ + Ω_Λ f_ID(z)]`
- z ≪ 1 投影函數近似常數則回復 GR；z ~ z_c 出現凹陷；z ≫ z_c 投影高度扭曲，對應暴脹

**L4 參數閉合層** — 凹陷深度由邊度嚟

- `ε = α/√N_horizon(z_c)`，α 為同步協議效率，理想同步 α = 1
- N_horizon 三種獨立幾何計數：1D 徑向殼層 `N₁ = D_C(z_*)/D_H(z_c) = 4.38`、2D 天區 `N₂ = N₁² = 19.2`、3D 視界體 `N₃ = N₁³ = 84.1`
- 幾何均 `√(N₂N₃) = 40.2` 對 CLT 值 `N_CLT = α²/ε² = 41.9` 差 4%；比值過 1 的點落在 z_c ≈ 0.6

**L5 觀測層** — 同數據對質

- 六項獨立檢驗（H(z)、BAO 距離、fσ₈、CMB R、CMB l_A、ISW）
- 四資料集聯合擬合（DESI DR2 BAO + Planck CMB + Pantheon SNe + H(z)），grid scan 後接三套獨立 MCMC
- CLASS v3.3.4 源碼級修改，把 f_ID(z) 寫入 `background_w_fld()`，做首個 Boltzmann 層 IDM 檢驗

**L6 延伸層** — 框架講到幾遠

- C3 質量即熵副作用：`m_i = m₀ · N_cluster · (1 + α · dS_local/dt)`。本版自認仍屬 labelling，唔係預測理論
- C4 暗物質即弱投影的熵梯度：有重力、無電磁作用，不需新粒子。本版只有概念 sketch
- 兩項皆明確標「未驗證」，唔升格為結果

## 項目結構

```
Information-Dynamics-Cosmology/         # 本 repo（v1.0 主體）
├── README.md                          # 本檔
├── core/                              # 模型核心（8 個 .py）
│   ├── cosmic_idm.py                  # 原始 IDM 膨脹模型
│   ├── cosmic_model.py                # 共用宇宙學函數
│   ├── derive_eps_from_N.py           # 由 N_horizon 反推 ε
│   ├── n_horizon_first_principles.py  # 1D／2D／3D 三種幾何計數
│   ├── f_secular_v3.py                # 長期同步消耗推導（1D 視界 + √N 標度）
│   ├── joint_final.py                 # 四資料集聯合擬合（quadrature 積分）
│   ├── joint_lambda_scan_fast.py      # λ 掃描
│   └── observational_data.py          # 觀測數據載入
├── mcmc/                              # 三套獨立 MCMC（4 檔）
│   ├── mcmc_quick.py                  # quadrature 實作（gold standard）
│   ├── mcmc_fast.py                   # 插值加速版
│   ├── mcmc_idm.py                    # 獨立重寫（D_V 公式已審）
│   └── mcmc_results_summary.md        # corner plot 與參數摘要
├── validation/                        # 13 個驗證腳本
│   ├── idm_comprehensive_tests.py     # 六項檢驗批次
│   ├── joint_validation.py            # 聯合驗證 v1
│   ├── joint_validation_v2.py         # 聯合驗證 v2
│   ├── joint_validation_v3.py         # 聯合驗證 v3（定稿）
│   ├── cmb_test.py                    # CMB shift parameter 檢驗
│   ├── f_secular_derivation.py        # 同步消耗推導初版
│   ├── verification_check.py          # 關鍵公式數值一致性
│   └── grid_search.py                 # 參數網格掃描
├── class_idm/                         # 修改版 CLASS 源碼（3 檔）
│   ├── background.h                   # enum 加 IDM、struct 加 eps_idm／zc_idm
│   ├── background.c                   # IDM 有效 w(z) 寫入 background_w_fld()
│   └── input.c                        # 解析 eps_idm、zc_idm
├── data/
│   └── fs8_data_compilation.txt       # 增長率數據編纂
├── results/                           # 5 個輸出
│   ├── mcmc_final.png                 # 4D corner plot
│   ├── mcmc_corner.png                # 參數相關性
│   ├── mcmc_final.npz                 # MAP 鏈
│   ├── mcmc_chain.npz                 # 完整鏈
│   └── idm_w_of_z.png                 # w(z) 曲線
├── docs/                              # 12 份文件（三代論文並存）
│   ├── IDM_paper_v1.0.md              # 本版主論文（英）
│   ├── IDM_paper_v1.0_zh.md           # 本版主論文（中）
│   ├── IDM_report_v3.1.md             # 完整報告 v3.1（英）
│   ├── IDM_report_v3.1_zh.md          # 完整報告 v3.1（中）
│   ├── IDM_paper_v4.0.md              # 聯合驗證論文 v4.0（英）
│   ├── IDM_paper_v4.0_zh.md           # 聯合驗證論文 v4.0（中）
│   ├── idm_final_conclusion.md        # 六項檢驗後的最終結論
│   ├── IDM_narrative_complete.md      # 敘事全文（緣起 → 過程 → 結論）
│   ├── IDM_narrative_part1.md         # 敘事一：緣起與動機
│   ├── IDM_narrative_part2.md         # 敘事二：探索過程
│   ├── IDM_narrative_part3.md         # 敘事三：結論
│   └── alerts.md                      # 六項批次運算的工作單
```

> 重構（2026-11-03）：原本內嵌的兩個子專案已分拆為平級獨立 repo —
> `../Entropy-Gradient-Validation/`（C3/C4 實測）與 `../Lagrangian-Reverse-Engineering/`（Horndeski 作用量逆推）。
> 三者各自帶 `.git` 與獨立 remote；本 repo 只保留 v1.0 主體（48 個 tracked 檔）。

## 狀態快照

```
頂層目錄：7（core／mcmc／validation／class_idm／data／docs／results）
文件總數：48（本 repo 主體）
markdown：12 份；Python：24 個（core／mcmc／validation／class_idm patch）

本體與動力層：三層框架完成，C3／C4 停在 labelling
參數閉合  ：ε 已閉合；z_c 未閉合（仍由數據定）
觀測層    ：六項檢驗完成；三碼 MCMC 一致；CLASS 已打通
子專案    ：已分拆為平級 repo（Entropy 4 測試完成；Lagrangian 4 階段完成）
未完成    ：α 的第一原理推導、z_c 的理論來源、音界自洽重算、C3/C4 定量化
```

> 本版屬 IDM 時期，未進入 O/P/M register 體系；`docs/` 內 v3.1 報告與 v4.0 論文係後續世代的產物，在此作為檔案保存。

## 已閉合關鍵結果（選錄）

- **2026-06-24（`docs/IDM_paper_v1.0.md`）**：三層框架成立；H(z) 39 點 χ² = 24.30 對 ΛCDM 30.14，Δχ² = +5.84，且少一個自由參數
- **2026-06-24（`docs/idm_final_conclusion.md`）**：六項獨立檢驗完成。H(z) Δχ² +5.84、CMB R pull −0.9σ（ΛCDM −1.9σ）、CMB l_A pull −23σ（ΛCDM −31.7σ）、ISW 差 −0.7% 皆 IDM 較好；BAO −1.3、fσ₈ −0.2 為平手或 ΛCDM 微勝
- **2026-06-25（`mcmc/mcmc_results_summary.md`）**：四資料集聯合 MCMC，H₀ = 69.17⁺⁰·⁵⁶₋₀·₅₇、Ω_m = 0.2973、ε = 0.1965、χ² = 1074.9、ΔAIC = +50.1；ε > 0 在 >2σ 確認
- **2026-06-25（`docs/IDM_paper_v4.0.md`）**：CLT 閉合。MAP ε = 0.1545 給 α = 1.025，對理想同步 α = 1 差 2.5%；ε_CLT = 1/√44 = 0.151 落在中位數 0.03σ 內
- **2026-06-25**：CMB 張力解除。χ²_CMB 由 37.9（ΛCDM，l_A 偏 5.5σ）跌至 1.3（IDM，偏 1.2σ）
- **2026-06-25**：三碼交叉驗證（mcmc_quick／mcmc_fast／mcmc_idm）中位數最大差 ΔH₀ < 0.03、Δε < 0.002、Δz_c < 0.05
- **2026-06-26（`docs/IDM_paper_v4.0.md` 4.8 節）**：CLASS v3.3.4 源碼級 IDM 打通。CLT 支（ε = 0.14、z_c = 0.6、H₀ = 69）χ²_CMB = 0.7，Deep 支（z_c = 1.3）5.2；合併 ΔAIC 為 −47 對 −34，ε–z_c 簡併被 CMB 打破
- **2026-06-26**：N_horizon 幾何計數與 CLT 自洽。√(N₂N₃) = 40.2 對 N_CLT = 41.9 差 4%，比值過 1 於 z_c ≈ 0.6；Deep 支比值 4.96，差約 5 倍（負結果，反證 CLT 支）
- **2026-06-28（`Lagrangian-Reverse-Engineering/docs/stage4_result.md`）**：Horndeski 逆推四階段完成。典範 quintessence 被 phantom crossing（z < 0.59）否證；α_B ≈ +0.15、α_T = 0（GW170817）、指數屏蔽 ξ(φ) = ξ₀e^(−λφ)
- **2026-06-27（`Entropy-Gradient-Validation/docs/FINAL_REPORT.md`）**：fσ₈ 15 點殘差 z-trend +2.81σ/z（p < 0.001），但 σ₈(0) 擬合值 0.600 低於 Planck 0.81，報告自標為 z-trend 偏差，未宣稱信號

## 驗證紀律

IDM 的嚴謹性分六層：獨立性 → 交叉實作 → 公式審計 → Boltzmann 層 → 可否證性 → 記錄紀律。任一層不過，結果只停留「數值候選」。

### 一、檢驗獨立性

| 檢驗 | 數據源 | 獨立性要求 |
|---|---|---|
| H(z) | 39 點 cosmic chronometers（Moresco+2016） | 與 BAO 為不同觀測量 |
| BAO 距離 | BOSS DR12 + eBOSS 10 點 | 用 D_V/r_s，不重用 H(z) |
| fσ₈ | RSD 13–15 點 | 增長率，非背景膨脹 |
| CMB R | Planck 2020 | 距離先驗，非功率譜 |
| CMB l_A | Planck 2020 | 聲學標度，與 R 不同投影 |
| ISW | Planck l < 40 | 低ℓ段，獨立於高ℓ |

> 六項須逐項報 Δχ²，禁只報總和。單一數據集大幅改善可能係吸收其他數據集的張力。

### 二、交叉實作

- 三套 MCMC 獨立編寫（`mcmc_quick`、`mcmc_fast`、`mcmc_idm`），共用同一 likelihood
- 判準：中位數最大差 ΔH₀ < 0.03、Δε < 0.002、Δz_c < 0.05
- grid scan 與 MCMC 須並存；grid 值（ε = +0.25）與 MCMC 值（ε = 0.1965）分歧須明寫為 coarse-grid artifact，唔可以揀好聽的報

### 三、公式審計

- `D_V = [cz(1+z)² D_A² / H]^(1/3)` 逐項核，`validation/idm_comprehensive_tests.py` 與 `mcmc_idm.py` 對照
- 關鍵恆等式做數值檢查（`verification_check.py`）：f_ID(0) = 1、f_ID(z → ∞) → 1、BAO D_V/r_s pull
- 純 numpy／scipy，無外部 fitting 套件介入

### 四、Boltzmann 層

- CLASS 修改只准動三處：`enum equation_of_state`、`background_w_fld()`、參數解析
- 三組對照必跑：ΛCDM（ε = 0）、CLT 支、Deep 支
- 報 CMB 差異須同時報 Planck S/N 倍數。低於 0.01× 的差異唔可以稱「可測」

### 五、可否證性

- 凹陷必須係 dip 唔係 bump：ε > 0，且 z_c ~ 0.6
- 高 z 無偏離：f_ID(z) → 1，CMB 主各向異性不受影響
- 明確 kill 條件：若 DESI／Euclid 在 z ~ 0.6 測不到 H(z) 凹陷，或 ε 轉號，模型即否證
- H₀ 須為中間值 69.17；若未來獨立測量確認 73.0 或 67.4，模型 disfavored

### 六、記錄紀律

- 每個結論須可追至具體檔案（論文節號／腳本名／輸出檔）
- 未完成項明寫「PENDING」或「labeling」，禁寫成結果
- 子專案自帶 `.git`，其 commit 歷史與本 repo 分離，引用時須註明來源 repo
- 版本差異係逆推對象改變，唔係假設改變；跨版本引用前先核歸一化

## O/P/M 系列全列表

> 本版屬 IDM 時期，未進入 O/P/M register 體系。下表以 **D 系列**（Deliverable，交付項）編號，每項對應本 repo 內實際存在的文件或模組，唔對應 register pid。

| 編號 | 標題 | 核心發現 | 狀態 |
|------|------|----------|------|
| D1 | 三層本體框架：OAS → 資訊動力 → 投影宇宙論 | 空間為 sync latency 投影、時間為位址更新 cycle、膨脹為熵釋放閥；三層職責分離 | 全部閉合 |
| D2 | 同步凹陷函數 f_ID(z) | f = 1 − ε(z/z_c)e^(−z/z_c)；z→0 與 z→∞ 皆回 1，僅 z ≈ z_c 有深度 ε/e 的凹陷 | 全部閉合 |
| D3 | ε 的 CLT 閉合 | ε = α/√N_horizon；MAP ε = 0.1545、α = 1.025，對 1/√44 = 0.151 差 2.5% | 全部閉合 |
| D4 | 六項獨立檢驗 | H(z) Δχ² +5.84、CMB R pull −0.9σ 勝；BAO −1.3、fσ₈ −0.2 平手或微敗 | 部分閉合 |
| D5 | 三碼 MCMC 交叉驗證 | 三套實作中位數差 ΔH₀<0.03、Δε<0.002、Δz_c<0.05；χ² = 1074.9、ΔAIC = +50.1 | 全部閉合 |
| D6 | CLASS Boltzmann 實作 | 改 background.h／background.c／input.c；CLT 支 χ²_CMB = 0.7、Deep 支 5.2，簡併被打破 | 全部閉合 |
| D7 | N_horizon 三種幾何計數 | N₁ = 4.38、N₂ = 19.2、N₃ = 84.1；幾何均 40.2 對 CLT 41.9 差 4% | 部分閉合 |
| D8 | C3：質量即熵副作用 | 形式 m = m₀·N_cluster·(1+α·dS/dt) 仍屬 tautology，cluster potential 未推導 | 開放 |
| D9 | C4：暗物質即弱投影熵梯度 | 只有概念 sketch；有重力無電磁的論證未定量，power spectrum 未計算 | 開放 |
| D10 | Lagrangian 逆推（Horndeski 四階段） | 典範 quintessence 被 phantom crossing 否證；α_B ≈ +0.15、α_T = 0、指數屏蔽 | 部分閉合 |
| D11 | fσ₈ 殘差 z-trend（子專案） | 15 點擬合 χ² = 4.5–5.1；σ₈(0) = 0.600 低於 Planck 0.81，trend 疑為偏差 | 部分閉合 |
| D12 | z_c 的理論來源 | z_c ≈ 0.6 仍由數據定；未由複合期因果視界與物質輻射比導出 | 開放 |