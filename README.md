# Information Dynamics Cosmology (IDM)

**A first-principles cosmological model where dark energy emerges from information synchronization between causal domains.**

| | Status |
|:---|:-------|
| Data | DESI DR2 BAO + Planck CMB + Pantheon SNe + H(z) CC + fσ₈ |
| Validation | Grid scan + 3-code MCMC cross-validation + CLASS Boltzmann solver |
| Result | ΔAIC = +50.1 (IDM vs ΛCDM) |
| Core prediction | ε = α/√N_horizon — CLT closure verified to 2.5% |

---

## Repository Structure

```
Information-Dynamics-Cosmology/
├── core/                    # Core model code (Python)
│   ├── cosmic_idm.py        # Original IDM model (v1)
│   ├── joint_final.py       # Audited joint validation pipeline
│   ├── n_horizon_first_principles.py  # First-principles N_horizon derivation
│   ├── f_secular_v3.py      # Secular sync consumption derivation
│   └── ...
├── mcmc/                    # MCMC implementations
│   ├── mcmc_quick.py        # Quadrature (gold standard)
│   ├── mcmc_fast.py         # Fast interpolation
│   ├── mcmc_idm.py          # Independent implementation
│   └── mcmc_results_summary.md
├── validation/              # Validation scripts
│   ├── idm_comprehensive_tests.py  # 6-test batch
│   └── joint_validation_v3.py
├── class_idm/               # Modified CLASS v3.3.4 C source
│   ├── background.h         # Added IDM equation of state
│   ├── background.c         # IDM w(z) implementation
│   └── input.c              # IDM parameter parsing
├── docs/                    # Papers and reports
│   ├── IDM_paper_joint_validation.md  # v4.0 — Main paper (648 lines)
│   ├── IDM_report_complete.md         # v3.1 — Full report
│   ├── idm_final_conclusion.md        # Final conclusions
│   └── IDM_narrative_complete.md      # Narrative: origin → exploration → conclusions
├── data/
│   └── fs8_data_compilation.txt       # Growth rate data
└── results/                 # Output figures
    ├── mcmc_final.png        # 4D corner plot
    └── mcmc_corner.png
```

## Core Concept

The Friedmann equation is modified by a sync-induced modulation of dark energy:

$$H^2(z) = H_0^2\left[\Omega_m(1+z)^3 + \Omega_r(1+z)^4 + \Omega_\Lambda \cdot f_{\rm ID}(z)\right]$$

$$f_{\rm ID}(z) = 1 - \varepsilon\cdot\frac{z}{z_c}\exp\!\left(-\frac{z}{z_c}\right)$$

The dip amplitude ε is predicted by the Central Limit Theorem:

$$\varepsilon = \frac{\alpha}{\sqrt{N_{\rm horizon}(z_c)}}$$

where N_horizon ≈ 44 causal domains, and α ≈ 1 (ideal sync efficiency). The MCMC MAP gives ε = 0.1545 — within **2.5% of the CLT prediction** (1/√44 ≈ 0.151).

## Key Results

| Parameter | ΛCDM | IDM (MCMC 68% CI) |
|:----------|:----:|:-----------------:|
| H₀ [km/s/Mpc] | 67.0 | 69.17⁺⁰·⁵⁶₋₀·₅₇ |
| Ω_m | 0.318 | 0.2973⁺⁰·⁰⁰⁵⁵₋₀·₀₀₅₁ |
| ε | 0 (fixed) | 0.1965⁺⁰·¹⁰⁰¹₋₀·⁰⁸⁸⁰ |
| MAP ε | — | **0.1545** (α=1.025, 2.5% from CLT) |
| χ²_total | 1129.0 | 1074.9 |
| ΔAIC | 0 | **+50.1** |

All four datasets independently prefer IDM over ΛCDM.

## Dependencies

```bash
pip install numpy scipy matplotlib emcee corner
# For CLASS Boltzmann integration:
pip install classy
```

## Citation

```bibtex
@misc{Ng2026IDM,
  author = {Ng, Lucifer},
  title = {Information Dynamics Cosmology: A First-Principles Sync-Induced Dark Energy Model},
  year = {2026},
  howpublished = {\url{https://github.com/LuciferNg/Information-Dynamics-Cosmology}}
}
```

## License

MIT
