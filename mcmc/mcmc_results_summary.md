## MCMC Corner Plot — 4-Parameter IDM Fit (DESI DR2 + CMB + SNe + H(z))  
**H₀ vs Ω_m vs ε vs z_c** — 24 walkers × 3500 steps (38s runtime)

![Corner plot](MEDIA:/home/wsl/cosmic_validation/mcmc_final.png)

| Parameter | MAP | Median ± 68% CI | Comments |
|:----------|:-----:|:----------------:|:---------|
| **H₀** | 69.1 | 69.20 ± 0.55 | Tightly constrained |
| **Ω_m** | 0.297 | 0.297 ± 0.005 | Tightly constrained |
| **ε** | 0.18 | 0.20 +0.10/−0.09 | **>2σ from 0** ✓ |
| **z_c** | 0.24 | 0.92 +0.41/−0.53 | Poorly constrained |

**χ²_min = 1074.6** | **Δχ² vs ΛCDM = +54.4** | **ΔAIC = +50.4**

### Degeneracies
- **H₀-Ω_m**: r = −0.989 (strong, expected — distance integrals)
- **H₀-ε**: r = +0.637 (moderate — more dip needs faster expansion to compensate)
- **ε-z_c**: r = +0.321 (weak — partially independent but coupled via dip shape)

### Key Findings
1. ε > 0 confirmed at >2σ → CLT sign correct ✓
2. z_c not well-constrained (degeneracy with ε absorbs uncertainty)
3. Need more BAO data at z=0.3–1.0 to break the degeneracy  
4. ΔAIC = 50.4 → IDM strongly preferred over ΛCDM
5. Full publication-quality contours now available

Generated: `mcmc_final.png` + `mcmc_final.npz` in `~/cosmic_validation/`