# Domain Adaptation under Hidden Confounding — Codebase

Reference implementation of the **Generative Invariance (GI)** estimator from  
*"Domain adaptation under hidden confounding"* by Carlos García Meixide and David Ríos Insua.

![Simpson's paradox illustration](simpson.pdf)

---

## Repository contents

This repository contains the GI estimator and the scripts to reproduce the main simulation studies and the real-data experiment discussed in the paper.

### Files

- **`main.R`**  
  Implements the Generative Invariance (GI) estimator and the full workflow (load data → fit GI → evaluate).

- **`aux_gi.R`**  
  Helper functions for efficient environment-combination selection aimed at minimizing asymptotic variance.

- **`generate_6_2.R`**  
  Data-generating process used in the simulation setup of Section 6.2.

- **`sect6_4.R`**  
  Replicates Model 25 from Section 6.2 of *Causal Dantzig* (Rothenhäusler et al.) with an additional test environment.  
  Also reproduces the comparisons in Section 6.4 of the *Generative Invariance* paper.  
  Compares GI, Causal Dantzig, and Instrumental Variables (IV) in terms of predictive accuracy in an unseen domain.

- **`sect6_6.R`**  
  Reproduces the analysis of the SPRINT trial in Section 6.6.  
  Requires access to the dataset available upon request at: <https://biolincc.nhlbi.nih.gov/login/?next=/requests/type/sprint/>

---

## Quick start

1. **Clone the repository**
   ```bash
   git clone https://github.com/<your-org>/<your-repo>.git
   cd <your-repo>
   ```

2. **Set up the R environment**
   - Use R ≥ 4.2
   - If an `renv.lock` file is present:
     ```r
     install.packages("renv")
     renv::restore()
     ```
   - Otherwise, install packages as needed (see each script's `library()` calls).

3. **Run the estimator**
   ```r
   source("aux_gi.R")
   source("main.R")
   ```

4. **Reproduce Section 6.2 simulations**
   ```r
   source("generate_6_2.R")
   ```

5. **Reproduce Section 6.4 comparisons (GI vs. Causal Dantzig vs. IV)**
   ```r
   source("sect6_4.R")
   ```

6. **Reproduce Section 6.6 (SPRINT trial)**
   - Request access to the dataset from BioLINCC.
   - Adjust the file paths in `sect6_6.R` to point to the downloaded data.
   ```r
   source("sect6_6.R")
   ```

---

## Method summary

**Generative Invariance (GI)** leverages cross-environment structure to enable domain adaptation under hidden confounding. The estimator:

- uses invariance constraints derived from environment shifts,
- selects environment combinations to minimize asymptotic variance (`aux_gi.R`),
- produces predictors that transfer effectively to unseen domains.

For full details and theoretical results, see the paper.

---

## Reproducibility notes

- Random seeds are set within scripts where relevant.
- Scripts produce results (tables/plots) comparable to those in the paper.
- The SPRINT analysis depends on restricted-access data and requires the dataset from BioLINCC.

---

## Dependencies

If not using `renv`, you will typically need:

```r
install.packages(c(
  "data.table", "Matrix", "matrixStats", "ggplot2",
  "glmnet", "MASS", "dplyr", "tidyr", "purrr"
))
```

---

## Citation

If you use this code, please cite:

```bibtex
@article{meixide_riosinsua_hiddenconfounding_2025,
  title   = {Domain adaptation under hidden confounding},
  author  = {Garc{\'\i}a Meixide, Carlos and R{\'\i}os Insua, David},
  journal = {...},
  year    = {2025},
  note    = {Code: https://github.com/<your-org>/<your-repo>}
}
```

(Update venue details when available.)

---

## Contact

- **Carlos García Meixide** — your.email@domain
- **David Ríos Insua** — contact@domain

---

## License

Please specify your preferred license in `LICENSE` (e.g., MIT, Apache-2.0).

---

## Acknowledgments

We acknowledge the authors of comparative baselines (Causal Dantzig, IV) and the custodians of the SPRINT dataset (BioLINCC/NHLBI).
