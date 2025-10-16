# aCSS-B

Code for replicating experiments in the paper *"Conditioning on posterior samples for flexible frequentist goodness-of-fit testing."*

## Simulation Workflow Instructions

This document provides instructions for running simulations in a SLURM-enabled environment.  
SLURM is required to execute the following commands since they use job scheduling and parallelization features.  
For users without SLURM, the corresponding `Run.R` files can be executed in a simple loop with appropriate seed values.

---

### Simulation 1: Logistic Regression

This experiment plots the power curves for **[aCSS]** and **aCSS-B** methods along with the oracle.

1. Navigate to the directory `logistic`:
   ```bash
   cd logistic
   ```

2. Run the following R scripts:
   ```bash
   Rscript run_aCSS_logistic.R
   Rscript run_BaCSS_logistic.R
   ```

3. Plot the power curve:
   ```bash
   Rscript plot.R
   ```

---

### Simulation 2: Mixture of Gaussians

This experiment plots the power curves for **[regularized aCSS]** and **aCSS-B** methods along with the oracle.  
The *p*-values used to generate the power curves for the oracle and **[regularized aCSS]** methods were provided by the authors of the **[regularized aCSS]** paper.  
The seeds used to run `Run.R` were chosen to match those used by the original authors, enabling a fair comparison of performance across methods.

1. Navigate to the directory `gaussian_mixture`:
   ```bash
   cd gaussian_mixture
   ```

2. Submit jobs using SLURM:
   ```bash
   sbatch --array=1002-2000:2 run.sh
   ```

3. Compile results:
   ```bash
   Rscript plot.R
   ```

---

### Simulations 3–5: Group Sparsity, Rank-One Matrix, and Linear Spline

The procedure for generating the power curves for these examples is the same.  
For brevity, we list the common steps below.

1. Navigate to the appropriate directory (`rank_one_matrix`, `group_sparsity`, or `linear_spline`):
   ```bash
   cd rank_one_matrix
   # or: cd group_sparsity
   # or: cd linear_spline
   ```

2. Submit jobs using SLURM:
   ```bash
   sbatch --array=1-500 run.sh
   ```

3. Compile results:
   ```bash
   Rscript plot.R
   ```

---

## References

Barber, R. F., & Janson, L. (2022). *Testing goodness-of-fit and conditional independence with approximate co-sufficient sampling.*  
*The Annals of Statistics, 50*(5), 2514–2544.  
[https://doi.org/10.1214/22-AOS2208](https://doi.org/10.1214/22-AOS2208)

Zhu, W., & Barber, R. F. (2023). *Approximate co-sufficient sampling with regularization.*  
*arXiv preprint* [arXiv:2309.08063](https://arxiv.org/abs/2309.08063)

[aCSS]: https://doi.org/10.1214/22-AOS2208
[regularized aCSS]: https://arxiv.org/abs/2309.08063
