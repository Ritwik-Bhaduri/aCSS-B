# aCSS-B
Code for replicating experiments in the paper -- Conditioning on posterior samples for flexible frequentist goodness-of-fit testing

## Simulation Workflow Instructions
This document provides instructions for running simulations in a SLURM-enabled environment. SLURM is required to execute the following commands as they use job scheduling and parallelization features. For users without SLURM, one can simply run the corresponding Run.R files in a loop with appropriate seeds. 

### Simulation 1: Logistic Regression
This code will plot the power curves for the aCSS and aCSS-B methods along with the oracle.
1. Navigate to the directory `logistic`:
   ```bash
   cd logistic
   ```
2. Run the following R scripts:
   ```bash
   Rscript run_aCSS_logistic.R
   Rscript run_BaCSS_logistic.R
   ```
3. Plot power curve:
   ```bash
   Rscript plot.R
   ```

### Simulation 2: Mixture of Gaussians
This code will plot the power curves for the regularized aCSS and aCSS-B methods along with the oracle.The p values used to plot the power curves for the oracle and regularized aCSS methods were obtained from the authors of the regularized aCSS paper. The seeds used to run the run.R file were used to match the seed used by the authors of the regularized aCSS paper so that the performance of the methods could be compared.

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
   
### Simulation 3-5: Group sparsity, Rank one matrix, Linear Spline
The steps to be executed to generate the power curves for all these examples are the same, so for brevity we write those down together.

1. Navigate to the directory `rank_one_matrix` or 'group_sparsity' or 'linear_spline':
   ```bash
   cd rank_one_matrix
   # change to cd group_sparsity, cd linear_spline for the other examples
   ```
2. Submit jobs using SLURM:
   ```bash
   sbatch --array=1-500 run.sh
   ```
3. Compile results:
   ```bash
   Rscript plot.R
   ```


