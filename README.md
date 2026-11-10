# idopNetwork-pipeline

R code for reconstructing trait/niche interaction networks from tri-variate
functional data. The pipeline takes three series of habitat–niche curves,
fits power-law trajectories, clusters the trajectories by shape, and then
infers a directed interaction network between the resulting modules.

## Pipeline overview

The three scripts are meant to be run **in order**; each stage writes an
intermediate `.RData` file that the next stage reads.

| Step | Script | Input | Main output |
| ---- | ------ | ----- | ----------- |
| 1 | `power_equation_fitting.R` | raw CSV of measured curves | `res1` (fitted power-law trajectories) |
| 2 | `tri_functional_clustering.R` | `res1` from step 1 | `cluster_result` (module assignment + BIC) |
| 3 | `network_reconstruction.R` | `cluster_result` from step 2 | module networks, edge/node tables, figures |

## Requirements

- R >= 4.0
- CRAN packages:

```r
install.packages(c(
  "idopNetwork", "glmnet", "orthopolynom", "mvtnorm", "pbapply",
  "ggplot2", "patchwork", "cowplot", "ggalt", "reshape2",
  "readxl", "readr", "dplyr", "Cairo"
))
```

`idopNetwork` (>= 0.1.2) provides the functional-clustering routines
(`bifun_clu_parallel`, `bifun_clu_plot`, `fun_clu_BIC`).
See <https://cran.r-project.org/package=idopNetwork>.

## Input data format

Step 1 reads a single CSV file:

- **column 1** — sample/module identifier (becomes the row names);
- **columns 2 – 94** — three consecutive blocks of numeric observations,
  corresponding to the three sample groups (20, 48 and 25 time points
  respectively). Adjust the `df1` / `df2` / `df3` column ranges at the
  bottom of `power_equation_fitting.R` if your data are partitioned
  differently.

The observations are assumed to be positive, log-scaled indices
(`Habitat Index`, `Niche Index`); the scripts internally apply `log10(x + 1)`.

## Usage

Set your working directory to this folder, then:

```r
# Step 1 - fit a power equation to every curve
source("power_equation_fitting.R")

# Step 2 - cluster the fitted curves into functional modules
#          (define the maximum number of candidate modules first)
n <- 10
source("tri_functional_clustering.R")

# Step 3 - reconstruct the module interaction network
source("network_reconstruction.R")
```

### Interactive file dialogs

All inputs, intermediate saves and figure/table exports use `file.choose()`,
so R will pop up a file dialog at each step. To run the pipeline
non-interactively (e.g. on a cluster), replace those calls with explicit
paths — for example:

```r
df <- read.csv("data/curves.csv")
save(res1, file = "out/res1.RData")
```

### Tuning

- `n` in step 2 — upper bound on the number of candidate modules. The final
  module count is chosen by BIC.
- `alpha` in `network_reconstruction.R` — ridge penalty of the ODE fit.
- `order <- 2L` — Legendre polynomial order used for the derivative fit.
- Step 3 selects the cluster solution stored as `cluster_result[["cluster2"]]`;
  change this key to analyse another candidate solution.

## Outputs

- fitted trajectories and parameters (`res1` from step 1);
- module assignment and BIC table, plus the clustering figure (step 2);
- network figures and two CSV tables per group (HC / UC / CD): `links`
  (from, to, absolute effect, sign) and `nodes` (independent effect,
  received effect) (step 3).

## Citation

If you use this code, please cite the accompanying manuscript and the
`idopNetwork` package:

> Chen, Y. et al. *idopNetwork: A Network Tool to Dissect Spatial Community
> Ecology*. R package version 0.1.2. <https://cran.r-project.org/package=idopNetwork>

## Contact

Questions about the code or the analysis: open an issue in this repository.
