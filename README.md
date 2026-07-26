# Meta-analysis with `glmmTMB`: simulation study and illustrative examples

This repository contains the code, simulated data, results and worked examples supporting the manuscript **“Meta-analysis with the `glmmTMB` R package”**.

**Supplementary webpage:**
https://coraliewilliams.github.io/equalto_sim_study/webpage.html

## Repository structure

```text
equalto_sim_study/
├── R/
│   ├── 01_sim_params.R      # Defines simulation conditions
│   ├── 02_functions.R       # Simulation functions
│   ├── 03_run_sims.R        # Runs simulations and fits models
│   ├── 04_sim_results.R     # Summarises results and produces figures
│   └── ...                  # Additional example and simulation scripts
├── data/                    # Simulation parameter grids
├── results/
│   ├── raw/                 # Raw simulation and model outputs
│   └── figures/             # Figures generated from the results
├── docs/
│   ├── webpage.qmd          # Source for the supplementary webpage
│   ├── webpage.html         # Rendered webpage
│   └── felix_data/          # Data used in the phylogenetic example
└── equalto_sim_study.Rproj
```

## Simulation study overview

The simulation study evaluates the `equalto` covariance structure in [`glmmTMB`](https://github.com/glmmTMB/glmmTMB), which allows a known sampling-error variance–covariance matrix to be incorporated directly into a model.

Results from `glmmTMB` are compared with equivalent models fitted using [`metafor`](https://wviechtb.github.io/metafor/).

Four effect-size measures are considered:

| Effect-size measure          | Outcome type | Main models compared                                       |
| ---------------------------- | ------------ | ---------------------------------------------------------- |
| Standardised mean difference | Continuous   | `metafor::rma.uni()` and Gaussian `glmmTMB` with `equalto` |
| Log response ratio           | Continuous   | `metafor::rma.uni()` and Gaussian `glmmTMB` with `equalto` |
| Log odds ratio               | Binary       | Two-stage models and binomial GLMMs                        |
| Log incidence rate ratio     | Count/rate   | Two-stage models and Poisson GLMMs                         |

The full design includes 20 studies per meta-analysis, null and moderate overall effects, three heterogeneity levels, moderate and rare event settings, and 1,000 repetitions per condition.

Performance is assessed using convergence, computation time, bias, root mean squared error, confidence-interval width, coverage, Type I error and power.

## Worked examples

The supplementary webpage includes examples of:

* multilevel meta-analysis with correlated sampling errors;
* bivariate meta-analysis;
* network meta-analysis;
* phylogenetic meta-analysis; and
* location–scale meta-analysis.

Where possible, equivalent models are fitted with `metafor` for comparison.

## R package requirements

The main packages required to run and analyse the simulations are the following:

install.packages(c(
  "glmmTMB",
  "metafor",
  "tidyverse",
  "furrr",
  "progressr",
  "data.table"
))

## Reproducing the simulation study

Clone the repository:

```bash
git clone https://github.com/coraliewilliams/equalto_sim_study.git
cd equalto_sim_study
```

Run the scripts from the repository root in the following order:

```r
source("R/01_sim_params.R")
source("R/03_run_sims.R")
source("R/04_sim_results.R")
```

`03_run_sims.R` automatically sources `02_functions.R`.

The full simulation is computationally intensive. For a preliminary test, reduce the number of repetitions in `01_sim_params.R` and adjust the number of parallel workers in `03_run_sims.R`.

## Supplementary webpage

The source file is:

```text
docs/webpage.qmd
```

The rendered version is available at:

https://coraliewilliams.github.io/equalto_sim_study/webpage.html

## Citation

This repository supports the associated manuscript:

> Williams, C. et al. *Meta-analysis with the `glmmTMB` R package*.
> Citation details and DOI will be added following publication.

## Contact

**Coralie Williams**
Email: [coralie.williams@outlook.com](mailto:coraliewilliams7@gmail.com)
