# GWAS Simulation Study: Power Lost to Genetic Model Misspecification

A simulation framework in R for measuring how much statistical power a genome-wide association study (GWAS) loses when a SNP is tested under the wrong genetic model.

Most GWAS pipelines test every SNP under an additive model. If the causal variant actually acts in a dominant, recessive or co-dominant way, that test is misspecified and needs a larger effect to reach genome-wide significance. This project quantifies that cost across a grid of realistic study settings.

> Ongoing research from my PhD in Applied Statistics at Baylor University. The code is public; results are still being written up and are not included here.

## Study design

For each scenario the simulation:

1. **Generates case–control genotypes** at a causal SNP under a chosen *true* model (dominant, co-dominant or recessive), given the minor allele frequency, disease prevalence, causal odds ratio and numbers of cases and controls.
2. **Re-codes and tests** the same data under each competing model: genotypic, dominant, co-dominant, recessive (chi-square tests) and additive (logistic regression).
3. **Estimates power** as the proportion of replicates significant at the genome-wide threshold of 5×10⁻⁸.
4. **Traces the power curve** over the odds ratio for every model, choosing odds ratios adaptively so that evaluations concentrate where the curve is changing.
5. **Summarises the loss of power** as the area between the correctly specified model's power curve and each misspecified model's curve, after smoothing both with monotone shape-constrained splines.

| Factor | Values |
|---|---|
| True genetic model | dominant, co-dominant, recessive |
| Models tested | genotypic, additive, dominant, co-dominant, recessive |
| Minor allele frequency | 0.10, 0.25, 0.40 |
| Disease prevalence | 0.1 to 0.5 |
| Sample size (cases + controls) | 1,000 · 5,000 · 10,000 |
| Significance level | 5×10⁻⁸ |

## Implementation notes

- **S3 object system.** Each genetic model is a constructor returning a classed object, with generics (`generateSNPsCausal`, `encodeSNPs`, `decodeSNPs`, `deepCopy`, `print`) dispatching on the model. Adding a model means adding one file under `scripts/models/`.
- **Adaptive Monte Carlo effort.** Each power estimate starts with 500 replicates, then adds only as many as needed to bring the interval width under a target margin of error, so easy points (power near 0 or 1) cost little.
- **Parallel execution.** Scenarios run across cores with `future` and `furrr`, with progress reporting through `progressr` and reproducible per-scenario seeds.
- **Shape-constrained smoothing.** Power curves are fitted with monotone increasing P-splines (`scam`) before computing the area between curves.

## Repository layout

```
scripts/
  scripts.R            simulation engine: experiment runners, adaptive OR search, power-loss summaries, plots
  models/              one file per genetic model (additive, co-dominant, dominant, genotypic, recessive)
runScriptInTerminal.R  runs the full grid for each true model and sample size
rds_files/             saved simulation output (one file per true model × sample size)
present plots.qmd      Quarto document that turns saved output into figures
Untitled.qmd           scratch notebook for exploring a single result file
```

## Running

```r
install.packages(c("tidyverse", "future", "furrr", "progressr", "gridExtra",
                   "here", "glue", "scam", "patchwork"))
```

```bash
Rscript runScriptInTerminal.R     # full grid; uses all but two cores
quarto render "present plots.qmd" # figures from the saved output
```

A single scenario can be run interactively:

```r
source(here::here("scripts", "scripts.R"))

model     <- ModelDominantCC()
modelComp <- list(ModelGenotypicCC(), ModelCoDominantCC(), ModelAdditiveCC(), ModelRecessiveCC())
model$setParms(n_case = 500, n_control = 500)

powerDf <- runFullExperiment(
  modelRef  = model,
  modelComp = modelComp,
  parmsGrid = tidyr::expand_grid(mar = 0.25, prevalance = 0.3)
)
```

## Skills demonstrated

R · simulation study design · statistical power analysis · statistical genetics · S3 object-oriented programming · functional programming with purrr · parallel computing (future/furrr) · shape-constrained regression · ggplot2 · Quarto

## Author

Nishith Ranjon Roy (Pritom) · [LinkedIn](https://www.linkedin.com/in/pritom65/) · [Google Scholar](https://scholar.google.com.pk/citations?user=HjqFwTUAAAAJ&hl=en)

## License

MIT, see [LICENSE](LICENSE).
