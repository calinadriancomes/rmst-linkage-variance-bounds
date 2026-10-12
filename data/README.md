# Data guide

The CSV and JSON files contain stored aggregate results, scenario definitions and population targets. Missing entries denote undefined or inapplicable quantities, not zero. Display reproduction reads these results without generating trial samples or fitting the Rotterdam analysis.

## Denominators

`denominators.csv` contains one row per scenario. `planned` is the number of original samples in the design; `attempted` counts samples attempted. `support_valid` counts samples with valid marginal support, including supported samples whose calculations failed. `completed` counts completed primary calculations. `support_failed`, `calculation_failed`, `not_attempted` and `running` describe the remaining outcomes. `completion_proportion` is completed/attempted. In each uncensored control, support_valid is 1,000, completed is 999, calculation_failed is 1 and completion_proportion is 0.999.

`bias`, `bias_MCSE` and `RMSE` describe the point estimator among completed calculations. `RMSE_MCSE_delta_method` is its reported Monte Carlo standard error and `RMSE_MCSE_denominator` is the completed count used for it. `true_structural_infeasibility_count` counts structurally infeasible samples. `support_includes_calculation_failures` states that valid support does not require successful computation.

In `method_summaries.json`, `completed` counts completed primary calculations and `method_completed` counts completed intervals for the specified method. Coverage conditions on method_completed. Each width, gain or ratio distribution carries its own finite-value count. Covered/planned is operational yield and differs from bootstrap completion/selected yield. `main_design_table.csv` uses `calculation_failed` for computation failures. `reported_calculation_failures.csv` identifies the two failed calculations retained in the reported aggregates.

## Rotterdam source

The breast-cancer data are publicly distributed with the R survival package. Obtain the exact survival 3.8-12 source archive from [CRAN](https://cran.r-project.org/src/contrib/survival_3.8-12.tar.gz) or its [version archive](https://cran.r-project.org/src/contrib/Archive/survival/survival_3.8-12.tar.gz), SHA-256 `ebb70d9a4e7e8ad7f886a563184d6f079e64c0591a0bce219ae90967f0bdf994`. Use this version to check the reported endpoint construction.

```sh
Rscript data/rotterdam_construct.R survival_3.8-12.tar.gz rotterdam_check.csv
```

All 2,982 persons are retained at 1,826.25 days. The progression-free endpoint uses recorded recurrence when present, an observed death before or at the recurrence follow-up boundary otherwise and that earlier boundary as a censored exit when neither applies. The death endpoint retains its own follow-up. All 43 deaths after earlier recurrence censoring remain included; recurrence follow-up is not extended to those deaths. The paired reference uses recorded correspondence rather than latent recurrence ground truth. No treatment or subgroup analysis is included.

For the optional Python analysis, verify the archive checksum above, run the constructor check and load the same archive in the R session that exports the source columns:

```r
archive <- "survival_3.8-12.tar.gz"
folder <- tempfile()
dir.create(folder)
untar(archive, exdir = folder, tar = "internal")
desc <- read.dcf(file.path(folder, "survival", "DESCRIPTION"))
stopifnot(desc[1, "Version"] == "3.8-12")
source_data <- new.env()
load(file.path(folder, "survival", "data", "cancer.rda"), envir = source_data)
z <- source_data$rotterdam
stopifnot(nrow(z) == 2982, length(unique(z$pid)) == 2982)
write.csv(z[c("pid", "rtime", "recur", "dtime", "death")],
          "rotterdam_source.csv", row.names = FALSE)
unlink(folder, recursive = TRUE)
```

The optional analysis then runs as `python analysis/rotterdam.py rotterdam_source.csv rotterdam_analysis.json`. It estimates the recorded-data results; it is not part of the aggregate-display reproduction route.

The constructor's count check is separate from this source CSV and from variance estimation. Person-level source data and the third-party archive are not distributed in this repository. Dataset attribution appears in the survival package manual and the article's references.

## Input requirements

The method requires complete unlinked endpoint records (R1) or full exact event-risk tables with the horizon checkpoint (R2) from the same known cohort. Sparse published Kaplan–Meier curves are insufficient. Scenario definitions and population targets describe the declared simulation inputs; the [analysis guide](../analysis/README.md) describes optional full execution.

The original simulations evaluate the complete full-history R1 and R2 releases. For a prespecified, data-supported fixed horizon, administrative censoring preserves the attainable upper empirical-variance endpoint using exact through-horizon information and the required horizon checkpoint. Under separate censoring, the lower endpoint may decrease; under genuinely common censoring, both endpoints are preserved. See [fixed-horizon release guidance](../HORIZON_RELEASE_NOTES.md) for the release variants and their conditions.
