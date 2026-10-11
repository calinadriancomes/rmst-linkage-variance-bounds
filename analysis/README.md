# Analysis scripts

`simulate.py` implements the declared sampling laws and time discretization. `estimation.py` computes marginal estimates, the paired reference and compatible variance endpoints. `bootstrap.py` implements the three reported bootstrap comparisons. `rotterdam.py` constructs recorded endpoints and evaluates separate endpoint follow-up at 1,826.25 days. `run.py` coordinates the full analysis, `storage.py` stores one outcome per sample and `summarize.py` reports support, completion and method-specific summaries.

## Reproduce the reported displays

Run `python code/check.py` and `python code/reproduce.py --output outputs` from the repository root. This route uses stored aggregate inputs and small deterministic examples.

## Optional full analysis

```sh
python analysis/run.py --run-full-study --output full_analysis
```

The declared design has 27 scenarios and 39,000 planned original-sample attempts, with seed 20261004. Each bootstrap calculation seeks 499 accepted draws within 4,990 attempts. Original samples are not replaced. An unresolved calculation failure stops execution for review.

The literal and support-conditioned order filters also reject a resampled death margin with no observed events. These attempts are counted as order rejections even when the curves satisfy the survival ordering. This extra rule is retained in the reported comparator; it is not a compatibility rule for R1 or R2.

The reported aggregates retain two failed calculations, one in each uncensored control at zero-based sample indices (0, 1) and (1, 62). The current implementation completes both in targeted checks. A new full execution therefore has its own outcomes and must use a separate output directory. `data/reported_calculation_failures.csv` identifies the reported failures.

New runs store computational state in `analysis_state/`, write `run_summary.json` and describe exported reports in `output_format.md`. `CALCULATION_FAILED` identifies a calculation failure. Valid marginal support is assessed before computation and remains distinct from completion.

The optional Rotterdam route requires source fields obtained locally from survival 3.8-12. See the [data guide](../data/README.md). The R constructor checks cohort counts and endpoint rules; the Python script performs the analysis. A valid paired analysis is preferable when correspondence is available. Neither marginal release provides a privacy guarantee.
