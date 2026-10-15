# Marginal survival tables and paired-reference calibration

These examples supplement the original simulation results. The 39,000 original attempts, their outcomes and the stored tables and figures are unchanged.

## Rotterdam marginal-table example

The files in `rotterdam_terminal_example/` were constructed from the publicly distributed Rotterdam breast-cancer data in R **survival** 3.8-12. Both margins refer to the same 2,982 patients at a fixed horizon of 1,826.25 days, with separate endpoint censoring. The progression/death and death tables contain 878 and 600 distinct pre-horizon event times. Their exact pre-event risk counts at the horizon are 1,571 and 2,084.

Each CSV gives the exact event times strictly before the horizon, the number at risk just before each event time and the number of events. The metadata specify the cohort and endpoint conventions. Event status at the horizon and subsequent follow-up are not included. This is a *terminal-symbol R2 release*, not a complete administrative event-risk table. The same-cohort membership and censoring conventions were verified from the original source when the tables were constructed; they cannot be established from the separate margins alone.

From the repository root, run:

```sh
python extensions/analyst_from_tables.py
```

The analyst uses only the two marginal CSV files and `metadata.json`. The attained upper empirical variance is 105.01687556404875 days squared and the associated same-center 95% interval is 214.20322028324568 to 254.37375689690083 days. These agree with the full-history R2 results. The lower variance bound is specific to the released information and is not generally preserved by terminal-symbol reduction. Reference values appear in `expected_results.json`.

The optional `produce_marginal_tables.R` constructs the two CSVs from the versioned **survival** source archive. That script runs on the source-data side and is not needed for the analyst calculation. This is a controlled release of margins derived from a publicly distributed linked dataset, not evidence of a routinely published unlinked clinical release or a privacy guarantee. See [third-party terms](../LICENSES/third_party.md) for dataset attribution.

## Additional paired-reference calibration

The separate diagnostic in `fresh_calibration_results.json` contains 20,000 samples in each of two n=500 settings, with seed 2026101317. Paired-reference coverage was 18,957/20,000 (94.785%) and 18,938/20,000 (94.690%). These settings were examined after the original 933/1,000 result. They do not replace the original samples or provide new R1 or R2 coverage estimates.

Run `python extensions/paired_calibration.py` to repeat the fixed-seed diagnostic. To test the example and its input checks, run `python extensions/test_extensions.py` (the numerical reference test solves the full matching problem). Python dependencies are listed in the root `requirements.txt`. R is needed only for the optional source-side construction.
