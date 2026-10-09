# Variance bounds for time alive after progression

This repository accompanies *Variance bounds for restricted mean time alive after progression from exact marginal survival information*. It contains the method implementation, analysis scripts, compact examples and stored aggregate results used to reproduce the article's main tables and figures.

The article considers two exact marginal information settings: complete unlinked endpoint records (R1) and complete event-risk tables from a known common cohort (R2). The reproduction notebook checks the method on small examples and rebuilds the reported displays from stored aggregate results.

## Reproduce the reported results

[Open the reproduction notebook in Colab](https://colab.research.google.com/github/calinadriancomes/rmst-linkage-variance-bounds/blob/main/notebooks/reproduce.ipynb) to run it in the browser, or run the notebook from a local checkout.

## Files

- `notebooks/reproduce.ipynb`: reproduction notebook.
- `src/`: marginal estimators and compatible variance bounds.
- `analysis/`: simulation, bootstrap and Rotterdam analysis scripts.
- `code/`: input checks and table and figure reproduction.
- `data/`: aggregate results, scenario definitions and population targets.
- `examples/` and `tests/`: small deterministic examples and unit tests.
- `checksums.sha256`: hashes for stable scientific inputs and stored results.

## Data

The stored data contain aggregate results and small synthetic examples. R1 and R2 require the specified complete margins from the same known cohort; ordinary sparse Kaplan–Meier output is insufficient. See the [data guide](data/README.md) for inputs, denominators and the optional Rotterdam source check.

## Local use

From the repository root:

```sh
python -m pip install -r requirements.txt
python code/check.py
python code/reproduce.py --output outputs
```

These commands reproduce the reported displays from stored results and small examples. The [analysis guide](analysis/README.md) describes the optional full analysis.

## Citation

Călin-Adrian Comes and Vasile Paul Bresfelean. *Variance bounds for restricted mean time alive after progression from exact marginal survival information*. Publication details will be added when available.

## Licenses

Author-written code is MIT licensed. Repository-authored documentation, generated figures and author-generated aggregate materials are CC BY 4.0. Third-party materials retain their upstream terms. See [LICENSE](LICENSE) and [third-party terms](LICENSES/third_party.md).
