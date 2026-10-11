# Fixed-horizon release variants

For a fixed supported horizon, the upper empirical variance endpoint of a truthful feasible full release is unchanged after administrative censoring. Use the original cohort size and retain every person. Events exactly at the horizon precede censoring and remain events in the administrative dataset. The function administratively_censor preserves the original endpoint-specific or common censoring relation.

from_horizon_records requires all declared n observations in each margin, with exits no later than the horizon. from_horizon_tables requires exact event-risk rows through the horizon and its exact pre-event risks. Original full-history functions keep their original contracts.

from_terminal_tables instead takes complete event-risk rows strictly before the horizon, plus exact counts X≥tau. from_terminal_records takes every exact exit strictly before the horizon and uses the declared n to derive the terminal counts. All remaining persons must truly have X≥tau; these are not incomplete or selectively dropped head rows. Event status at the boundary is intentionally unreleased. Terminal censoring cells are optimization representatives, not statements that original boundary events were censored. The upper endpoint agrees with the feasible full release. Its lower endpoint belongs to the coarser release.

Separate-censoring lower endpoints can decrease; common-censoring projected score-pairing sets and both extrema agree. No finite-sample nominal coverage improvement or anonymity guarantee follows. Fixed-horizon pointwise inference still needs independent censoring, support, correct cohort membership and regularity. A paired variance or score cross-product is preferable if available.

These adapters retain the original nonnegative-time compatibility model with no supplied hard latent upper bound, delayed entry or weights. A known bounded latent domain needs a correctly clipped compatibility graph. The production bootstrap event filter is a different procedure and is unchanged.

Run python -m unittest discover -s tests -v for all 26 tests, including eight horizon tests. The original code/check.py verifies the 26 stable scientific checksums and its original deterministic suites. No production study is required by these adapters. Stored scientific results and original source files are unchanged.