# Data contract

Data are organized by lifecycle and are never copied into miscellaneous folders.

- `registry/`: source/version/licence records.
- `manifests/`: immutable checksums, row provenance, exclusions, and schemas.
- `raw/`: immutable downloads; ignored by Git.
- `external/`: locked benchmarks; ignored by Git.
- `interim/`: restartable normalized material; ignored by Git.
- `processed/`: small permitted pilot tables and analysis-ready metadata.
- `features/`: rebuildable arrays; ignored by Git.
- `splits/`: immutable split assignments and leakage diagnostics.

No acquisition command overwrites an existing raw snapshot. A source is usable only after its release identifier, official URL, licence/terms, byte size, SHA-256, acquisition time, and processing command have been recorded.

The Gate-0/1 pilot admits only S0/S1 mappings. S0 means the activity and co-complex are the same curated record/construct. S1 means an exact assay construct/sequence and site are mapped to another structure with explicit evidence. UniProt or target-name agreement alone is insufficient.

The canonical bounded pilot is `processed/gate01/rcsb-kd-pilot-strict-v2.tsv`; its paired receptor features are `processed/gate01/rcsb-kd-pocket-views-strict-v2.tsv`. Use only with `manifests/gate01-pilot-strict-v2.json`. The earlier `*-v1` candidate tables are retained because BindingDB reconciliation discovered 12 censored values whose relations were not represented in the RCSB query result. They are non-canonical.

For the completed Gate-2 interaction-feasibility study, the canonical strict data
are `processed/gate02/rcsb-kd-strict-v3.tsv` and
`processed/gate02/rcsb-kd-pocket-views-strict-v3.tsv`, governed by
`configs/gate02/feasibility-effective-v3.yaml`. S2 predicted and S3
pocket-unoccupied experimental views are sensitivity subsets, not substitutes for
the complete S1 primary view. `splits/gate02/prefit-union-components-v3.tsv` was
created before outcome fitting and must be used for any comparison to the reported
Gate-2 result. See `reports/gate02/reproducibility.md` for canonical hashes.

The terminal Gate-3 cohort is `processed/gate03/gate03-strict-kd-v3.tsv`; its
immutable component/scaffold assignments are
`splits/gate03/gate03-component-scaffold-splits-v3.tsv`. The final evaluation
fold ledger is `splits/gate03/gate03-evaluation-folds-v1.tsv`, governed by
`configs/gate03/evaluation-effective-v1.yaml` and amendments 04--05. Final output
hashes and the exact container command are in
`reports/gate03/reproducibility.md`. These results are a terminal hypothesis
no-go and must not be reused as a tuning set for the same architecture proposal.
