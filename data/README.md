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
