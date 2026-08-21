# Gate-3 pre-fit data freeze

**Frozen:** 2026-08-21, before any Gate-3 efficacy model was fit

**Protocol:** `configs/gate03/interaction-identifiability-v1.yaml`

**Canonical labels:** `data/processed/gate03/gate03-strict-kd-v3.tsv`

**Canonical splits:** `data/splits/gate03/gate03-component-scaffold-splits-v3.tsv`

## Gate result

The preregistered dataset gate passes. The frozen cohort contains 34 independent
protein/site components, 661 unique exact-Kd observations, and 34 measurement
publications. Twenty-six series contain at least 10 ligands. Eleven series support
22 eligible scaffold-held-out folds under the strict scaffold graph. The largest
component contains 13.01% of observations, and the observed pKd interval is
2.226--10.699.

The selected components contain zero edges under the frozen union leakage graph.
That graph includes sequence, site, global/pocket US-align structure, ligand
identity/similarity, Bemis--Murcko scaffold, and shared-publication relations.
US-align was pinned at source revision
`fa4376bd99fa17a123d05d7ea47cf6574c80d64f` and binary SHA-256
`b5b44c885e61bba1352a20fd8726a4d73f1bd9798e13dc926e315c72b550e4d3`.

## Provenance and label audit

Every component has a strict protein-sequence to structure/entity mapping and an
explicit reference-ligand-defined site. Mapping thresholds were applied before
selection. Exact numeric, positive Kd values are retained; censored values,
ambiguous units, multi-chain target definitions, missing stable publications, and
conflicting duplicate units are excluded by construction.

A stratified primary-source audit covers 12 publications and 24 consensus
measurements. Twenty-two measurements (91.67%) were verified directly, satisfying
the preregistered 90% threshold, and every publication contributes at least one
verified measurement. One CYP51 ketoconazole record disagreed by tenfold with its
primary figure and was quarantined before this freeze. One CBL-B endpoint remains
explicitly unresolved in the accessible primary display; it is retained and
flagged rather than silently treated as verified.

## Machine-readable evidence

- Dataset gate: `reports/gate03/evidence/dataset-gate-audit-v3.json`
- Independent selection: `reports/gate03/evidence/final-independent-selection-audit-v2.json`
- Structural leakage screen: `reports/gate03/evidence/structural-similarity-audit-v2.json`
- Structure mapping: `reports/gate03/evidence/structure-mapping-audit-v2.json`
- Source verification: `reports/gate03/evidence/primary-source-verification-v4.json`
- Label quarantine: `configs/gate03/primary-source-quarantine-v1.yaml`
- Source downloads and checksums: `data/manifests/gate03-primary-sources-v4.json`

Raw BindingDB, coordinate, and publication files remain in ignored `data/raw/`
and `data/interim/` paths. Tracked manifests give their source URLs, sizes, and
hashes; processed labels and split assignments are versioned.

## Boundary after this commit

The next operations are representation generation and preregistered evaluation.
No outcome-dependent dataset, split, endpoint, model-capacity, or progression
criterion change is permitted. Any necessary correction must be recorded as a
dated amendment before the affected fit and reported as a sensitivity analysis.
