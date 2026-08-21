# Gate-4A data contract

Data are organized by lifecycle and never copied into miscellaneous folders.

- `raw/gate4a/`: immutable publisher or repository downloads; ignored by Git.
- `manifests/gate4a/`: versioned URLs, terms, byte sizes, SHA-256 hashes, schemas,
  and acquisition commands.
- `interim/`: restartable normalized material; ignored by Git.
- `processed/gate4a/`: small, auditable normalized metadata only.
- `features/gate4a/`: rebuildable conformer and encoder arrays; ignored by Git.
- `splits/gate4a/`: immutable scaffold/target-component assignments and fold
  ledgers created before fitting.
- `external/`: locked confirmatory panels; ignored by Git and access-controlled.

No acquisition or processing command may overwrite an existing snapshot. A
source is unusable until its official URL, release identifier, byte size,
SHA-256, acquisition time, terms, and processing command have been recorded.

The PubChem name-resolution snapshot remains immutable candidate evidence. The
final tracked identity ledger is
`processed/gate4a/davis-compound-identity-final-v2.tsv`: 69 parent structures are
admitted after completed adjudication, technical validation and explicit project-
owner authorization. BMS-345541 and JNJ-28312141 follow the explicit publication-
parent/counterion-removal rule; no other multicomponent record is silently desalted.
The earlier 69-row review packet remains immutable historical provenance but is
superseded; no additional external/manual chemical signature is required.

The completed isolated experiment's rebuildable arrays are
`features/gate4a/delta3d-ligand-v1.npz`, raw free/control ensembles are
`interim/gate4a/delta3d-free-ensembles-v1.jsonl.gz`, and OOF predictions are
`interim/gate4a/delta3d-ligand-oof-v1.npz`. They remain Git-ignored; tracked
manifests under `reports/gate4a/evidence/` record hashes, dimensions and provenance.

Coordinate provenance is tracked in
`processed/gate4a/alphafold-pocket-admission-v1.tsv` and the corresponding manifest.
The 85 KLIFS columns include explicit alignment gaps: 288 estimands have exact
coordinates for every non-gap column, while 50 fail closed. Apo/non-holo decisions
are in `processed/gate4a/apo-view-admission-v1.tsv`. Pairwise structural similarities,
frozen structural edges, and their transitive receptor components live in
`processed/gate4a/alphafold-pocket-structural-similarity-v1.tsv` and
`splits/gate4a/`.

## Affinity observations

Every pair must retain its source-cell semantics:

- `exact`: a positive numeric Kd with its original unit;
- `right_censored_kd`: Kd is greater than a documented assay limit;
- `missing`: not assayed or genuinely unavailable;
- `invalid`: present but unparseable or contradictory.

In Davis 2011 and Karaman 2008, a blank affinity cell means weak binding above
10 micromolar or no detection in the 10 micromolar primary screen. It is
therefore right-censored Kd, equivalently left-censored pKd at 5, and must never
be stored as an exact 10,000 nM observation.

## Governance

Davis and Karaman are development/technical-replication resources. They cannot
serve as the sole confirmatory test because of target, chemotype, and assay
overlap. A distinct locked panel is required for a progression decision. The
terminal Gate-3 cohort is quarantined by policy and absent from the active tree.
The OKL and KIRHub raw workbooks are confirmation-candidate evidence only;
aggregate audits are tracked, but no external cell-level label is copied into Git
or released to a modeling process. Outcome-blind eligibility ledgers under
`processed/gate4a/confirmation/` are frozen and currently release zero strict pairs.
