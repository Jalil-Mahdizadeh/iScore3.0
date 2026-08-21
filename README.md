# iScore3.0

iScore3.0 is now a bounded research programme testing whether independently
generated free-ligand conformers add target-dependent affinity signal beyond a
strong ligand-2D baseline. No query-ligand crystal coordinates, docked poses,
protein-ligand contacts, docking teachers, or ligand-to-pocket correspondences
are allowed.

## Current phase

Gate-4A is an interaction-identifiability gate, not a model-building campaign.
It separates three predictive effects:

1. `delta_3d_ligand`: free-ligand 3D beyond frozen gMolAI 2D;
2. `delta_pocket_additive`: receptor information without interaction;
3. `delta_3d_x_pocket`: free-ligand-3D-specific interaction beyond a
   capacity-matched 2D-by-pocket interaction control.

The phase begins with source and censoring audits, deterministic invariant
conformer descriptors, and one frozen pretraining-only ligand-3D representation.
No trainable E(3) encoder is permitted unless Gate-4A passes its preregistered
criteria.

## Start here

- [`docs/gate4a/PROTOCOL.md`](docs/gate4a/PROTOCOL.md) — scientific design and
  progression/termination rules.
- [`docs/gate4a/DATA_GOVERNANCE.md`](docs/gate4a/DATA_GOVERNANCE.md) — fresh-data,
  censoring, split, and test-access policy.
- [`configs/gate4a/protocol-v1.yaml`](configs/gate4a/protocol-v1.yaml) — machine-readable
  frozen settings.
- [`data/README.md`](data/README.md) — lifecycle and provenance contract.
- [`reports/history/GATE3_CLOSURE.md`](reports/history/GATE3_CLOSURE.md) — why the
  former hypothesis was closed and how its full history can be recovered.

## Repository boundary

The active tree retains only reusable gMolAI, ligand mapping, pocket, apo-view,
structure-view, and information-boundary utilities. The terminal Gate-0--3 tree
is recoverable from Git tag `gate3-terminal-2026-08-21`; its data and outcomes
must not be used for Gate-4A architecture selection.

Large source files, container images, external repositories, generated
conformers, features, and checkpoints are ignored by Git. Versioned manifests
record their official URLs, byte sizes, hashes, licences/terms, and processing
commands.

## Status

Gate-4A initiation is complete through the initial source audit and a
candidate-only conformer technical preflight. Davis contains 9,424 exact and 22,400
right-censored tested pairs; none of the censored cells are exact 10,000 nM
labels. Of 72 named compounds, 71 have unaccepted PubChem candidates and one is
quarantined as an unspecified derivative. Deterministic conformer generation
succeeded for 69 candidates; two salt/multicomponent records were correctly
rejected pending a frozen standardization policy.

Exact assay-construct boundaries/sequences are not supplied by the source, all
71 compound mappings still require manual verification, and a locked external
panel has not been selected. Therefore dataset materialization and model fitting
remain forbidden. See [`reports/gate4a/GATE4A_INITIATION_REPORT.md`](reports/gate4a/GATE4A_INITIATION_REPORT.md).
There are no Gate-4A performance results and no model suitable for prospective,
clinical, or regulatory use.
