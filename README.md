# iScore3.0

iScore3.0 Gate-4A is closed with a **no-go** decision. Under strict ligand-component
OOD evaluation, independently generated free-ligand conformers did not add
reproducible predictive information beyond frozen gMolAI 2D. No query-ligand
crystal coordinates, docked poses, protein-ligand contacts, docking teachers, or
ligand-to-pocket correspondences were used.

## Current phase

Gate-4A was an interaction-identifiability gate, not a model-building campaign.
It separates three predictive effects:

1. `delta_3d_ligand`: free-ligand 3D beyond frozen gMolAI 2D;
2. `delta_pocket_additive`: receptor information without interaction;
3. `delta_3d_x_pocket`: free-ligand-3D-specific interaction beyond a
   capacity-matched 2D-by-pocket interaction control.

The isolated `delta_3d_ligand` experiment is complete and failed the joint
preregistered progression rule. `delta_pocket_additive` and
`delta_3d_x_pocket` were not trained and remain blocked on Davis. No trainable
E(3), interaction, or cross-attention model is authorized.

## Start here

- [`docs/gate4a/PROTOCOL.md`](docs/gate4a/PROTOCOL.md) — scientific design and
  progression/termination rules.
- [`reports/gate4a/GATE4A_DELTA3D_LIGAND_REPORT.md`](reports/gate4a/GATE4A_DELTA3D_LIGAND_REPORT.md)
  — final result, controls, limitations, and no-go recommendation.
- [`docs/gate4a/DATA_GOVERNANCE.md`](docs/gate4a/DATA_GOVERNANCE.md) — fresh-data,
  censoring, split, and test-access policy.
- [`configs/gate4a/protocol-v1.yaml`](configs/gate4a/protocol-v1.yaml) — machine-readable
  frozen settings.
- [`reports/gate4a/GATE4A_DATASET_ADMISSION_REPORT.md`](reports/gate4a/GATE4A_DATASET_ADMISSION_REPORT.md)
  — original dataset-admission decision.
- [`reports/gate4a/GATE4A_PROVENANCE_CLOSURE_REPORT.md`](reports/gate4a/GATE4A_PROVENANCE_CLOSURE_REPORT.md)
  — historical pre-experiment provenance decision, now superseded where noted.
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

The 69 Davis ligand identities are final based on the completed adjudication,
technical validation, and explicit project-owner authorization. The historical
secondary-QA packet is retained but superseded; no additional external/manual
chemical sign-off is required.

Canonical AlphaFold v6 pocket coordinates pass exact sequence/coordinate checks
for 288/338 standardized receptor estimands. Fifty fail closed because an exact
KLIFS-to-canonical mapping was not established. Thirty-four structures pass the
strict global-zero-nonpolymer apo definition and 39 pass the audited binding-site-
unoccupied tier. The predeclared structural-leakage union produces only 10 target
components and places 323/338 targets in one component; this makes a credible
double-cold Davis interaction experiment non-identifiable under the frozen rule.
Outcome-blind OKL/KIRHub eligibility ledgers are frozen but release zero strict
confirmation pairs. Public repeat evidence still cannot support numeric metric-
specific practical-equivalence margins.

The isolated experiment evaluated 51 fully nested, censor-aware models across 66
ligand leakage components and five independent conformer seeds. Deterministic 3D
had pooled NLL gain -0.0026 versus gMolAI (95% component-bootstrap interval
-0.0178 to 0.0136); frozen Uni-Mol v1 had -0.0136 (-0.0311 to 0.0041). Actual
Uni-Mol coordinates were worse than destroyed-coordinate controls. The free-
conformer pivot therefore terminates under the frozen rule. The separate future
kinase-selectivity design is archived as design-only and is not authorization to
train.
