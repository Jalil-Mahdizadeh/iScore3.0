# Gate-3 configuration

`interaction-identifiability-v1.yaml` is the pre-fit scientific contract for the
bounded Gate-3 phase. It was frozen before any Gate-3 label/model fit. Dataset
census, chemistry validation, structure mapping, and manual source verification
are data curation, not model-outcome inspection.

The audit trail intentionally retains every dated decision:

- `protocol-amendment-01.yaml` fixes the deterministic, bounded maximum
  independent-set selection procedure used after structural leakage screening.
- `protocol-amendment-02.yaml` defines scaffold-OOD coverage reporting without
  weakening the preregistered scaffold graph or its fold-local isolation.
- `manual-primary-source-audit-v*.yaml` records the successive, pre-fit source
  strata, including candidates that could not be verified from accessible
  primary material.
- `primary-source-verification-v4.yaml` is the final stratified audit registry.
- `primary-source-quarantine-v1.yaml` excludes the one confirmed BindingDB versus
  primary-publication label discrepancy before any efficacy fit.
- `protocol-amendment-03.yaml` makes gMolAI a complete-case corroborating
  sensitivity after four disconnected ion-pair inputs fail unchanged encoding.
- `protocol-amendment-04.yaml` replaces an invalid component-grouped scaffold
  inner split with training-only within-series scaffold centring.
- `protocol-amendment-05.yaml` transparently excludes only a view/encoder-specific
  scaffold sensitivity fold when no independent inner fold exists. Two aborted
  runs wrote no result artifacts and exposed no model outcomes.

The canonical pre-fit cohort is `data/processed/gate03/gate03-strict-kd-v3.tsv`.
Its component/scaffold assignments are frozen in
`data/splits/gate03/gate03-component-scaffold-splits-v3.tsv`. See
`reports/gate03/PREFIT_DATA_FREEZE.md` for the signed-off census, evidence paths,
and the boundary between curation and efficacy modelling.

Gate-3 has three terminal states:

- `DATA-NO-GO`: the minimum depth, independence, mapping, or audit contract is
  not met; no efficacy models are fit.
- `HYPOTHESIS-NO-GO`: the data gate passes, but low-capacity interactions do not
  meet the joint statistical, practical, replication, and specificity criteria.
- `GO`: all predeclared criteria pass; only then may a separate committee-approved
  phase propose a larger architecture.

No cross-attention, atom--residue message passing, docking, ligand conformer, or
ligand-coordinate model is permitted in this phase.

Gate-3 completed as `HYPOTHESIS-NO-GO`. The current no-pose ligand--pocket
interaction architecture hypothesis is terminated. See
`reports/gate03/GATE_3_REPORT.md` and `reports/gate03/reproducibility.md`; no
larger architecture was trained.
