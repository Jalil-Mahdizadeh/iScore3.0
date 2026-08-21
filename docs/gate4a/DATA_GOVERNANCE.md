# Gate-4A data and test governance

## Dataset roles

- **Davis 2011:** development source for dense kinase-panel censor-aware methods.
- **Karaman 2008:** technical replication of source semantics; not an independent
  final test because of shared assay technology, targets, and compounds.
- **Locked confirmatory panel:** must be selected and approved before development
  outcomes are reviewed. Prefer a prospective matrix or a distinct protein
  family with non-overlapping scaffolds and receptor structural components.
- **Gate-3:** permanently quarantined for selection.

## Required source audit

For every matrix cell preserve the original publication, table, row/column,
assay technology, construct, relation, value, unit, censor limit, and whether the
cell was assayed. Derivative ML matrices are comparison aids only and cannot be
the source of truth.

Davis Supplementary Table 4 states that blanks are tested pairs with Kd above
10 micromolar or no detection in a 10 micromolar screen. Karaman Supplementary
Table 2 uses the same convention. Both are encoded as `right_censored_kd`.

The initial Davis audit found 31,824 tested pairs: 9,424 exact numeric Kd values
and 22,400 right-censored cells. No numeric source cell equals 10,000 nM; the
largest exact value is 9,900 nM. These counts are frozen in
`reports/gate4a/evidence/davis-source-audit-v1.json`.

## Qualification before feature generation

Freeze and report:

- exact/censored/missing cell counts and matrix completeness;
- unique standardized ligands, Bemis-Murcko scaffolds, stereochemical ambiguity,
  salts, protomers, and tautomers;
- exact protein sequence/construct and mutation mapping;
- binding-site sequence and structural components;
- rotatable-bond, stereocentre, macrocycle, shape, and conformer-diversity
  distributions;
- assay batches, publications, duplicate measurements, and disagreement;
- effective scaffold and target-component sample sizes and a power analysis.

A kinase panel with insufficient ligand-3D or pocket diversity is non-diagnostic
for a general free-conformer hypothesis.

## Current qualification blockers

The Davis metadata provide 442 assay labels linked to 384 RefSeq accessions, but
not the actual recombinant construct sequences or residue boundaries. There are
54 mutant rows and separate phosphorylation-state labels. A full-length RefSeq
sequence is therefore only a candidate reference, not exact assay-construct
provenance.

The publication provides compound names and drawn structures, not machine-readable
SMILES, InChI, or CIDs. PubChem name resolution yielded 71 candidates; all remain
pending manual comparison to the published structure and independent identity
evidence. `BIBF-1120 (derivative)` is quarantined because the derivative is not
identified. BMS-345541 and JNJ-28312141 resolve to multicomponent records and
require a preregistered salt/parent policy.

Davis Table 3 names ruxolitinib's development code `INCB018424`, whereas the
corresponding Table 4 affinity column is `INCB18424`. The row/column position and
alternative name strongly suggest a typographical alias, but the pipeline does
not silently equate them; the manual review ledger flags the discrepancy.

No target structure, pocket, ligand feature, split, or model may be materialized
until these identity decisions and the locked-confirmation governance are frozen.

## Structure boundary

Primary receptor views are true apo structures or sequence-predicted structures
for the exact construct/state. Pockets are defined by family alignment or a
query-independent reference-site policy. Query-ligand coordinates, query-ligand
contact radii, ligand-derived pharmacophores, and complex-derived exclusion
volumes are forbidden.

Holo receptors may appear only in a separately labelled sensitivity analysis and
must never define the pocket used by the primary model.

## Access control

Development and test manifests are separate. The locked test labels should be
held by an independent custodian or encrypted/unavailable to the modeling
process. Test access is one-shot after signed hashes for code, config, features,
splits, checkpoints, and the analysis plan are recorded.
