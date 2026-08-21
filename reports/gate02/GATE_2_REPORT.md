# iScore3.0 Gate-2 Bounded Interaction Feasibility Report

**Decision date:** 21 August 2026
**Scope:** frozen shallow-control test of whether a receptor pocket adds
reproducible ligand–pocket interaction signal to SMILES-only ligand features and
a fixed protein-sequence encoder. No full iScore3.0 model was trained.

## Executive decision

**NO-GO for full iScore3.0 architecture training under the frozen Gate-2
criterion.** The primary ECFP–ESM-2–S1-pocket tensor model improved pooled
component-OOD RMSE by only 0.002 pKd relative to matched additive fusion; its
95% paired component-bootstrap interval crosses zero. The predeclared practical
margin was 0.10 pKd and gMolAI did not corroborate the interaction direction.

This does not show that pockets are irrelevant. Real-pocket models improve by
about 0.05 pKd over no-pocket and permuted-pocket controls. Rather, the observed
benefit is explainable by additive pocket information, not a reproducible
ligand–pocket interaction term. Cross-attention or a larger geometric model would
therefore be unjustified capacity escalation.

## Gate ledger

| Question | Result | Evidence-based interpretation |
|---|---|---|
| Strict scale and mapping quality | **PASS** | 271 exact uncensored Kd observations, 105 exact constructs, and 21 frozen union components. |
| BindingDB provenance | **PASS** | PDB, ligand, UniProt, exact Kd, and publication fields reconcile; all 10 retained top-score ties were equivalence-audited. |
| Structural leakage graph | **PASS** | 181 validated US-align construct edges were included before fitting. |
| S2/S3 receptor views | **PASS, sensitivity only** | S2: 260 observations/21 components; S3: 144 observations/12 components. |
| gMolAI and ESM-2 adapter contracts | **PASS with exposure caveat** | Hashes, dimensions, determinism, mapping, and label quarantine passed; exact pretraining identity overlap is not enumerable. |
| Sequence and nearest-neighbour controls | **PASS** | ESM-2, sequence identity, ESM cosine, US-align, and descriptor KNN controls used the same nested component protocol. |
| Primary interaction beyond matched additive fusion | **FAIL** | Tensor minus additive RMSE −0.0020 pKd; 95% CI [−0.0074, +0.0044]. |
| Progress to full architecture | **NO** | Frozen progression check is `FAIL`; no full architecture was fitted. |

## Data, views, and information boundary

The canonical strict table is `data/processed/gate02/rcsb-kd-strict-v3.tsv`
(SHA-256 `124a0b342698a8bf2ddb0e64695f5fedd902adf547095d105f4b749f2bc2b5de`).
It has 271 supervised Kd records and 105 label-quarantined historical site
references: 267 unique InChIKeys, 190 Bemis–Murcko scaffolds, 86 UniProt
accessions, and 105 exact construct sequences. pKd spans 2.678–11.301 (median
5.572); resolution spans 0.95–2.80 Å (median 1.80 Å).

The union graph was frozen before any outcome fit. It combines exact
construct/scaffold, ECFP4 Tanimoto ≥0.35, full-sequence identity ≥0.30,
site-sequence identity ≥0.50, shared structure/measurement publication, and
validated structural-similarity relations. It has 21 components but is
materially imbalanced: the largest has 203 observations and 74 constructs.
Pooled RMSE is therefore accompanied by component-macro metrics, component
bootstrap intervals, and largest-component-held-out RMSE.

| View | Receptor source | Coverage | Role |
|---|---|---:|---|
| S0 | query-holo receptor and frozen reference site | 271 observations | sensitivity only; query-holo conformation is privileged |
| S1 | historical holo receptor/site fixed per construct | 271 observations, 21 components | primary view |
| S2 | AlphaFold DB receptor with transferred site | 260 observations, 21 components | predicted-structure sensitivity |
| S3 | mappable pocket-unoccupied experimental receptor | 144 observations, 12 components | incomplete-case apo sensitivity |

No model received query ligand coordinates, conformers, docking poses,
query-derived contacts, or historical reference affinity labels. The pocket input
is a transparent 52-dimensional receptor-only descriptor, not a learned
iScore3.0 pocket encoder.

## Model and evaluation contract

The primary candidate is deliberately low capacity: ECFP4, frozen ESM-2 150M,
and S1 pocket features are additive; ligand and pocket blocks are standardized
and reduced within each fitting fold to 8 and 4 PCs; their 32 outer-product terms
are appended to Ridge regression. All scalers, PCAs, and Ridge/KNN choices are
refit inside nested union-component splits. A deterministic construct-level S1
pocket derangement is the negative control.

Controls include means, nuisance, ligand-only ECFP/gMolAI, pocket-only, ESM-2
sequence-only, additive fusion, ECFP/gMolAI KNN, sequence identity/ESM cosine,
US-align global/pocket KNN, and pocket-descriptor KNN. Sequence and structural
KNN labels are target-balanced.

## Primary S1 results

All values are leave-one-frozen-union-component-out pooled RMSE in pKd.

| Model | RMSE | Component macro | Largest held-out component |
|---|---:|---:|---:|
| global training mean | 1.4490 | 1.2429 | 1.4552 |
| nuisance Ridge | 1.4108 | 1.0412 | 1.4583 |
| S1 pocket Ridge | 1.4165 | 1.1081 | 1.4645 |
| gMolAI Ridge | 1.4378 | 1.1623 | 1.4252 |
| ECFP + ESM-2 additive | 1.5320 | 1.2284 | 1.5411 |
| ECFP + ESM-2 + S1 pocket additive | 1.4843 | 1.1667 | 1.5009 |
| **ECFP + ESM-2 + S1 tensor** | **1.4822** | **1.1651** | **1.4986** |
| ECFP + ESM-2 + permuted S1 tensor | 1.5335 | 1.2313 | 1.5417 |
| gMolAI + ESM-2 + S1 pocket additive | 1.4699 | 1.1034 | 1.4906 |
| gMolAI + ESM-2 + S1 tensor | 1.4704 | 1.1006 | 1.4978 |

No nearest-neighbour shortcut is competitive: ECFP KNN RMSE is 1.6412, ESM
cosine KNN 1.5625, global US-align KNN 1.5614, and pocket US-align KNN 1.5667.

## Interaction and conformation tests

Paired deltas are candidate RMSE minus comparator RMSE; negative favours the
candidate. Intervals use 10,000 component-bootstrap replicates.

| Test | Delta RMSE | 95% interval | Interpretation |
|---|---:|---:|---|
| primary S1 tensor vs matched S1 additive | −0.0020 | [−0.0074, +0.0044] | no incremental interaction evidence |
| primary S1 tensor vs ECFP+ESM-2 no-pocket additive | −0.0498 | [−0.0946, −0.0418] | real pocket information helps, below −0.10 margin |
| primary S1 tensor vs permuted-pocket tensor | −0.0513 | [−0.1094, −0.0405] | real pairing helps, below −0.10 margin |
| gMolAI S1 tensor vs gMolAI S1 additive | +0.0006 | [−0.0352, +0.0075] | corroboration direction fails |
| ECFP S0 tensor vs S0 additive | −0.0018 | [−0.0065, +0.0040] | same null incremental result |
| S2 tensor vs S2 additive | −0.0021 | [−0.0081, +0.0039] | directional point estimate only |
| S3 tensor vs S3 additive | −0.0014 | [−0.0130, +0.0059] | incomplete sensitivity; inconclusive |

S2 retains all 21 components and therefore satisfies the predeclared directional
replication requirement, but its effect is again approximately 0.002 pKd. S3 has
only 12 components and cannot overturn the primary decision.

## Leakage findings and recommendation

Primary folds have zero component and zero construct overlap. Maximum held-out to
training ECFP Tanimoto is 0.3415 and maximum full-sequence identity 0.2724. A raw
pocket TM-score can reach 0.7079 without an edge because the structural rule also
requires the prescribed alignment mode, coverage, aligned-residue count, and RMSD.

Random-pair splitting is strongly optimistic and non-headline: its tensor RMSE is
1.0005 versus 1.4822 component-OOD, while target mean is 1.0678 versus 1.4490.
Each random fold shares 8–14 union components and 35–37 construct groups across
train/test.

Do not train cross-attention, atom–residue message passing, or another
higher-capacity iScore3.0 model from these data. Preserve this as a bounded
negative feasibility result. Any future proposal needs a new committee-approved
phase and new pre-fit freeze, with broader disconnected component coverage while
retaining the exact provenance, structural leakage graph, conformation sensitivity,
matched additive comparator, and permuted-pocket control.

## Canonical evidence

- `configs/gate02/feasibility-effective-v3.yaml`
- `data/processed/gate02/rcsb-kd-strict-v3.tsv`
- `data/splits/gate02/prefit-union-components-v3.tsv`
- `reports/gate02/evidence/bindingdb-tie-equivalence-audit-v1.json`
- `reports/gate02/evidence/gmolai-adapter-audit-v3.json`
- `reports/gate02/evidence/esm2-adapter-audit-v3.json`
- `reports/gate02/evidence/structural-similarity-audit-v3.json`
- `reports/gate02/evidence/leakage-diagnostics-final-v1.json`
- `reports/gate02/evidence/baseline-metrics-v1.json`
- `data/manifests/gate02-baseline-experiment-v1.json`

See [reproducibility.md](reproducibility.md) for commands and expected hashes.
