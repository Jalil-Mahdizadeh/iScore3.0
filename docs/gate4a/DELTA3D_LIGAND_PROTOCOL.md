# Gate-4A isolated Delta3D-ligand protocol

Status: frozen before outcome generation on 2026-08-21. Only the ligand main effect is released. Receptor-additive and ligand-receptor interaction fitting on Davis remain blocked.

## Scientific estimand

Each ligand-only predictor emits the same value for all 338 targets for a given ligand. It therefore tests whether free-conformer information improves prediction of a held-out ligand's overall Davis pan-kinase potency/promiscuity distribution. It cannot predict target preference and cannot establish ligand-pocket complementarity. A positive result does not reverse the receptor structural-component collapse or unblock the failed no-pose interaction hypothesis.

The 69 admitted ligands are final under the project owner's explicit authorization. The prior secondary-QA packet is retained as historical provenance but its pending external signature is superseded; no further external/manual sign-off is required.

## Labels and split

The source matrix is parsed without converting inactivity censoring to exact Kd. Numeric cells are exact pKd; blank cells are right-censored Kd at 10,000 nM, equivalently pKd <= 5. The primary score is held-out-cell Gaussian Tobit negative log likelihood. Exact-only RMSE/MAE and within-target held-out-ligand Spearman/NDCG are supportive. Within-ligand target ranking is omitted because every ligand-only model is constant over targets.

All rows for a ligand remain together. The immutable ligand graph uses exact Bemis-Murcko equality or Morgan-radius-2/2048 Tanimoto >= 0.60, followed by transitive closure. Ten label-blind outer folds and five inner folds are assigned at component level. Random pair splits are forbidden. Model selection occurs only inside the outer training set.

## Frozen representations and capacity

- `M2D`: released frozen gMolAI 384-vector, strictly generated from the final parent SMILES.
- `M2D+Det3D`: gMolAI plus deterministic invariant shape, pharmacophore-distance, energy and ensemble-diversity features. The topology-only pharmacophore-presence block is excluded.
- `M2D+UniMol3D`: gMolAI plus frozen Uni-Mol v1 molecule-all-H CLS representations generated from exactly the same free conformers. Uni-Mol is never affinity-finetuned.
- 3D-only versions are diagnostics, not replacements for the preregistered nested contrasts.

Uni-Mol v1 is used instead of the previously tentative, disabled Uni-Mol2 entry. This substitution is frozen before feature extraction because the official Uni-Mol2 implementations have unresolved cross-interface/batching reproducibility reports, whereas the mature v1 API accepts explicit atom-coordinate inputs. The exact upstream commit, checkpoint and dictionary hashes must be recorded before inference.

The pre-outcome GPU audit found exact repeats for an identical batch, but maximum absolute CLS differences of 0.0047 when the same conformers were inferred singly rather than batched and 0.0034 after a rigid transform plus re-centring. These are floating-point execution effects, not a change of molecular distances. Extraction therefore freezes batch size 64 and ledger/conformer order. A 0.006 maximum-absolute numerical tolerance is recorded; no post-outcome batching changes are permitted.

Each raw branch is independently median-imputed, standardized and projected to 32 dimensions using training-fold data only. PCA signs use the largest-absolute-loading-positive rule. The common 32-dimensional 3D branch keeps augmentation capacity matched across actual and negative-control conditions. A linear Gaussian Tobit head with a globally fitted scale and nested, group-held-out L2 selection is intentionally low capacity.

The baseline design is 32-dimensional and each augmented design is 64-dimensional (32 gMolAI plus 32 candidate/control dimensions). Thus an augmented-versus-baseline gain alone is insufficient: actual, destroyed and topology-fake branches have exactly the same added capacity. L2 is selected by minimum inner held-out cell NLL from the frozen seven-value grid; exact ties prefer the larger penalty. The intercept is unpenalized and Gaussian scale is bounded broadly to 0.05--5 pKd for numerical stability.

## Geometry controls

Five complete ETKDGv3/MMFF94s generations use seeds 20260821--20260825. For each seed:

- actual ensemble uses every converged conformer within 10 kcal/mol of the minimum;
- single uses only the minimum-energy conformer;
- energy permutation applies a deterministic nonidentity permutation between the existing coordinate and energy lists;
- coordinate destruction replaces coordinates by seeded isotropic Gaussian clouds, centred and scaled to each actual conformer's heavy-atom radius of gyration, preserving atom identity/count, ensemble size, energy multiset and gross scale while destroying covalent geometry;
- topology fake 3D embeds unweighted graph-shortest-path distances by deterministic classical multidimensional scaling. It contains no free-conformer geometry.

For Uni-Mol, actual and coordinate-destroyed inputs use the identical frozen encoder and aggregation. Single versus ensemble and repeat-seed stability are evaluated. Any apparent gain reproduced by fake/destroyed geometry is a capacity, topology, or regularization effect rather than evidence for free 3D.

No bound/crystallographic ligand coordinates, receptor data, docking, bioactive-conformer selection, complex contacts, ligand-pocket correspondence, or test-guided conformer/checkpoint selection may enter this experiment.

## Decision rule

No numeric practical-equivalence margin is available from defensible repeat data, so none is invented. Reproducible incremental free-3D information requires all of the following:

1. The actual-ensemble addition has a positive primary NLL gain over gMolAI with a paired ligand-component bootstrap 95% interval above zero.
2. The primary gain is positive for at least four of five independently generated conformer seeds.
3. Actual deterministic geometry outperforms coordinate destruction and topology-fake 3D on the primary metric.
4. Frozen Uni-Mol actual-coordinate representations outperform the identically encoded coordinate-destruction control and corroborate the direction.
5. Exact-error and within-target-ranking results show no material directional contradiction.

Failure of any requirement is reported as no reproducible Delta3D-ligand information for this gate. Because no noise-derived equivalence margin exists, a pass establishes reproducible predictive information, not necessarily practical importance.

## Locked receptor decision

On Davis, `Delta_pocket-additive` and `Delta3D x pocket` remain `BLOCKED`. The existing broad receptor-structural-novelty graph and its 323/338 dominant component are immutable. They will not be relaxed or redefined after outcomes.
