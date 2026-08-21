# iScore3.0 Gate-3 Interaction Identifiability Report

**Decision date:** 21 August 2026

**Scope:** bounded, preregistered test of low-capacity no-pose ligand--pocket
interaction signal. No cross-attention or full iScore3.0 architecture was fit.

## Executive decision

**HYPOTHESIS-NO-GO. Terminate the current no-pose ligand--pocket interaction
architecture hypothesis; do not escalate model capacity.**

The data, provenance, representation, and leakage gates passed, but the primary
bilinear interaction did not show statistically supported and practically useful
gains jointly across absolute component-OOD affinity and within-target
scaffold-OOD ranking. The result did not replicate across the interpretable
pocket descriptor and pretrained ESM-IF1 structure encoder, failed the real-versus-
permuted pairing criterion, failed one of two S2 non-holo replications, and
reversed under the gMolAI complete-case corroboration. FiLM was unstable in
component-OOD prediction and was preregistered as unable to rescue a failed
bilinear primary result.

This is not evidence that receptor structure is universally irrelevant. Some
descriptor/S2 and small S3 sensitivities are favourable. They are inconsistent
across encoders, endpoints, and receptor views and therefore do not identify a
reproducible ligand-specific pocket interaction.

## Gate ledger

| Question | Result | Evidence-based interpretation |
|---|---|---|
| Dataset depth and independence | **PASS** | 661 exact-Kd ligands in 34 disconnected protein/site components; 26 series have at least 10 ligands. |
| Assay and primary-source provenance | **PASS with flags retained** | 22/24 stratified measurements verified; one tenfold discrepancy quarantined and one unresolved display retained as flagged. |
| Structural and chemical leakage controls | **PASS** | Zero cross-component edges in the frozen union graph; all 120 evaluation folds passed their boundary checks. |
| Receptor and ligand adapters | **PASS** | S1 34/34, S2 32/34, S3 7/34; ESM-2, centred ESM-IF1 v2, and gMolAI audits passed. |
| Absolute component-OOD interaction gain | **FAIL** | Neither S1 receptor encoding has a supported RMSE/MAE gain over its matched additive control. |
| Within-target scaffold-OOD interaction/ranking gain | **FAIL** | Descriptor gains are uncertain; ESM-IF1 significantly worsens centred RMSE. |
| Pocket-pairing specificity | **FAIL** | The real-versus-permuted direction is not jointly favourable for either S1 encoding. |
| Non-holo replication | **FAIL** | S2 descriptor is favourable, but S2 ESM-IF1 is not; S3 absolute and ranking directions conflict. |
| gMolAI corroboration | **FAIL** | Absolute interaction deltas reverse for both encodings, significantly for ESM-IF1. |
| Full architecture progression | **NO** | Frozen progression state is `HYPOTHESIS-NO-GO`; no full model was fit. |

## Frozen cohort and views

The cohort contains 661 unique exact-Kd observations, 34 measurement publications,
and 34 series/components. pKd spans 2.226--10.699 (median 6.252); the largest
component is 13.01% of observations. Eleven series provide 22 eligible scaffold
folds and 158 held-out observations. The strict union graph includes sequence,
site, global/pocket US-align structure, exact ligand, ECFP4, Bemis--Murcko, and
shared-publication edges; no retained component pair is connected by any rule.

| View | Series/components | Component-OOD observations | Scaffold folds / observations | Role |
|---|---:|---:|---:|---|
| S1 historical fixed holo | 34 | 661 | 22 / 158 | primary |
| S2 AlphaFold DB non-holo | 32 | 628 | 21 / 156 | required replication |
| S3 strict pocket-unoccupied X-ray | 7 | 88 | 4 / 12 | sensitivity only |

The two frozen receptor representations are a 116-dimensional rigid-invariant
pocket geometry/chemistry descriptor and a 1,024-dimensional residue-level
ESM-IF1 site summary. ECFP4 is the primary ligand representation; gMolAI covers
657/661 ligands and is complete-case corroboration. ESM-2 supplies the competitive
sequence representation. No query-ligand 3D coordinates, conformers, contacts,
docking poses, or co-complex query geometry enter any model.

## Essential controls

S1 absolute results are leave-one-frozen-component-out. Scaffold metrics are
macro-averaged within target after training-only centring.

| S1 control | Absolute RMSE / MAE | Scaffold RMSE / Spearman | Interpretation |
|---|---:|---:|---|
| global training mean | 1.824 / 1.508 | 1.070 / not defined | intercept reference |
| global-structure KNN | **1.800** / 1.459 | 1.070 / 0.067 | best non-gMolAI absolute control |
| ESM-2 sequence Ridge | 1.869 / 1.480 | 1.070 / -0.258 | competitive sequence control |
| ECFP Ridge | 2.001 / 1.667 | 1.395 / 0.028 | ligand-only linear control |
| ECFP KNN | 2.039 / 1.624 | **1.006** / 0.057 | best centred-error scaffold control |
| pocket descriptor Ridge | 2.178 / 1.809 | 1.070 / -0.159 | pocket-only control |
| ESM-IF1 Ridge | 2.042 / 1.709 | 1.070 / -0.159 | structure-only control |
| ECFP + ESM-2 additive, no pocket | 1.900 / 1.545 | 1.392 / 0.041 | no-pocket fusion control |
| gMolAI Ridge, complete case | 1.843 / 1.475 | 1.096 / 0.181 | secondary ligand-only control |
| gMolAI + ESM-2 + ESM-IF1 additive | **1.729** / 1.369 | 1.134 / 0.270 | best S1 absolute model, complete case |

Absolute performance is generally weak: even the primary ECFP interaction
models are worse than the global mean and global-structure KNN. This reinforces
the no-go but is not used as a substitute for the matched interaction test.

## Primary interaction tests

Deltas are bilinear candidate minus its capacity-matched additive control.
Negative RMSE/MAE and positive rank deltas favour interaction. Intervals are
10,000 paired bootstrap replicates over frozen components (absolute) or assay
series (scaffold).

| S1 receptor | Endpoint | RMSE delta (95% interval) | MAE or Spearman delta (95% interval) | Result |
|---|---|---:|---:|---|
| descriptor | absolute | +0.049 [-0.071, +0.179] | MAE +0.055 [-0.077, +0.187] | no support; wrong point direction |
| descriptor | scaffold | -0.056 [-0.142, +0.030] | Spearman +0.034 [-0.218, +0.263] | practical rank point gain, unsupported |
| ESM-IF1 | absolute | -0.077 [-0.317, +0.204] | MAE -0.050 [-0.301, +0.214] | unsupported |
| ESM-IF1 | scaffold | **+0.256 [+0.026, +0.467]** | Spearman -0.192 [-0.701, +0.368] | supported harm in centred RMSE |

The corresponding S1 scaffold models are also qualitatively inconsistent:
descriptor bilinear has centred RMSE 1.011 but near-random Spearman -0.007 and
concordance 0.487; ESM-IF1 bilinear has RMSE 1.339, Spearman -0.191, and
concordance 0.405. FiLM component-OOD RMSE is 7.148 (descriptor) and 5.967
(ESM-IF1), demonstrating severe instability rather than a credible rescue.

## Replication, permutation, and gMolAI

- S2 descriptor is favourable: absolute RMSE delta -0.020
  [-0.121, +0.093] and scaffold RMSE delta -0.126
  [-0.272, -0.019], although its Spearman interval still crosses zero.
- S2 ESM-IF1 reverses the absolute and scaffold error directions (+0.108 and
  +0.134 RMSE); required dual-encoding non-holo replication therefore fails.
- S3 has only seven components and four scaffold folds. Both encodings improve
  absolute RMSE by about 0.263, but both worsen scaffold RMSE (+0.703 descriptor,
  +0.236 ESM-IF1). It cannot support interaction identifiability.
- S1 real-versus-permuted pairing fails the predeclared joint direction for both
  receptor encodings. A favourable error contrast in the descriptor scaffold
  view is accompanied by worse ranking, while ESM-IF1 is inconsistent.
- gMolAI reverses the absolute S1 interaction direction for both encodings:
  +0.100 RMSE for descriptor and +0.301 [+0.103, +0.507] for ESM-IF1. Its
  scaffold descriptor point estimate is favourable but unsupported; ESM-IF1 is
  harmful. The required corroboration fails.

## Leakage diagnostics and uncertainty

All 34 S1 and 32 S2 component folds have zero component and series overlap.
Scaffold folds intentionally share one assay series but share no frozen scaffold
cluster across the boundary; their maximum ECFP4 similarity is 0.3448, below the
0.35 connectivity threshold. Maximum component-OOD ECFP4 similarity is 0.3455.
Every fold passes the encoded boundary rule. Transformations and target means are
fit on outer-training rows only; scaffold hyperparameter selection uses only
independent inner scaffold clusters after amendment 04.

Two result-free runs exposed implementation edge cases. Amendment 05 records
them transparently: partial fits occurred, but no predictions, metrics, or result
artifacts were written or inspected. The final revision omits only a gMolAI
complete-case scaffold sensitivity fold when no independent inner fold exists;
it does not use endpoint values or change primary S1/S2/S3 fold counts.

## Risks and interpretation limits

- Only 11 of 34 series support scaffold-OOD evaluation; S3 ranking rests on 12
  observations. These sensitivities are underpowered and cannot rescue failure.
- Exact Kd and source mapping are unusually strict, but 34 publications still
  contain laboratory and assay-context heterogeneity not fully representable in
  public metadata.
- Frozen pretrained encoders cannot provide an exhaustive enumerable audit of
  all pretraining identity exposure. Primary ECFP conclusions do not depend on
  gMolAI, and all label-bearing downstream fits obey strict folds.
- The experiment rejects this architecture hypothesis at this evidence level;
  it does not prove a mathematical impossibility for every no-pose formulation.

## Terminal recommendation

Do not train cross-attention, atom--residue message passing, a larger FiLM model,
or another capacity escalation under the current no-pose ligand--pocket
interaction hypothesis. Archive Gate-3 as a bounded negative result. Any future
work must begin from a materially different, committee-approved scientific
hypothesis and a new prefit contract; it must not be presented as tuning or
rescuing iScore3.0 after this gate.

Canonical machine-readable evidence is
`reports/gate03/evidence/evaluation-metrics-v1.json`,
`reports/gate03/evidence/leakage-diagnostics-v1.json`, and
`data/manifests/gate03-evaluation-v1.json`. See
[reproducibility.md](reproducibility.md) for the exact command and hashes.
