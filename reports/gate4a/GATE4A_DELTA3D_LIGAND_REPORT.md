# Gate-4A isolated Delta3D-ligand report

**Decision: `NO-GO — no reproducible incremental free-ligand 3D information`.**
Stop the free-conformer pivot under the preregistered termination rule. Do not
train Davis receptor-additive, ligand-pocket interaction, trainable E(3), or full
cross-attention models.

## Admission and execution

| Item | Decision | Evidence |
|---|---|---|
| Davis ligand identity | **PASS / FINAL** | 69/69 owner-authorized identities; exact SMILES/InChIKey round trips; no additional manual/external sign-off required |
| Label semantics | **PASS** | 6,581 exact Kd cells; 16,741 blank cells treated as right-censored at pKd <= 5, never substituted as exact values |
| Ligand split | **PASS** | 66 frozen scaffold/similarity components; ten outer folds of 6–7 ligands; zero component crossings |
| Feature provenance | **PASS** | gMolAI and Uni-Mol source/checkpoint hashes verified; 69/69 ligands across five conformer seeds; no receptor, bound-ligand, docking or contact inputs |
| Result integrity | **PASS** | 51 models, 510 converged outer fits, all likelihoods independently reconstructed; all protocol/artifact hashes match |
| `delta_3d_ligand` | **FAIL / NO-GO** | seven of eight individual checks failed; the joint progression rule failed |
| Davis `delta_pocket_additive` | **BLOCKED, unchanged** | frozen structural graph has a dominant 323/338 receptor component |
| Davis `delta_3d_x_pocket` | **BLOCKED, unchanged** | same structural identifiability failure; no leakage rule was relaxed |

## Frozen experiment

The test is ligand-only. Each model emits one value per ligand and repeats it over
338 targets, so the estimand is a held-out ligand's marginal pan-Davis potency and
promiscuity—not target preference or ligand-pocket complementarity. Outer and inner
folds keep the frozen Bemis-Murcko/Morgan2 components intact. Training-fold-only
PCA projects gMolAI and each candidate/control branch to 32 dimensions. A linear
Gaussian Tobit model uses nested L2 selection; uncertainty is a paired 10,000-
replicate bootstrap over ligand components.

Five complete ETKDGv3/MMFF94s generations used seeds 20260821–20260825. The 201-
feature deterministic branch and frozen Uni-Mol v1 all-H branch were evaluated as
actual ensembles, coordinate destruction, topology-derived fake 3D, single minimum-
energy conformers and energy permutations. Negative branches have the same added
projection capacity as actual 3D.

## Results

Positive gain means lower NLL than the named reference.

| Model/contrast | NLL or gain | 95% component-bootstrap interval | Exact RMSE | Mean within-target Spearman |
|---|---:|---:|---:|---:|
| gMolAI 2D baseline | 0.96129 | — | 2.6116 | -0.0012 |
| gMolAI + deterministic actual 3D, five-seed mean | 0.96386 | — | 2.7109 | 0.0469 |
| deterministic actual vs gMolAI | **-0.00257** | [-0.01780, 0.01359] | — | — |
| deterministic actual vs destroyed 3D | **-0.00201** | [-0.00566, 0.00098] | — | — |
| deterministic actual vs topology-fake 3D | 0.02155 | [-0.00797, 0.07412] | — | — |
| gMolAI + Uni-Mol actual 3D, five-seed mean | 0.97486 | — | 2.7435 | -0.1494 |
| Uni-Mol actual vs gMolAI | **-0.01356** | [-0.03111, 0.00410] | — | — |
| Uni-Mol actual vs destroyed 3D | **-0.01536** | [-0.03180, -0.00006] | — | — |
| Uni-Mol actual vs topology-fake 3D | 0.01704 | [-0.01237, 0.06168] | — | — |

Only one of five deterministic seeds improved NLL over gMolAI; none of five
Uni-Mol seeds did. Deterministic exact RMSE worsened despite a small positive mean
ranking change. Uni-Mol worsened NLL, exact error and ranking, and its destroyed-
coordinate control outperformed actual coordinates with a component interval below
zero. Single/ensemble and energy-permutation contrasts also included zero. High NDCG
(approximately 0.94) is not treated as evidence because it coexists with near-zero
or negative Spearman and is insensitive on this all-positive exact-value subset.

## Leakage and reproducibility diagnostics

- All 69 ligands stay in exactly one of 66 leakage components; outer fold counts
  are 7/7/7/7/7/7/7/7/7/6 and crossings are zero.
- All predictions and likelihoods are finite. Fitted sigma is 1.594–2.144 pKd, far
  from the broad numerical bounds.
- Independent reconstruction differs from reported NLL by at most 4.5e-16 and from
  stored per-ligand loss by at most 8.0e-13.
- gMolAI exact pretraining-entity exposure is unknown and likely for public
  compounds, but labelled affinity pretraining was not detected. This limits claims
  about a fully cold representation baseline; it does not preferentially rescue
  any failed 3D control.
- Uni-Mol same-batch inference is exactly repeatable. Its pre-outcome GPU audit
  recorded <=0.0047 cross-batch/rigid-transform numerical variation and froze batch
  order/size; this is far smaller than the observed control contrasts.

## Interpretation and recommendation

This result does **not** prove molecular 3D is universally uninformative. It is a
69-ligand/66-component, target-blind pan-kinase marginal estimand; no defensible
repeat-derived practical-equivalence margin is available. It does show that neither
interpretable conformer geometry nor a strong frozen 3D encoder provides robust
incremental prediction over gMolAI here, and geometry-destructive controls do not
support a geometry-specific explanation.

The scientific recommendation is therefore no-go for the current free-conformer
pivot. The separately written kinase-selectivity interaction design is archived as
design-only. A future project would require explicit new authorization, fresh
outcome-sealed data and prospective family/scaffold admission; it must not be
presented as continuation justified by this failed gate.

Machine-readable evidence is in
[`delta3d-ligand-results-v1.json`](evidence/delta3d-ligand-results-v1.json),
[`delta3d-result-validation-v1.json`](evidence/delta3d-result-validation-v1.json)
and [`delta3d-feature-manifest-v1.json`](evidence/delta3d-feature-manifest-v1.json).
