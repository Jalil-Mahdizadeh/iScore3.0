# iScore3.0 Gate-0/1 Scientific Feasibility Report

**Decision date:** 20 August 2026<br>
**Scope:** bounded literature, adapter, provenance, pocket, leakage, and shallow-control study; no full iScore3.0 architecture was trained.

## Executive decision

**NO-GO for full architecture training at this gate.** Continue only a narrow remediation phase for data provenance, independent component count, validated pocket-structure similarity, and stronger protein encoders/baselines. Re-present the project to the committee before cross-attention or other high-capacity architecture work.

This is not a conclusion that receptor geometry is useless. It is the narrower finding that the current evidence does not justify the proposed full model:

1. protein-3D plus ligand-2D/no-pose affinity prediction is established prior art, including PSG-BAR, CASTER-DTA, 3DProtDTA, AttentionMGT-DTA, PLMCA, HoloProt, and Graph_RG;
2. the strict pilot contains only 61 measurements in five independent construct components;
3. the predeclared S1 pocket+gMolAI control improved pooled RMSE by only 0.040 pKd relative to gMolAI-only, versus the required 0.10 pKd, with a 95% component-bootstrap interval of `[-0.152, +0.080]` for candidate minus comparator; and
4. a permuted-pocket negative control performed slightly better than the real S0 pocket fusion, with an interval that includes no difference.

The most valuable Gate-0/1 outcomes are therefore the detected qualifier/provenance failure, a tested gMolAI atom adapter, a rebuildable strict pilot, and an honest negative baseline result.

## Gate ledger

| Question | Result | Evidence-based interpretation |
|---|---|---|
| Is the broad modality novel? | **FAIL** | PSG-BAR, PLMCA, AttentionMGT-DTA and earlier direct methods already combine structural protein features with a graph/SMILES ligand. |
| Is a narrower contribution still plausible? | **CONDITIONAL PASS** | Strict end-to-end ligand-3D exclusion, provenance, leakage control, conformation uncertainty, and demonstrated incremental geometry could form a contribution. None is yet a performance claim. |
| Is gMolAI `node_z` technically usable? | **PASS** | Dimensions, atom bijections, preprocessing, hashes, repeatability, equivalent-SMILES behaviour, order independence, and label quarantine passed. |
| Is gMolAI provenance fully resolved? | **FAIL / OPEN** | No root licence was detected and the 223,180,699-graph pretraining corpus lacks an identity ledger, so exact pilot-compound exposure cannot be determined. |
| Is the original 80-row RCSB candidate table trustworthy as labelled? | **FAIL** | BindingDB 202608 showed that 12 apparent exact Kd values were censored lower bounds whose `>` qualifiers were absent from the RCSB query output. |
| Is the corrected strict pilot fit for bounded diagnostics? | **PASS** | 61 exact, uncensored Kd measurements have exact PDB, ligand, UniProt, and measurement-publication provenance; five reference labels remain quarantined. |
| Is the pilot sufficient for external/general performance claims? | **FAIL** | Five independent construct components and no held-out external benchmark are insufficient. |
| Did the required essential controls run? | **PASS** | Global/target means, nuisance-only, ligand-only ECFP/gMolAI, sequence-only, pocket-only, simple fusion, KNN, two-sided KNN, and permuted-pocket controls ran on identical splits. |
| Did receptor-pocket information pass the progression criterion? | **FAIL** | The predeclared S1 gMolAI+pocket improvement was below 0.10 pKd and its confidence interval crossed zero. |
| Is progression to the full iScore3.0 architecture scientifically justified now? | **NO** | Data and geometry-signal risks remain larger than the evidence for added architectural capacity. |

## Competitive and novelty update

The review now contains 82 publication/resource records and 47 software audits. Three additions are decisive:

- **PSG-BAR (2022)** is a direct full-protein-3D plus 2D-ligand affinity model. Its cold-both degradation also reinforces the need for joint protein/ligand novelty evaluation.
- **AttentionMGT-DTA (2024)** uses structural pocket graphs, ESM features, and a 2D ligand graph. Its released script selects and schedules on the test set, so reported results cannot substitute for an iScore3.0-controlled reproduction.
- **PLMCA (2026)** combines PLMs, protein geometry/physicochemistry, ligand graphs, cross-attention, and assay conditions. No official source was found, leaving its complete feature boundary unverified but its architectural prior-art relevance clear.

BlendNet is complex-free at inference but uses PLIP complex-interaction labels for teacher training. MMPD-DTA, AlignNet, and LigUnity are non-strict because released pathways use ligand or complex coordinates. AttentionSiteDTI, BindingSite-AugmentedDTA, PGraphDTA, and CSCo-DTA are important binary/site/contact-map precedents.

**Permissible future novelty wording:** “A leakage-audited and provenance-first receptor-structure/SMILES affinity system that forbids ligand 3D information throughout the declared strict track and quantifies incremental pocket-geometry value under conformation-stressed, similarity-component, and external evaluation.” The modality combination, cross-attention, or “docking-free” label is not novel by itself.

## gMolAI-v2.0 adapter audit

The adapter pins upstream revision `29ced8352886692c00a792d24a183f05e5de0059`, checkpoint SHA-256 `02f49a2a94ddfc9dc780cc3d5f1a3df54306ae0fdc5d4b3767e3fd2e7f27b05e`, calibrator/config hashes, source-file hashes, and the container digest. The audited contract is:

| Object | Dimension |
|---|---:|
| atom input | 48 |
| bond/edge input | 15 |
| atom `node_z` | 128 |
| raw graph `graph_z` | 256 |
| released molecule vector | 384 |

On 80 unique candidate ligands comprising 1,971 atoms:

- canonical-to-input and input-to-canonical atom maps were bijective and aligned with `node_z` rows;
- atom symbols, graph edges, features, and output hashes were recorded per molecule;
- repeat runs, reversed batch order, and equivalent SMILES were exactly reproducible on CPU;
- CPU/GPU maximum absolute difference was at most `9.54e-06`, within the frozen `2e-4` tolerance;
- atom maps, isotopes, input conformers, and CXSMILES extensions are rejected at the interface; and
- historical site-reference rows have blank labels and were not encoded.

Two caveats remain release-blocking: no explicit upstream licence was found, and exact pretraining identity overlap is unknowable without the withheld training-entity ledger. The latter is representation exposure, not evidence of affinity-label leakage, but all gMolAI results must remain paired with non-pretrained ECFP controls.

## Strict S0/S1 pilot

The canonical pilot was created by scanning 3,234,499 BindingDB 202608 rows and reconciling every RCSB candidate against exact PDB, ligand, UniProt, Kd relation/value, and measurement-publication fields.

| Property | Strict-v2 value |
|---|---:|
| exact uncensored supervised Kd observations | 61 |
| label-quarantined historical site references | 5 |
| exact construct / union components | 5 |
| supervised component sizes | 19, 15, 10, 9, 8 |
| receptor feature rows | 122 (61 S0 + 61 S1) |
| pKd range | 3.6253–8.3010 |
| mean pKd | 5.2361 |

The retained constructs are nitric oxide synthase oxygenase, Peregrin, BRD4, p53, and mycocyclosin synthase. Twelve of 80 candidate labels failed exact-uncensored reconciliation. The remaining WDR5 group then had only seven measurements and was dropped under the predeclared minimum of eight.

S0 excludes query-ligand coordinates from features but uses the query-holo receptor conformation. S1 uses one frozen, earlier exact-sequence reference receptor and a transferred 6-Angstrom residue set for all ligands in a construct group. Both therefore retain disclosed holo privilege; S1 removes query-specific conformation but is not an apo test.

The BindingDB audit’s top-level status is intentionally `FAIL` because the candidate dataset failed. Its nested strict-dataset status is `PASS`; the corrected strict-v2 table is canonical. The v1 candidate and result files remain only as an audit trail.

## Leakage diagnostics

Union components include exact construct, exact Bemis–Murcko scaffold, ECFP4 Tanimoto at least 0.35, global sequence identity at least 0.30, local-site identity at least 0.50, shared structure publication, and shared measurement publication. Leave-one-component-out predictions have zero component overlap by construction.

No similarity or provenance edge joined two retained constructs. The largest cross-construct similarities were ECFP Tanimoto 0.272, full-sequence identity 0.295, and local-site identity 0.333, all below their thresholds. A validated Foldseek/TM-align pocket-structure edge was not available; this is an explicit limitation, not a passed check.

Other unresolved channels are gMolAI pretraining identity exposure, S0 induced-fit information, and lack of independent source-paper table extraction. Random pair cross-validation is reported only as a shortcut diagnostic.

## Essential baseline results

Primary results use union-component leave-one-out. All 25 controls use the same 61 observations; preprocessing and hyperparameter selection occur within training data. RMSE/MAE are in pKd units.

| Model | View | RMSE | MAE | Pearson | Target-macro Spearman |
|---|---|---:|---:|---:|---:|
| ECFP + sequence ridge | sequence/ligand | **1.227** | 1.040 | -0.212 | 0.135 |
| ECFP ridge | ligand only | 1.240 | 1.060 | -0.250 | 0.209 |
| global training mean | none | 1.263 | 1.054 | -0.681 | n/a |
| nuisance ridge | metadata only | 1.267 | 1.059 | -0.620 | -0.025 |
| ECFP + pocket ridge | S0 | 1.269 | 1.108 | -0.317 | 0.179 |
| sequence ridge | protein only | 1.275 | 1.057 | -0.582 | n/a |
| ECFP + pocket ridge | S1 | 1.298 | 1.138 | -0.357 | 0.167 |
| ECFP KNN | ligand only | 1.400 | 1.184 | -0.085 | 0.136 |
| pocket-structure ridge | S0 | 1.405 | 1.241 | -0.401 | -0.035 |
| pocket-structure ridge | S1 | 1.441 | 1.270 | -0.649 | n/a |
| permuted pocket + gMolAI ridge | S0 negative control | 1.552 | 1.366 | -0.316 | 0.207 |
| sequence-identity KNN | sequence only | 1.556 | 1.266 | -0.451 | n/a |
| gMolAI + pocket ridge | S1 | 1.567 | 1.374 | -0.294 | 0.196 |
| gMolAI + pocket ridge | S0 | 1.587 | 1.417 | -0.340 | 0.246 |
| gMolAI ridge | ligand only | 1.607 | 1.396 | -0.232 | 0.156 |

All pooled leave-component-out correlations are weak or negative. This is expected when transferring among only five chemically and biologically different series and reinforces that pooled random-pair results are misleading. For example, ECFP+sequence RMSE is 0.704 on the non-headline random split but 1.227 on the primary split.

Paired component bootstrap results (`candidate RMSE - comparator RMSE`, negative favours candidate):

| Candidate vs comparator | Point delta | 95% interval | Interpretation |
|---|---:|---:|---|
| gMolAI+pocket S1 vs gMolAI | -0.040 | [-0.152, +0.080] | below required -0.10; includes no improvement |
| gMolAI+pocket S0 vs gMolAI | -0.021 | [-0.099, +0.131] | inconclusive |
| ECFP+pocket S0 vs ECFP | +0.028 | [-0.066, +0.143] | inconclusive, point estimate worse |
| ECFP+pocket S1 vs ECFP | +0.058 | [-0.037, +0.192] | inconclusive, point estimate worse |
| permuted-pocket gMolAI vs real-pocket S0 | -0.035 | [-0.084, +0.055] | negative control is not distinguishable and is slightly better at the point estimate |

The shallow 52-dimensional pocket descriptor is a diagnostic, not the planned geometric neural encoder. Its failure does not test the full architecture; it tests whether there is enough robust pocket signal to justify that architecture now. The answer is no.

## Principal risks

1. **Data mapping:** the qualifier-loss incident proves that aggregator fields cannot be trusted without source reconciliation.
2. **Statistical independence:** five components are too few for stable cross-target inference.
3. **Structural leakage:** no validated local pocket-structure similarity edge was available.
4. **Holo privilege:** S0 can encode query-induced conformation; S1 remains historical holo rather than apo/predicted.
5. **Encoder provenance:** gMolAI identity overlap and redistribution rights remain unresolved.
6. **Baseline scope:** PSG-BAR/BANANA/AttentionMGT/modern sequence baselines were audited but not fully reproduced on this tiny pilot; the implemented controls are deliberately shallow.
7. **Signal risk:** current pocket features fail to beat ligand/sequence controls and the permuted-pocket test.

## Recommendation and bounded next gate

Do not train cross-attention, residue GVP, atom-level protein graphs, or a large multi-task model yet. Authorize a second bounded feasibility tranche only if the committee accepts the following deliverables:

1. grow to at least 20–30 independent, high-confidence construct/site components without weakening exact measurement provenance;
2. independently re-extract a stratified subset of values from source publications and add qualifier/unit regression tests at acquisition;
3. add validated Foldseek/TM-align and local pocket-alignment union edges;
4. add apo/predicted receptor views and report conformation sensitivity;
5. resolve gMolAI licensing and obtain or reconstruct a pretraining identity ledger, or replace it for confirmatory work;
6. reproduce one direct strict comparator (BANANA for activity or a clean-room PSG-BAR/3DProtDTA-style affinity model) and one modern sequence model on the same components; and
7. repeat the predeclared geometry-value test with a frozen invariant pocket encoder.

Progression to full iScore3.0 is justified only if S1/apo structural fusion improves the matched ligand+sequence baseline by the committee-approved practical margin with component-level uncertainty excluding no benefit, while the permuted-pocket control does not.

## Canonical evidence

- `data/manifests/gate01-pilot-strict-v2.json`
- `reports/gate01/evidence/gmolai-adapter-audit-v1.json`
- `reports/gate01/evidence/bindingdb-provenance-audit-v2.json`
- `reports/gate01/evidence/leakage-diagnostics-strict-v2.json`
- `reports/gate01/evidence/baseline-metrics-strict-v2.json`
- `data/manifests/gate01-baselines-strict-v2.json`
- `third_party/reviewed_sources.tsv`

The full reproduction commands and immutable hashes are in [reproducibility.md](reproducibility.md).

Final verification passed 18/18 project tests, 67/67 tests at the pinned upstream gMolAI revision, Python byte-compilation, JSON/TSV schema checks, manifest hash checks, and label-boundary assertions.
