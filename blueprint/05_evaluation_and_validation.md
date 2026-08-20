# Evaluation and Validation Plan

## 1. Evaluation objective

The evaluation must answer a causal design question as closely as retrospective data allow:

> Does receptor-pocket geometry add useful, reproducible information beyond ligand chemistry, target identity, sequence, assay, and dataset priors when the query ligand has no coordinates?

A high correlation on a conventional PDBbind split does not answer this question. The protocol therefore treats leakage controls, matched baselines, conformation stress, external targets, and uncertainty as part of the model—not post hoc analysis.

## 2. Predeclared hypotheses

### Primary hypotheses

**H1 — pocket value:** On the locked dual-novel test, iScore3.0-v1 improves target-macro affinity/ranking performance over the strongest ligand-only model using the same labels and split.

**H2 — geometry value:** iScore3.0-v1 improves over matched models that replace receptor 3D with sequence/identity or bag-of-pocket-residue features.

**H3 — translational geometry:** The improvement in H1/H2 remains detectable on matched apo or predicted receptor structures and is not confined to ligand-conditioned holo structures.

### Secondary hypotheses

**H4 — atomic resolution:** A local receptor heavy-atom branch improves over a residue-only equivariant graph.

**H5 — interaction fusion:** Cross-attention improves over additive, concatenation, and low-rank bilinear fusion.

**H6 — gMolAI value:** Frozen gMolAI atom/global features improve sample efficiency or external generalization over a task-trained ligand GNN and ECFP baseline.

**H7 — endpoint modelling:** endpoint-specific censor-aware heads improve likelihood/calibration without degrading pKd/pKi ranking.

**H8 — uncertainty:** calibrated intervals achieve predeclared coverage on held-out target components, and error increases with predicted uncertainty/out-of-domain score.

### Exploratory hypotheses

Pair-state, pharmacophore-field, surface, top-k-pocket marginalization, large-scale pocket pretraining, and ligand fine-tuning are exploratory until v1 is frozen.

## 3. Experimental governance

Before full-model training:

1. assign immutable observation, ligand, target, structure, and site identifiers;
2. freeze raw and processed manifest hashes;
3. generate splits once from a versioned configuration;
4. publish or escrow test identifiers and nearest-training similarities;
5. hide final external labels from architecture developers where practical;
6. define primary endpoint, metric, baseline, and minimal important effect;
7. define seed count and compute budget;
8. register exclusions and failure handling; and
9. freeze model-selection rules.

A result cannot be moved from exploratory to confirmatory after its test performance is seen. New data or corrected labels create a new benchmark version.

## 4. Split suite

All models, including literature baselines, use identical processed observations.

| Split | Purpose | Headline? | Main threat addressed |
|---|---|---:|---|
| random pair | pipeline sanity and historical comparability | no | none; expected optimistic |
| ligand scaffold/component | novel chemistry | secondary | analogue memorization |
| protein/pocket component | novel targets/sites | secondary | target memorization |
| dual-novel union component | novel ligand and pocket simultaneously | primary | cross-modal transitive leakage |
| temporal | future measurements/structures | primary if dates adequate | retrospective curation leakage |
| LP-PDBbind official | published hard split comparison | secondary | recognized protein/ligand overlap |
| CleanSplit official | PDBbind redundancy benchmark | secondary | structural/chemical redundancy |
| holo→apo/predicted matched | conformation robustness | primary stress test | induced-fit/experimental structure bias |
| target-series holdout | medicinal-chemistry ranking | primary secondary-endpoint | global-label priors |
| CASP16/BayesBind/approved target series | external or blinded screening/ranking | external confirmation | dataset-specific optimization |

The dual-novel split is based on connected components of the union of ligand and protein/pocket similarity graphs. A second audit recomputes similarity using independent methods to catch split-implementation blind spots.

## 5. Baseline suite

### 5.1 Non-neural and nearest-neighbour

- global mean and target mean where target-known use is legitimate;
- ligand nearest neighbour with ECFP and gMolAI similarity;
- pocket nearest neighbour with sequence and structural similarity;
- paired ligand/pocket k-nearest-neighbour regression/classification;
- ridge/elastic net on ECFP and pocket descriptors;
- random forest or gradient boosting on concatenated transparent features; and
- simple physicochemical/assay covariate models.

Nearest-neighbour baselines must retrieve only training examples and report the retrieved identities/similarities.

### 5.2 Ligand-only

- ECFP shallow learner;
- gMolAI global embedding with a matched prediction head;
- gMolAI atom/global model with no protein input;
- task-trained GINE with parameter count matched where practical; and
- optional SMILES-transformer baseline.

### 5.3 Protein-only

- target/family/sequence identity model;
- pooled protein-language-model representation;
- pocket composition/global descriptors;
- residue graph with coordinates;
- residue+atom pocket encoder; and
- detector score/structure-quality metadata only.

A protein-only model can exploit target-specific label distributions. Its performance is evidence of dataset structure, not ligand recognition.

### 5.4 Paired baselines

- DeepDTA/GraphDTA or PSICHIC-style sequence+ligand;
- 3DProtDTA-style late concatenation;
- BANANA reproduction;
- additive pocket+ligand predictions;
- low-rank bilinear fusion; and
- CASTER-DTA-equivalent clean-room model if direct code cannot be used.

Literature-reported numbers are context, never a substitute for retraining on the same split.

## 6. Ablation matrix

Run a staged matrix rather than every combination.

| Question | Reference | Ablation/alternative | Required split |
|---|---|---|---|
| does receptor input help? | full v1 | ligand-only | dual-novel, temporal, external |
| does geometry help? | residue GVP | remove coordinates; sequence/contact proxy | pocket-novel, apo/predicted |
| do receptor atoms help? | residue+atom | residue only | dual-novel, chemistry strata |
| does gMolAI help? | frozen gMolAI | ECFP; task GINE | learning curves, dual-novel |
| does atom-level ligand fusion help? | atom cross-attention | global-only fusion | dual-novel |
| does cross-attention help? | cross-attention | concatenate/bilinear | all core splits |
| does pocket focus help? | supplied pocket | whole protein; shuffled remote pocket | target-novel |
| does site detection work? | supplied site | P2Rank/fpocket top-1/top-k | external receptor-only |
| does endpoint separation help? | multi-head | pooled label; single endpoint | endpoint-stratified |
| does censor handling help? | interval NLL | discard; boundary-as-point | censored subset |
| is holo signal privileged? | holo | matched apo/predicted | conformation set |
| is pretraining useful? | pretrained pocket/gMolAI | random/task-only | sample-size curves |
| is uncertainty useful? | ensemble+calibration | single uncalibrated | all locked tests |

Only ablations relevant to an accepted module are required; do not use external tests to select among them.

## 7. Negative and invariance controls

### 7.1 Information-boundary tests

- fail preprocessing if any ligand conformer or coordinate property exists;
- inspect serialized ligand artifacts for coordinate arrays and 3D descriptors;
- block imports/calls to conformer-generation functions in strict pipelines;
- trace feature lineage and verify every cross-modal feature is derived after independent encoding;
- scan configuration and dependency calls for pose/contact/distance labels; and
- use synthetic ligands whose atom ordering changes while topology does not.

### 7.2 Model negative controls

- permute ligands across pockets within endpoint/source strata;
- permute pockets across ligands while preserving target-frequency strata;
- replace the supplied pocket with a remote non-pocket region of comparable size;
- zero or randomize receptor coordinates while preserving residue identities;
- use only pocket size, detector rank, PDB resolution, and publication year;
- permute labels within target/assay blocks;
- train with random split and compare the inflation relative to union components; and
- remove close training neighbours at progressively stricter thresholds.

A pocket permutation that leaves performance unchanged is evidence the fusion pathway is being ignored.

### 7.3 Geometric invariances

For the same pocket:

- arbitrary global rotation and translation must leave predictions unchanged;
- atom/residue input permutation must leave predictions unchanged;
- reflections should follow the declared model symmetry and chirality policy;
- small coordinate perturbations should change predictions smoothly; and
- padding/batching should not affect predictions.

Use numerical tolerances fixed by precision mode.

## 8. Metrics

### 8.1 Quantitative affinity

Report per endpoint:

- mean absolute error in pK units;
- RMSE, with sensitivity to outliers;
- Pearson correlation for legacy comparability;
- Spearman rho and Kendall tau;
- concordance index, including a censor-aware form where applicable;
- negative log likelihood or interval log score;
- calibration slope/intercept;
- empirical coverage and width for 50%, 80%, 90%, and 95% intervals; and
- fraction outside the calibrated applicability domain.

Global correlations can be driven by target-specific mean differences. The primary ranking metric is the macro-average within-target Kendall/Spearman statistic on targets with sufficient comparable ligands. The primary absolute-error metric is target-cluster-macro MAE.

### 8.2 Virtual screening

Report:

- area under the precision–recall curve with class prevalence;
- ROC AUC only as a secondary metric;
- enrichment at fixed screened fractions with confidence intervals;
- BEDROC or an explicitly parameterized early-recognition metric;
- BayesBind EF^B on its benchmark;
- top-k hit rate and number needed to test;
- Brier score, log loss, reliability diagram, and expected calibration error; and
- ligand-only versus full-model delta.

For constructed decoys, report results separately and disclose decoy generation.

### 8.3 Ranking series

For CASP16, D3R/SAMPL, and medicinal-chemistry series:

- Kendall tau and Spearman rho per system;
- pairwise ordering accuracy with experimental uncertainty;
- MAE after any predeclared target-specific offset, reported alongside absolute MAE;
- top-compound recovery;
- uncertainty-aware rank probability; and
- macro-average across systems.

Do not pool all pairs and allow large series to dominate.

### 8.4 Operational metrics

- preprocessing success and warning rates by stratum;
- time per receptor/site and per ligand;
- cached compounds per second;
- peak CPU/GPU memory and storage;
- batch-size scaling;
- deterministic-repeat discrepancy; and
- energy/cost where available.

## 9. Statistical analysis

### 9.1 Resampling unit

Use paired bootstrap over independent target/pocket components or whole target series, not individual protein–ligand rows. For a small number of external systems, show every system and use exact/permutation or hierarchical estimates rather than a misleading narrow interval.

### 9.2 Paired comparisons

Compare models on the same resampled clusters and report:

- point estimate;
- 95% confidence or compatible interval;
- probability/direction of improvement where a Bayesian analysis is used;
- standardized and domain-scale effect;
- win/tie/loss across target clusters; and
- compute/failure-rate trade-off.

A statistically nonzero but negligible improvement is insufficient. The committee should freeze a minimal important difference before test labels are opened. Provisional values for discussion are 0.10 pK MAE and 0.03 macro rank correlation, but these are not final until justified against experimental noise and project utility.

### 9.3 Multiple comparisons

Use a hierarchical testing order:

1. H1 pocket value;
2. H2 geometry value;
3. H3 apo/predicted persistence;
4. H4–H8 secondary hypotheses; and
5. exploratory extensions.

Control the family-wise error or false-discovery rate within the secondary family. Report all attempted models, not only the winner.

### 9.4 Seeds and variability

Use at least five independent seeds for final small/medium models unless a compute-power analysis supports another number. Separate variability from initialization, data resampling, receptor conformation, and pocket detection. A single favourable seed cannot define the final checkpoint.

## 10. Leakage audit

For every test prediction, attach:

- exact-duplicate indicators;
- maximum ECFP and gMolAI similarity to training ligands;
- scaffold overlap;
- maximum protein sequence/domain identity;
- maximum pocket structural similarity;
- nearest joint train neighbour;
- same publication/assay-series indicator;
- structure and measurement date gaps; and
- known pretraining exposure.

Plot performance and full-minus-ligand-only improvement against these variables. Repeat primary analyses after removing the closest similarity bins. If the benefit disappears in the genuinely novel region, the scope of the claim must shrink.

Use independent software/feature implementations for at least one ligand and one protein similarity audit. Split-generation features can contain bugs.

## 11. Conformation and site evaluation

### 11.1 Matched conformation test

For each matched system, score the same ligand list against aligned holo, apo, and predicted pockets. Do not reselect the best structure per ligand. Report:

- prediction/rank change by conformation;
- pocket RMSD and local side-chain displacement;
- confidence/missingness;
- full-versus-ligand-only delta; and
- whether uncertainty increases when geometry degrades.

### 11.2 Pocket detection

Evaluate in two stages:

1. detector recall/overlap of the biological site at top-1/top-3/top-k; and
2. affinity/screening performance conditional on supplied, correct detected, and incorrect detected pockets.

For top-k aggregation predeclare max, mean, learned mixture, or uncertainty-weighted marginalization. A max over many pockets can inflate scores and false positives.

### 11.3 Robustness perturbations

Apply receptor-only perturbations calibrated to real preprocessing variability:

- pocket radii/context shells;
- alternative side-chain repairs;
- protonation feature variants;
- removal/retention of water, metals, or cofactors;
- coordinate noise;
- low-confidence residue masking;
- chain/assembly alternatives; and
- detector/version changes.

Report both prediction stability and accuracy. Stability around a consistently wrong prediction is not success.

## 12. Calibration and abstention

Fit calibration only on held-out validation target components. Evaluate calibration marginally and by endpoint, affinity range, target family, chemical class, structure source, and applicability bin.

Predeclare abstention policies such as:

- unsupported chemistry/structure: mandatory failure;
- joint-domain score below threshold: abstain;
- prediction interval wider than threshold: low-confidence flag;
- large conformation/pocket disagreement: request expert review.

An acceptable applicability score should order errors: retained high-confidence subsets improve in accuracy as coverage decreases. Plot risk–coverage curves and guard against subgroup exclusion.

## 13. External evaluation protocol

1. freeze source, preprocessing, models, ensembles, and calibrators;
2. hash all artifacts and environment;
3. acquire or unlock external labels only after freeze;
4. run one scripted pass with no manual record deletion;
5. publish all failures and warnings;
6. calculate predeclared metrics and target-cluster intervals;
7. permit only clearly labelled post hoc diagnosis; and
8. never overwrite the confirmatory prediction file.

CASP16, BayesBind, and any approved target series serve different endpoints; no single aggregate “external score” is created. LIT-PCBA is included only as a red-team diagnostic or after an independently audited rebuild, not as unmodified confirmatory evidence.

## 14. Acceptance criteria

### Data and implementation

- all strict information-boundary and invariance tests pass;
- zero train/test component leakage at frozen thresholds;
- complete provenance/licence fields for used observations;
- less than a committee-approved failure rate on in-domain inputs;
- full reproducibility from immutable manifests; and
- no external-label access during model selection.

### Scientific

- H1 improvement over the best ligand-only comparator exceeds the frozen minimal effect and its target-cluster interval excludes no meaningful benefit;
- H2 shows that observed 3D geometry, not merely sequence/composition, contributes;
- H3 shows a useful effect on apo or predicted structures;
- at least one locked external benchmark improves over non-structural and nearest-neighbour baselines;
- uncertainty is calibrated within predeclared tolerances; and
- gains persist after close-neighbour removal.

A failure of H1–H3 is a useful negative scientific result and triggers scope redesign; it must not be hidden by switching to the random split.

## 15. Reporting package

Every reported table should be generated from immutable prediction files and include:

- data/split/model/configuration hashes;
- number of targets, pockets, ligands, observations, and failures;
- endpoint and mapping-tier composition;
- model parameter count and training/inference compute;
- point estimates and target-cluster intervals;
- nearest-neighbour leakage strata;
- seed-level results;
- calibration and applicability; and
- exact baseline parity conditions.

The model card must state intended use, exclusions, data rights, receptor/ligand information boundary, uncertainty limitations, and known failure modes.
