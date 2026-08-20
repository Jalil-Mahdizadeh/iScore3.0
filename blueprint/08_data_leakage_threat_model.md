# Data Leakage Threat Model and Audit Protocol

**Priority:** release-blocking  
**Evidence cut-off:** 20 August 2026

## 1. Position

Data leakage is the central validity risk for iScore3.0. It is broader than exact duplicate rows and cannot be solved by one protein sequence cutoff, one Bemis–Murcko split, or one dataset labelled “leak-proof”.

The project will not report one binary “leak-free” flag. It will:

1. define the intended deployment regime;
2. enumerate every information path from future/test data into training or model selection;
3. prevent literal contamination;
4. quantify continuous chemical, protein, pocket, assay, and provenance similarity;
5. evaluate several novelty regimes;
6. run models designed to expose shortcuts; and
7. make the complete per-test nearest-neighbour audit available with results.

LP-PDBbind and PDBbind CleanSplit are required standardized protocols, but neither is treated as a universal proof of no leakage.

## 2. Why this problem is unusually severe

Protein–ligand datasets are bipartite and highly redundant:

- one ligand appears against many targets;
- one target appears with congeneric ligand series;
- homologous proteins can share binding profiles;
- globally dissimilar proteins can have similar pockets;
- one biological target has many structures, constructs, mutants, and conformations;
- one experiment is copied into multiple databases;
- one publication/patent contributes a coherent medicinal-chemistry series;
- PDBbind releases and CASF benchmarks are historically intertwined;
- labels from Kd, Ki, IC50, and pChEMBL may describe related but nonidentical quantities; and
- huge pretrained encoders may already have seen evaluation entities.

A random pair split therefore allows a model to reconstruct a test result from either side of the interaction, or from source/assay artifacts, without learning general pocket–ligand compatibility.

## 3. Important distinction: leakage versus legitimate similarity

Similarity is not automatically misconduct or even undesirable. In lead optimization, training on earlier analogues for the same target is the intended use. In target deorphanization, the ligand may be known but the protein novel. In scaffold hopping, the pocket is familiar but the chemistry must be novel. In broad virtual screening, both may be unfamiliar.

The methodological error is to claim one regime while evaluating another, or to allow information that would not exist at the declared prediction time.

iScore3.0 will report four deployment regimes:

| Regime | Protein/pocket novelty | Ligand novelty | Example use | Appropriate split |
|---|---|---|---|---|
| I: series interpolation | known | close analogues allowed | lead optimization | target/assay time split with series chronology |
| II: scaffold hopping | known | novel scaffold/low similarity | new chemotype for known site | ligand novelty tiers within held-out series |
| III: target transfer | novel pocket/family | known or related chemistry | new target/site | target/pocket cold |
| IV: dual novelty | novel | novel | broad generalization | union-component dual cold plus temporal/external |

Random pair results are an optimistic diagnostic and carry no deployment claim.

## 4. Evidence from the debate

### 4.1 Chemical memorization and decoy bias

[Wallach and Heifets, 2018](https://doi.org/10.1021/acs.jcim.7b00403) showed that training–validation redundancy in ligand classification benchmarks correlates strongly with performance and introduced AVE bias analysis. [Sieg, Flachsenberg, and Rarey, 2019](https://doi.org/10.1021/acs.jcim.8b00712) demonstrated that structure-based virtual-screening benchmarks can contain learnable chemical biases. [Chen et al., 2019](https://doi.org/10.1371/journal.pone.0220113) attributed strong DUD-E CNN enrichment substantially to analogue and decoy bias rather than molecular-recognition generalization.

Consequence: DUD-E is diagnostic only; every screening benchmark is accompanied by ligand-only and nuisance-feature models.

### 4.2 PDBbind/CASF overlap

[Su et al., 2020](https://doi.org/10.1021/acs.jcim.9b00714) studied how scoring power depends on training-set similarity. [Gusev et al., 2025](https://doi.org/10.1038/s42256-025-01124-5) reported substantial PDBbind-to-CASF contamination, introduced PDBbind CleanSplit using combined structural/chemical filtering, and showed large drops for retrained models after leakage reduction. Their [GEMS repository](https://github.com/camlab-ethz/GEMS) releases split JSON, filtering code, similarity matrices, and a nearest-neighbour search baseline.

Consequence: original PDBbind→CASF performance cannot support iScore3.0’s generalization claim. CleanSplit is run exactly and independently audited.

Independent studies reinforce the need for single-side controls. [Predicting or Pretending, 2020](https://doi.org/10.3389/fphar.2020.00069) found ligand-only and protein-only models could retain much of apparent PDBbind performance, [Volkov et al., 2022](https://doi.org/10.1021/acs.jmedchem.2c00487) found nearest ligand/protein and single-side descriptors dominated tested deep complex models, and [Latent Biases, 2023](https://doi.org/10.1021/acsomega.2c06781) linked PDBbind selection to a shortage of weak/nonbinding examples.

Consequence: full-model accuracy is never credited to learned interaction unless it improves over ligand-only, protein-only, and two-sided retrieval controls on the relevant novelty strata.

### 4.3 LP-PDBbind

[Leak Proof PDBbind, 2026](https://doi.org/10.1021/acs.jpcb.5c08598) reorganizes PDBbind by protein sequence/structure and ligand similarity and evaluates retrained classical and machine-learning scoring functions. The [source metadata](https://github.com/THGLab/LP-PDBBind) exposes official split and cleanliness levels.

Consequence: LP-PDBbind is a required published protocol, but its selected thresholds describe one novelty definition. Results are stratified by the authors’ tiers and by iScore3.0’s independent continuous similarity audit.

### 4.4 The counterargument

[Li et al., 2021](https://doi.org/10.1093/bib/bbab225) argued that machine-learning scoring functions trained on dissimilar complexes can still outperform classical counterparts and that similar training examples are legitimate in prospective settings when they truly exist at prediction time. This is an important correction to a simplistic “any similarity is leakage” position.

Consequence: iScore3.0 reports both deployment-matched interpolation and hard novelty. It does not discard useful historical data from the operational model merely to create an artificially difficult training set; instead, it separates model training from claim-specific evaluation.

### 4.5 Two-dimensional split methodology

[DataSAIL, 2025](https://doi.org/10.1038/s41467-025-58606-8) formalizes similarity-aware splitting for interaction data, where both sides of a bipartite pair can leak. Some interactions cannot be retained while keeping both entity groups isolated.

Consequence: DataSAIL is evaluated against a transparent union-component implementation. Any discarded edge/entity and class/endpoint distortion is reported.

### 4.6 LIT-PCBA audit

The reproducible 2025 preprint [Data Leakage and Redundancy in the LIT-PCBA Benchmark](https://arxiv.org/abs/2507.21404) and its [audit source](https://github.com/sievestack/LIT-PCBA-audit) report exact train–validation duplicates, query compounds present in train/validation, and extensive analogue redundancy. The paper shows that a simple memorization/similarity procedure can exploit those artifacts.

Consequence: because this is preprint evidence, iScore3.0 will independently reproduce the audit on a hash-pinned official archive. Unmodified LIT-PCBA is downgraded from external confirmation to a red-team diagnostic. A confirmatory LIT-PCBA-derived test would require reconstruction from underlying assay records with isolated queries, identities, analogues, and provenance.

### 4.7 Pocket similarity below global homology

Binding pockets can be structurally/chemically similar even when proteins have low global sequence or fold similarity. Pocket-comparison studies such as [APoc](https://pmc.ncbi.nlm.nih.gov/articles/PMC3582269/) demonstrate this explicitly.

Consequence: global sequence identity is never the only protein-side leakage control. Pocket structure, local sequence/chemistry, domain/family, and target binding-profile relationships are audited.

### 4.8 Current 2026 claims

The June 2026 bioRxiv preprint [Identifying and Addressing Systematic Data Leakage in Protein–Ligand Affinity Benchmarks](https://www.biorxiv.org/content/10.64898/2026.06.29.735309v1) argues that sequence-threshold splits can leave “target mirroring” and releases the [Novelty-Tiered Affinity Benchmark](https://doi.org/10.5281/zenodo.19665374), binned by maximum ligand similarity after a time split. [HonestAffinity](https://arxiv.org/abs/2606.03422) reports split-conditioned reversals in the apparent value of protein/pocket priors.

These are timely preprints, not settled peer-reviewed conclusions. They strengthen the case for novelty tiers, ligand-only controls, and target-profile audits; they do not replace peer-reviewed LP-PDBbind, CleanSplit, or prospective experiments.

## 5. Leakage taxonomy

The release audit covers all rows below.

| ID | Leakage/shortcut channel | Example | Prevention/detection |
|---|---|---|---|
| L-01 | exact observation duplication | same Kd imported through PDBbind, BindingDB, and BioLiP2 | source-lineage graph and experiment fingerprint |
| L-02 | ligand identity aliases | salt, tautomer, stereochemical or canonical-SMILES variants cross folds | original plus standardized identity layers |
| L-03 | close chemical analogues | same medicinal-chemistry series split row-wise | fingerprints, scaffold/network, publication/series components |
| L-04 | exact target aliases | UniProt aliases, isoforms, constructs, chain IDs | normalized sequence/construct identity |
| L-05 | protein homology | close homologues split across folds | multi-threshold sequence/domain clusters |
| L-06 | pocket homology | similar local binding sites in globally dissimilar proteins | receptor-only pocket alignment/similarity graph |
| L-07 | target mirroring | different targets have highly correlated activity profiles for shared compounds | cross-target shared-compound correlation audit |
| L-08 | repeated structures | same complex/ligand in alternate PDB entries or biological assemblies | structure/ligand/provenance linkage |
| L-09 | repeated assay/publication | one series or assay plate appears in train and test | assay/document/patent components |
| L-10 | temporal leakage | later measurements, structures, or curation inform earlier prediction | earliest-public-date snapshots and time-aware acquisition |
| L-11 | benchmark overlap | CASF record or close neighbour is in PDBbind training | benchmark-first exclusion plus global overlap audit |
| L-12 | preprocessing leakage | scaler, PCA, imputation, feature selection, calibration fit on all rows | fit artifacts on train only and record parent split |
| L-13 | label-derived features | target mean, threshold, assay rank, or replicate consensus uses test labels | feature lineage and fold-scoped aggregation |
| L-14 | hyperparameter leakage | repeated test/leaderboard use selects architecture | locked labels, evaluation owner, immutable first run |
| L-15 | pretrained entity exposure | gMolAI/ESM/pocket pretraining saw test entities | exposure inventory and nearest-neighbour sensitivity |
| L-16 | pretrained label exposure | supervised pretraining contains the same assay/affinity | global observation-lineage deduplication |
| L-17 | pocket-definition leakage | query ligand coordinates define its own site/crop | independent site manifest; transferred/receptor-only site |
| L-18 | pose/complex leakage | contact/distance labels or ligand 3D enter a nominally pose-free model | strict interface and serialized-artifact inspection |
| L-19 | holo privileged signal | receptor conformation is induced by the exact test ligand | matched apo/predicted evaluation and disclosure |
| L-20 | decoy/source bias | actives and inactives come from distinguishable generators | source classifier, confirmed inactives, PU learning |
| L-21 | metadata/source shortcuts | resolution, year, lab, target count, missingness predict labels | nuisance-only models and source-held-out analysis |
| L-22 | retrieval/network leakage | inference searches a mutable database containing test labels | offline frozen train-only indexes |
| L-23 | split-transitivity failure | A resembles B chemically; B resembles C by pocket; A/C cross folds | connected components of union relation graph |
| L-24 | post-split filtering | exclusions after splitting reconnect or distort components | process before split or regenerate and re-audit |
| L-25 | calibration leakage | test labels used for uncertainty scaling/threshold selection | target-component calibration set only |
| L-26 | human-in-the-loop leakage | manual “bad” test cases removed after prediction | frozen exclusions; publish all failures |

## 6. Canonical identity and provenance graph

Before any split, build a normalized knowledge graph with nodes for:

- raw database record;
- experimental observation;
- assay;
- publication/patent;
- ligand form and parent connectivity;
- target accession, construct, sequence, and domain;
- receptor structure, assembly, chain, conformation;
- biological site and pocket instance; and
- source database/release.

Edges encode “same experiment”, “derived from”, “same parent ligand”, “same construct”, “homologous”, “same site”, “alternate structure”, and “imported by”. All records connected to the same experimental observation receive one split assignment.

The raw record ID is retained even when deduplicated. Conflicts are observations, not aliases.

## 7. Independent similarity views

No learned split is trusted to audit itself.

### 7.1 Ligand similarity

Compute and retain:

- exact standardized stereochemical identity;
- connectivity-only identity;
- Bemis–Murcko scaffold and scaffold-network relations;
- ECFP/Morgan similarities at at least two radii;
- atom-pair/topological-torsion similarity;
- maximum common substructure or series membership for audits;
- gMolAI embedding similarity; and
- publication/assay-series membership.

Thresholds are not selected to maximize model performance. Report continuous nearest similarity and sensitivity under several predeclared thresholds. The 2026 NTAB bins—below 0.35, 0.35–0.5, 0.5–0.7, 0.7–below 1, and exact—are included for comparability, not assumed universally optimal.

### 7.2 Protein and pocket similarity

Compute independently:

- exact construct sequence;
- MMseqs2 global and domain/local sequence identity and coverage;
- Pfam/CATH/SCOP or equivalent family/domain linkage;
- Foldseek/TM-align global structural similarity;
- receptor-only pocket residue alignment;
- pocket shape/chemical similarity using a pinned method;
- local pocket sequence/geometry fingerprints;
- shared biological-site annotation; and
- training-only learned pocket embedding similarity as an additional diagnostic.

A low global sequence identity does not overrule a high local pocket similarity.

### 7.3 Assay and target-profile similarity

Link observations sharing document, patent family, assay identifier, laboratory/source where available, endpoint/context, or compound series. For shared-compound targets, compute binding-profile correlations using training-era measurements only. High target-profile similarity creates a mirroring edge or at minimum a reporting stratum.

### 7.4 Complex similarity used only for curation

For structural datasets, split construction/audit may use native co-complex information such as aligned pocket/contact patterns or ligand RMSD, as CleanSplit does. These data are quarantined in the split-builder environment and never enter strict model features. Benchmark curation and model information boundaries are distinct.

## 8. Split construction

### 8.1 Benchmark-first quarantine

For fixed external tests:

1. freeze test identities and raw-label access;
2. construct all ligand/protein/pocket/provenance relations to candidate train data;
3. remove/quarantine training observations exceeding predeclared relations;
4. inspect transitive components;
5. freeze resulting train/validation manifests; and
6. regenerate nearest-neighbour reports independently.

This avoids choosing a favourable test set after seeing the training distribution.

### 8.2 Union-component dual split

Create graphs G_ligand, G_protein, G_pocket, and G_provenance. Their union connects examples if any forbidden relation holds. Assign entire connected components to one fold.

A giant component is expected at some thresholds. Responses, in order:

1. keep the locked test and remove connected training records;
2. use DataSAIL-style optimization with explicit edge/entity loss;
3. relax only the deployment-inappropriate relation through a committee-approved amendment; or
4. report that the desired dual-novel evaluation is underpowered.

Never split a giant component arbitrarily and call it independent.

### 8.3 Time-forward split

Use earliest public availability of the measurement, not current database ingestion date. Structures, site annotations, and labels must also be available at cutoff. All ligand/protein/provenance overlap is still audited: a time split alone does not prevent analogue or target mirroring.

### 8.4 Nested novelty tiers

Within each test, label—not select after performance—each example by:

- ligand novelty tier;
- protein sequence/domain novelty tier;
- pocket structural novelty tier;
- assay/provenance novelty;
- structure-conformation source; and
- joint novelty.

This produces performance surfaces rather than one misleading score.

## 9. Leakage red-team models

Every benchmark release must be attacked by models that should not solve true molecular recognition:

1. ligand-only ECFP and gMolAI nearest neighbour;
2. target/pocket nearest neighbour;
3. target ID/family/sequence-only model;
4. pocket size/composition-only model;
5. source/database/publication/year/assay-metadata model;
6. missingness, PDB resolution, detector rank, and structure-source model;
7. label lookup/mean by target, scaffold, publication, and nearest train group;
8. combined two-sided k-nearest-neighbour model;
9. full model with pocket permuted within target/family strata;
10. full model with receptor coordinates destroyed but sequence retained; and
11. a train/test membership classifier over inputs and metadata.

Interpretation:

- strong ligand-only performance means the set does not isolate receptor value;
- strong protein-only performance means target label priors remain;
- strong nuisance-only performance reveals source confounding;
- unchanged full performance after pocket permutation/coordinate destruction means the structural branch is not causally used;
- high membership classification suggests distribution/source separation and possible dataset bias.

These tests diagnose shortcuts; none alone proves intentional leakage.

## 10. Preprocessing and software controls

All data-dependent operations are fitted inside the training fold:

- standardization statistics;
- descriptor/feature selection;
- vocabulary/category construction;
- PCA/dimensionality reduction;
- nearest-neighbour indexes;
- class weights and sampling;
- replicate aggregation that uses labels;
- calibration/conformal scores;
- decision/abstention thresholds; and
- learned pocket/ligand clustering used by the model.

Structure and SMILES normalization rules can be globally fixed if they do not learn from the dataset. Their software versions remain identical across folds.

The pipeline records read dependencies. Test-label files are mounted read-only only in the evaluation job. Training jobs fail if an external-label path is present.

## 11. Pretrained-model leakage policy

For every pretrained component record:

- pretraining corpora and date;
- objective and whether affinity/activity/pose labels were used;
- entity identifiers available for overlap checks;
- source/checkpoint revision;
- whether test ligands, proteins, pockets, or complexes are exact/near matches; and
- whether embeddings were fitted/calibrated on project data.

Tracks:

- **entity-exposed, label-free:** allowed with disclosure and novelty sensitivity;
- **same supervised label/assay exposed:** prohibited for confirmatory test;
- **unknown corpus:** allowed only as a separately labelled baseline or after committee risk acceptance;
- **foundation model released after test data:** time-forward claim must state this or use an earlier checkpoint.

For gMolAI, chemistry exposure is expected; the key comparison is whether receptor information adds value beyond gMolAI and whether performance persists in low-similarity bins.

## 12. Benchmark portfolio and what each proves

| Benchmark/protocol | Leakage role | What it can support | What it cannot support alone |
|---|---|---|---|
| random PDBbind | optimistic diagnostic | implementation comparability | generalization |
| official LP-PDBbind | published protein/ligand-controlled split | standardized hard comparison | universal “leak-free” claim |
| official CleanSplit + independent CASF subsets | PDBbind/CASF structural/chemical filtering | CASF comparability with reduced contamination | prospective deployment |
| project union-component split | exact intended information boundary | internal confirmatory hypotheses | protection from unknown upstream errors |
| temporal target-series split | deployment chronology | lead-optimization simulation | target/chemical novelty by itself |
| NTAB | recent time+ligand novelty tiers | ligand-memorization stress test | settled peer-reviewed protein/pocket benchmark |
| BayesBind | pockets selected away from BigBind training | hard virtual screening versus KNN | quantitative affinity |
| CASP16/D3R/SAMPL | blinded labels | strongest external evidence | broad statistical coverage |
| LIT-PCBA | experimentally measured assay compounds but documented preprint audit failures | red-team memorization/analogue diagnostic, or a newly rebuilt split | unmodified external confirmation |
| DUD-E | known decoy/analogue bias | red-team diagnostic | headline screening performance |

## 13. Acceptance criteria

### Literal contamination: zero tolerance

- no exact experimental observation crosses train/validation/test;
- no canonical ligand/target/structure alias reconnects forbidden folds;
- no test label contributes to preprocessing, model selection, calibration, or exclusion;
- no query pose or ligand coordinate enters the strict model;
- no external test record is present in supervised pretraining;
- no post hoc test deletion or best-run selection occurs.

### Similarity-induced shortcut risk: quantified, not hidden

- every test example has continuous nearest ligand, sequence, structure, pocket, assay, and provenance measures;
- results are shown by novelty tiers and target clusters;
- official LP-PDBbind and CleanSplit manifests are reproduced exactly before extensions;
- union-component split passes an independent implementation audit;
- all discarded observations/components and distribution shifts are reported;
- ligand-only, protein-only, nuisance-only, and two-sided KNN baselines are present;
- performance persists, with uncertainty, in the predeclared novelty region needed for the claim.

### Benchmark governance

- test IDs/labels and run count are access-controlled;
- first confirmatory predictions are immutable;
- preprocessing failures are included in denominators;
- pretraining exposure is disclosed;
- benchmark and source versions/hashes are recorded; and
- corrections create a new named benchmark version.

## 14. Leakage audit artifacts

Each dataset build must produce:

- identity_aliases.parquet
- observation_lineage.parquet
- ligand_similarity_edges.parquet
- protein_similarity_edges.parquet
- pocket_similarity_edges.parquet
- assay_provenance_edges.parquet
- target_mirroring_audit.parquet
- union_components.parquet
- split_assignments.parquet
- external_overlap_report.parquet
- pretraining_exposure.tsv
- per_test_nearest_neighbors.parquet
- split_distribution_report.json
- leakage_gate_report.json
- human-readable leakage_report.md

Only permitted metadata and hashes are versioned in Git; restricted payloads remain in controlled storage.

## 15. Required sensitivity analyses

1. vary ligand, protein, and pocket similarity thresholds over a frozen grid;
2. compare union components with DataSAIL’s two-dimensional split;
3. remove nearest training neighbours progressively;
4. stratify by exact, close, moderate, and remote ligand similarity;
5. stratify by global versus local pocket similarity disagreement;
6. remove same-publication/patent/assay-series records;
7. compare before/after cross-source observation deduplication;
8. compare gMolAI-exposed versus lowest-similarity chemistry;
9. evaluate holo, apo, and predicted structures separately;
10. compare ligand-only/full deltas rather than only absolute full-model scores; and
11. reproduce official LP-PDBbind/CleanSplit nearest-neighbour baselines.

If conclusions change materially across reasonable thresholds, report the dependency and narrow the claim.

## 16. Governance rule

Any team member may file a leakage incident. A credible incident freezes the affected benchmark/model-release state until:

1. exact scope is identified;
2. affected artifacts and conclusions are enumerated;
3. corrected data/splits receive new hashes and versions;
4. confirmatory status is reset if necessary; and
5. the committee/evaluation owner approves the correction note.

Leakage discovered after publication or release requires a visible erratum, not a silent artifact replacement.

## 17. Bottom line

iScore3.0 will not seek one impressive benchmark number. It will make structural value falsifiable by measuring what can be predicted from chemistry, target identity, source, and nearest neighbours before crediting the 3D pocket. The strongest result would be a reproducible full-minus-ligand-only gain that persists across low chemical similarity, low local-pocket similarity, time-forward data, apo/predicted receptors, and a genuinely locked external experiment.
