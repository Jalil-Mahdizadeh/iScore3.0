# Literature and Software Review

**Evidence cut-off:** 20 August 2026  
**Question:** What published methods can predict protein–small-molecule affinity or screening activity from protein information and a coordinate-free ligand, and what do their implementations teach us about encoding a 3D binding site?

## 1. Review method

The review used combinations of “pose-free”, “docking-free”, “binding pocket”, “protein 3D graph”, “ligand 2D graph”, “SMILES”, “affinity”, “drug–target affinity”, “virtual screening”, “surface”, and known benchmark names. Searches covered Crossref/DOI pages, PubMed/PMC, journal sites, arXiv/bioRxiv, conference proceedings, GitHub, official dataset sites, and repository references. Forward/backward citation chasing was applied to direct competitors and protein encoders.

A method was treated as a direct precedent when it:

1. predicts affinity or binding/activity;
2. consumes receptor sequence or structure;
3. consumes a small-molecule representation; and
4. does not require a docked query complex as its primary input.

A second, stricter classification asks whether ligand coordinates or pose-derived labels are absent throughout training and inference. The distinction is crucial: “docking-free” does not imply “ligand-3D-free”.

Source repositories were inspected, when public, for actual data interfaces, coordinate use, model modules, checkpoints, dependencies, and licence files. A repository URL alone was not counted as reproducibility. Detailed records are in [evidence/publications.tsv](evidence/publications.tsv) and [evidence/software.tsv](evidence/software.tsv).

The publication registry contains 82 evaluated methods, studies, benchmarks, and scientific resources. The Gate-0/1 update added PSG-BAR, PLMCA, AttentionMGT-DTA, BlendNet, AttentionSiteDTI, BindingSite-AugmentedDTA, PGraphDTA, MMPD-DTA, AlignNet, LigUnity, and CSCo-DTA. The infrastructure/resource papers cover CrossDocked2020, APObind, scPDB, KLIFS, GPCRdb, AlphaFold DB, RCSB PDB, SIFTS, Guide to PHARMACOLOGY, PubChem BioAssay, Drug Target Commons, MUV, DEKOIS 2.0, PDBFlex, CATH, and the public FEP benchmark. Their record-level use decisions are synthesized in [03_dataset_strategy.md](03_dataset_strategy.md).

This is a deep, structured scoping review rather than a registered systematic review. It is designed for architecture and governance decisions; the search log preserves queries and unresolved items.

## 2. Classification

| Class | Protein input | Ligand input | Cross-complex geometry | Core examples | Relevance |
|---|---|---|---|---|---|
| strict pose-free structural affinity | receptor/pocket 3D | graph/SMILES only | none | PSG-BAR, BANANA, CASTER-DTA, 3DProtDTA, AttentionMGT-DTA, PLMCA, Graph_RG, HoloProt affinity branch | closest precedents |
| pose-free inference with privileged 3D training | pocket 3D or sequence | graph/SMILES at inference | no query pose at inference; pose-derived targets or pretrained ligand geometry in training | PLANET, BlendNet | strong architecture precedents, outside strict claim |
| implicit pose/complex or ligand-3D models | pocket/complex 3D | ligand graph plus coordinates/conformer | present or reconstructed | GRIPHIN, T-ALPHA, TankBind, DrugCLIP, PLANTAIN, PocketDTA, DTBind, GEMS, MMPD-DTA, AlignNet, LigUnity | design inspiration and non-core comparators |
| sequence/2D affinity | protein sequence | graph/SMILES | none | DeepDTA, GraphDTA, PSICHIC | essential controls |
| protein representation learning | protein 3D | none | not applicable | GVP, GearNet, ProNet, ProteinWorkshop, MaSIF, dMaSIF, ScanNet | pocket-encoder candidates |
| receptor-only pocket finding | receptor 3D | none | not applicable | P2Rank, fpocket, DeepPocket, GRaSP | secondary Mode B pipeline |

## 3. Principal novelty finding

The generic concept “protein 3D plus ligand 2D without docking” is prior art. In particular:

- BANANA takes a pocket PDB and ligand SMILES and performs pocket–ligand activity prediction.
- PSG-BAR explicitly combines a full-protein 3D residue graph and a 2D ligand graph for affinity regression.
- CASTER-DTA explicitly combines an equivariant 3D protein graph with a standard 2D molecular graph and residue–atom cross-attention.
- 3DProtDTA uses a residue-level protein structure graph with a ligand graph or fingerprint.
- AttentionMGT-DTA combines a receptor pocket graph, ESM features, and a SMILES-derived molecular graph.
- PLMCA combines protein language-model, 3D geometric, and physicochemical features with a ligand molecular graph and cross-attention.
- HoloProt combines hierarchical structural protein representations with a molecular graph for affinity prediction.
- Graph_RG separately encodes a protein pocket and ligand for pose-free affinity prediction.
- PLANET combines a 3D pocket graph and 2D ligand graph, although its auxiliary training targets use complex geometry.

The publishable research gap is therefore narrower: a machine-enforced coordinate-free ligand pathway, including pretraining; high-confidence assay-to-construct-to-site provenance; multi-scale pocket chemistry and geometry whose incremental value is demonstrated beyond ligand/sequence/nearest-neighbour controls; explicit receptor/pocket uncertainty; and leakage-resistant, conformation-stressed evaluation. These components may be individually precedented. Novelty must be framed as a validated scientific system and evidence package, never as ownership of the protein-3D plus ligand-2D modality.

## 4. Direct and near-direct methods

### 4.1 BANANA and BigBind

**Primary sources:** [BigBind paper, Journal of Chemical Information and Modeling, 2024](https://doi.org/10.1021/acs.jcim.3c01211); [BigBind source](https://github.com/molecularmodelinglab/bigbind); [BANANA source](https://github.com/molecularmodelinglab/banana); [BayesBind follow-up](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980085/).

BigBind maps ChEMBL activity records to binding pockets derived from CrossDocked2020, yielding hundreds of thousands of activities across more than one thousand pockets without requiring an experimental pose for each assayed molecule. Its 3D pocket split was explicitly designed to reduce target leakage. BANANA is the associated pose-free classifier: residues are nodes positioned by C-alpha coordinates, nearby residues form the protein graph, the ligand is a 2D molecular graph, and two message-passing branches are fused through an outer-product-like interaction representation.

The public BANANA interface accepts a pocket PDB and a SMILES. Inspection of the model and featurization confirms that the query ligand does not require Cartesian coordinates. Both repositories carry MIT licences at the reviewed revision.

**Strengths**

- Exact input match to the strict iScore3.0 concept.
- A realistic many-assay-to-pocket mapping pipeline.
- Pocket-held-out evaluation and mandatory ligand-only comparisons.
- Lightweight, reproducible source suitable as the first structural baseline.

**Limitations**

- It is primarily a binary activity/virtual-screening model, not a calibrated multi-endpoint affinity model.
- Putative inactives inherit uncertainty from thresholding and negative construction.
- Coarse residue graphs omit much atomic chemistry.
- BayesBind subsequently showed that nearest-neighbour baselines remain extremely difficult to beat on carefully selected novel targets.

**Consequence for iScore3.0:** reproduce BANANA or create an architecture-matched implementation before claiming value from a more elaborate pocket encoder.

### 4.2 CASTER-DTA

**Primary sources:** [Briefings in Bioinformatics paper](https://doi.org/10.1093/bib/bbaf554); [source repository](https://github.com/rachitk/caster-dta).

CASTER-DTA is the clearest recent direct competitor. It represents the protein as an equivariant 3D graph, the compound as a conventional 2D graph, and uses residue–atom cross-attention for affinity prediction. It evaluates on BindingDB, Davis, KIBA, and Metz.

The repository exposes separate protein-structure and molecular-graph preprocessing and does not require a ligand pose for the described core pathway. At review time its licence text imposed research/academic or non-commercial conditions rather than a standard permissive open-source licence; reuse therefore needs a legal check.

**Strengths**

- Nearly identical modality boundary.
- Explicit geometric protein encoding and cross-modal attention.
- Multiple DTA benchmarks and useful ablations.

**Limitations**

- Whole-protein structure encoding is not the same as a known-pocket model.
- Conventional DTA splits can retain protein and ligand analogue leakage.
- Davis and KIBA are narrow or non-thermodynamic benchmarks.
- The custom licence may prevent direct code incorporation.

**Consequence:** iScore3.0 must benchmark against the paper and, if necessary, a clean-room equivalent. Pocket focus, strict splits, endpoint semantics, uncertainty, and apo/predicted testing are differentiators—not the broad architecture.

### 4.3 3DProtDTA

**Primary sources:** [RSC Advances paper](https://doi.org/10.1039/D3RA00281K); [source repository](https://github.com/receptor-ai/3d-prot-dta).

3DProtDTA combines a residue-level protein graph derived from 3D structure with either a ligand molecular graph or an ECFP representation. Source inspection shows separate protein and compound branches followed by representation fusion; ligand coordinates are not required. The reported experiments use Davis and KIBA.

**Strengths**

- Simple strict precedent and useful late-fusion baseline.
- Demonstrates that receptor geometry can be introduced without a complex pose.

**Limitations**

- Narrow benchmark regime and limited evidence of dual-novel generalization.
- Little pocket-specific chemistry or explicit atom–residue compatibility.
- The reviewed repository was small and did not expose a clear root licence file.

**Consequence:** include a 3DProtDTA-style concatenation baseline; do not infer that cross-attention is valuable unless it beats this simpler control.

### 4.4 HoloProt

**Primary sources:** [NeurIPS 2021 paper](https://proceedings.neurips.cc/paper_files/paper/2021/hash/d494020ff8ec181ef98ed97ac3f25453-Abstract.html); [source repository](https://github.com/vsomnath/holoprot).

HoloProt learns a hierarchical protein representation from surface patches, residue structure, and sequence. Its affinity branch combines the protein representation with a ligand molecular graph. It used a filtered PDBbind 2019 refined set and included scaffold and protein-identity splits. The source is MIT-licensed.

**Strengths**

- Strong precedent for multi-scale protein representation.
- Ligand graph does not require a query pose in the affinity branch.
- Scaffold/protein-aware evaluation is more relevant than random pair splitting.

**Limitations**

- Protein preprocessing depends on an operationally heavy surface/electrostatics toolchain including MSMS, PDB2PQR/APBS, DSSP, and related geometry software.
- A whole-protein hierarchy can dilute site-specific signal.
- The underlying affinity labels and holo structures remain vulnerable to PDBbind redundancy.

**Consequence:** borrow the multi-scale principle, not necessarily the complete surface stack. A surface branch must earn its cost through ablation.

### 4.5 Graph_RG

**Primary sources:** [Proteins paper](https://doi.org/10.1002/prot.70010); [CASP16 assessment](https://doi.org/10.1002/prot.70061); [CASP16 affinity data](https://zenodo.org/records/16762332).

Graph_RG separately represents the receptor pocket and ligand and was the leading pose-free affinity approach in the CASP16 pharmaceutical-category analysis. The assessment involved 140 measured affinities across five systems and reported a maximum weighted Kendall rank correlation around 0.42, with experimental noise placing a lower ceiling on achievable agreement than perfect rank correlation.

Targeted searches of the paper, title, author pages, and GitHub did not identify an official public source repository by the evidence cut-off. This is recorded as “not found”, not proof that code is unavailable.

**Strengths**

- Direct pocket/ligand task and unusually valuable blinded evaluation.
- Highlights ranking under real experimental noise rather than retrospective random splits.

**Limitations**

- Reproduction is limited without source.
- Public descriptions do not yet provide enough implementation detail to treat it as an executable baseline.
- A small number of target series produces wide system-level uncertainty.

**Consequence:** use CASP16 as a locked external benchmark and compare at the protocol level. Contact authors if exact reproduction becomes a committee requirement.

### 4.6 PLANET and PLANET v2.0

**Primary sources:** [PLANET paper](https://doi.org/10.1021/acs.jcim.3c00253); [source repository](https://github.com/ComputArtCMCG/PLANET); [PLANET v2.0 preprint](https://arxiv.org/abs/2601.07415); [PLANET v2.0 site](https://www.pdbbind-plus.org.cn/planetv2).

PLANET uses an E(3)-equivariant pocket encoder, a graph-attention ligand encoder, and a protein–ligand interaction module. It jointly predicts affinity, protein–ligand contacts, and intraligand distances. Source inspection identified separate ProteinEGNN, LigandGAT, and interaction components. The model’s ligand message passing is graph based, but complex geometry supplies auxiliary contact/distance supervision during training. Its public workflow also commonly uses structure files for molecules.

PLANET v2.0 introduces probabilistic contact/distance modelling and a distance-energy formulation. It is highly current and should be rechecked during implementation.

**Strengths**

- Closely relevant equivariant pocket and pairwise fusion design.
- Multi-task supervision can improve data efficiency.
- Broad benchmarking across affinity and virtual-screening tasks.

**Limitations for the strict iScore3.0 claim**

- Training receives ligand/complex 3D-derived privileged information.
- It can learn geometry unavailable to a truly coordinate-free ligand pipeline.
- No clear root licence file was visible in the reviewed PLANET repository, so code reuse is not assumed.

**Consequence:** compare PLANET as a privileged-information upper or external comparator. Do not copy its pose-derived auxiliary losses into the strict core.

### 4.7 PSG-BAR

**Primary sources:** [Pandey et al., Molecules 2022](https://doi.org/10.3390/molecules27165114); [source repository](https://github.com/diamondspark/PSG-BAR).

PSG-BAR is an unambiguous direct precedent. It encodes a protein’s folded structure as a residue graph, a ligand as a DeepChem two-dimensional molecular graph, and predicts continuous affinity through residual graph-attention branches and an interaction-scoring module. The paper reports PDBbind, BindingDB, KIBA, and Davis experiments under warm, cold-drug, cold-protein, and cold-both settings. The cold-both drop is scientifically more informative than its warm score: for example, the reported BindingDB MSE/Pearson moves from approximately 0.651/0.864 in the warm setting to 2.102/0.515 in cold both, and PDBbind from 1.660/0.762 to 2.100/0.599.

The reviewed source was pinned at `6540ab8ccccdb543ada7bc8b51a01a171d5c3786`. It confirms a 2D ligand pathway but raises three reproducibility issues: the paper describes five nearest protein neighbours whereas `GraphProcessing.py` sets `k=3`; raw Cartesian coordinates are concatenated into node features processed by a conventional GAT, so rigid-motion invariance is not guaranteed; and no root licence was visible. These issues limit direct code reuse, not PSG-BAR’s status as prior art.

**Consequence:** PSG-BAR must appear in every novelty comparison and should be reproduced only through a clean-room, invariant implementation if its protocol is used as a baseline.

### 4.8 PLMCA and AttentionMGT-DTA

**PLMCA:** [Journal of Medicinal Chemistry, 2026](https://doi.org/10.1021/acs.jmedchem.5c03431).

PLMCA integrates two protein language models, three-dimensional geometric and physicochemical protein features, a ligand molecular graph, and cross-attention. It jointly addresses pocket identification and affinity and includes assay-condition variables for its ChEMBL experiments. The paper reports random, unseen-ligand, and unseen-protein PDBbind21 tests for Kd/Ki; ChEMBL_mini R-squared values of 0.531, 0.635, and 0.519 for IC50, Kd, and Ki; and pocket AUPR up to 0.655. Targeted title, author, DOI, and GitHub searches found no official public implementation by the cut-off. Because the accessible interface description does not expose the full preprocessing and pretrained-data lineage, its strict no-ligand-3D status remains provisional rather than assumed.

**AttentionMGT-DTA:** [Neural Networks, 2024](https://doi.org/10.1016/j.neunet.2023.11.018); [source repository](https://github.com/JK-Liu7/AttentionMGT-DTA).

AttentionMGT-DTA constructs receptor binding-pocket graphs from AlphaFold/PDB structures, combines residue geometry with ESM representations, represents the drug as a SMILES-derived graph, and learns protein/ligand interactions with attention. The pinned source (`e94a28dad3642abab82f353f799aa1246e7ab0dc`) confirms the strict query-ligand modality. It also reveals important evaluation limitations: the training loop evaluates the test set each epoch, uses test MSE for the learning-rate scheduler and best-epoch selection, and relies on conventional Davis/KIBA folds. In `protein_process.py`, multiple ConvexHull pockets are iterated but only the final bounding-box values are subsequently used. No root licence was visible.

**Consequence:** both methods eliminate any residual broad architecture claim. Their assay-context and multimodal ideas are useful, but iScore3.0 must distinguish itself through verified lineage and leakage-controlled evidence.

### 4.9 Adjacent methods that refine the boundary

- [BlendNet](https://doi.org/10.1093/bib/bbae712) predicts affinity from pocket/protein sequences and a 2D compound graph at inference. Its teacher is trained on PLIP atom–residue interaction labels extracted from experimental complexes, so it is a privileged-training comparator, not strict evidence.
- [AttentionSiteDTI](https://doi.org/10.1093/bib/bbac272) is a direct receptor-site-graph plus ligand-graph precedent for binary DTI. [BindingSite-AugmentedDTA](https://doi.org/10.1093/bib/bbad136) then uses such binding-site information to augment DTA predictors, but its reviewed repository is primarily data tables rather than a complete executable pipeline.
- [PGraphDTA](https://arxiv.org/abs/2310.04017) and [CSCo-DTA](https://doi.org/10.1093/bib/bbad512) add protein contact-map structure proxies to 2D-ligand affinity models. They belong in structure/sequence-control discussions, but contact maps do not establish value from a chemically detailed experimental pocket.
- [MMPD-DTA](https://doi.org/10.1021/acs.jcim.4c01528) includes a pocket–drug graph with intermolecular edges and therefore requires complex/pose geometry; its source consumes precomputed graphs without releasing their construction.
- [AlignNet](https://doi.org/10.1093/bioinformatics/btag599) aligns ESM/GearNet and MolFormer/GraphMVP modalities, but the released pipeline reads ligand MOL2 coordinates and builds 5-Angstrom protein–ligand edges. “Structure-agnostic” in that work does not mean ligand-3D-free.
- [LigUnity](https://doi.org/10.1016/j.patter.2025.101371) is a high-priority docking alternative and shared pocket–ligand foundation model, but its Uni-Mol pipeline generates or consumes ligand conformers and its case-study preparation uses bound ligand coordinates to define sites.

These cases demonstrate why input classification must follow executable preprocessing, pretrained lineage, and auxiliary labels rather than titles such as “complex-free”, “structure-agnostic”, or “docking-free”.

### 4.10 Summary of direct methods

| Method | Protein 3D | Ligand coordinates needed by core claim? | Query complex distances? | Affinity or activity | Public source | Strict comparator |
|---|---:|---:|---:|---|---:|---:|
| PSG-BAR | whole protein | no | no | affinity | yes | yes, clean-room due source caveats |
| BANANA | pocket | no | no | binary activity | yes | yes |
| CASTER-DTA | whole protein | no | no | affinity | yes | yes |
| 3DProtDTA | protein graph | no | no | affinity | yes | yes |
| AttentionMGT-DTA | predicted/experimental pocket | no | no | affinity | yes | yes, with evaluation caveats |
| PLMCA | protein multimodal 3D | described as graph-only ligand | no in stated interface | affinity and pocket | not found | provisional paper-level |
| HoloProt | hierarchical protein | no for ligand graph branch | no | affinity/function | yes | yes, with preprocessing caveats |
| Graph_RG | pocket | described as pose-free | no | affinity/ranking | not found | protocol-level |
| PLANET | pocket | no at inference claim; 3D privileged training | learned contact/distance targets | affinity/activity | yes | no, privileged comparator |
| BlendNet | pocket sequence | no at inference | PLIP-derived teacher targets | affinity | yes | no, privileged comparator |
| PLANET v2.0 | pocket | inspect final implementation before classification | probabilistic geometry | affinity | site/preprint | surveillance |

## 5. Related methods that do not satisfy the strict boundary

### 5.1 GRIPHIN

**Sources:** [Journal of Cheminformatics paper](https://doi.org/10.1186/s13321-026-01203-8); [MIT-licensed source](https://github.com/molinfo-vienna/GRIPHIN).

GRIPHIN encodes a protein pocket as a ten-channel GRAIL pharmacophore grid and a ligand with a graph-attention network. Its name and abstract can look strictly pose-free, but the paper and source use ligand Cartesian positions/Fourier encodings and a ligand-centred spatial construction. It is therefore not a strict comparator. Its receptor pharmacophore channels are valuable inspiration if reconstructed solely from the receptor and centred independently of the query ligand.

### 5.2 T-ALPHA

**Sources:** [JCIM paper](https://doi.org/10.1021/acs.jcim.4c02332); [MIT-licensed source](https://github.com/gregory-kyro/T-ALPHA).

T-ALPHA combines protein, ligand, and complex channels using surface/equivariant/sequence features and hierarchical cross-attention. It uses ligand and co-complex 3D information. Its censored-label treatment, multimodal ablations, LP-PDBbind/BDB2020+ evaluation, and uncertainty design are relevant; its input boundary is not.

### 5.3 TankBind

**Sources:** [NeurIPS 2022 paper](https://proceedings.neurips.cc/paper_files/paper/2022/hash/2f89a23a19d1617e7fb16d4f7a049ce2-Abstract-Conference.html); [MIT-licensed source](https://github.com/luwei0917/TankBind).

TankBind represents protein blocks and ligand atoms, predicts interatomic distance maps, and optimizes a ligand pose. Ligand conformer/internal-distance information and pose generation make it a binding-structure predictor rather than the target strict affinity system. Its pairwise tensor and triangle-update concepts motivate a later pose-free compatibility module.

### 5.4 DrugCLIP and ProFSA

**Sources:** [DrugCLIP NeurIPS paper](https://proceedings.neurips.cc/paper_files/paper/2023/hash/8bd31288ad8e9a31d519fdeede7ee47d-Abstract-Conference.html); [DrugCLIP source](https://github.com/bowen-gao/DrugCLIP); [ProFSA source](https://github.com/bowen-gao/ProFSA).

DrugCLIP contrastively aligns pockets and molecules for screening. Its distributed molecule data contain RDKit conformers and its Uni-Mol lineage is three-dimensional. ProFSA pretrains on millions of fragment–surrounding pseudo-complexes and is valuable for pocket representation learning, but its ligand/fragment geometry and pseudo-complex construction put it outside the strict pathway.

### 5.5 PLANTAIN

**Source:** [PLANTAIN paper and source links](https://pmc.ncbi.nlm.nih.gov/articles/PMC10402188/).

PLANTAIN begins from a 3D pocket and 2D ligand graph, then predicts distance-dependent interaction fields and optimizes a ligand pose. It is useful evidence that pairwise compatibility can be built before coordinates are known, but its endpoint includes pose reconstruction and is not the proposed core task.

### 5.6 PocketDTA

**Sources:** [Bioinformatics paper](https://doi.org/10.1093/bioinformatics/btae594); [MIT-licensed source](https://github.com/ZhaoLongSYSU/PocketDTA).

PocketDTA selects predicted pockets, encodes protein information with ESM/GVP-like components, and applies bilinear attention to drug features. Its reported ligand representation uses GraphMVP, whose pretraining incorporates 3D molecular information. It is valuable for unknown-pocket aggregation but fails the strict pretrained-information test. Its DoGSite3 dependency also requires scrutiny for academic/commercial use.

### 5.7 DTBind

**Sources:** [Research paper](https://doi.org/10.34133/research.1022); [source](https://github.com/liqy09/DTBind).

DTBind uses protein-surface and interaction-oriented components with MSMS/PyMesh/PLIP-style preprocessing. Its affinity mode consumes complex structure information. It informs surface encoding and engineering risks, not strict scoring.

### 5.8 GEMS

**Sources:** [Nature Machine Intelligence paper](https://doi.org/10.1038/s42256-025-01124-5); [source](https://github.com/camlab-ethz/GEMS); [artifacts](https://doi.org/10.5281/zenodo.14260170).

GEMS analyses structural and chemical redundancy in PDBbind and introduces CleanSplit. Its own geometry-enhanced interaction representation uses a bound complex, so it is not a strict model. Its most important contribution here is methodological: apparent performance can be dominated by train–test similarity, and rigorous split metadata can reverse model rankings.


### 5.9 Data-leakage and shortcut literature

Leakage findings materially change how every method above should be interpreted. [Wallach and Heifets](https://doi.org/10.1021/acs.jcim.7b00403), [Sieg et al.](https://doi.org/10.1021/acs.jcim.8b00712), and the [DUD-E hidden-bias study](https://doi.org/10.1371/journal.pone.0220113) show that chemical redundancy and decoy construction can reward memorization. [Predicting or Pretending](https://doi.org/10.3389/fphar.2020.00069), [Volkov et al.](https://doi.org/10.1021/acs.jmedchem.2c00487), and [Latent Biases](https://doi.org/10.1021/acsomega.2c06781) show why ligand-only, protein-only, and nearest-neighbour controls are indispensable for PDBbind-like affinity data.

[LP-PDBbind](https://doi.org/10.1021/acs.jpcb.5c08598), CleanSplit/GEMS, [DataSAIL](https://doi.org/10.1038/s41467-025-58606-8), PLINDER, BigBind/BayesBind, and the recent NTAB/target-mirroring preprint offer complementary remedies; they do not agree on one universal threshold. [Li et al.](https://doi.org/10.1093/bib/bbab225) provide an important counterpoint: similar training data can be legitimate and useful when the intended deployment is interpolation. The correct response is deployment-matched evaluation plus continuous novelty reporting, not a ritual claim of being “leak-free”.

A 2025 reproducible [LIT-PCBA audit preprint](https://arxiv.org/abs/2507.21404) reports exact cross-split/query leakage and extensive analogue redundancy, so unmodified LIT-PCBA is a diagnostic rather than confirmatory benchmark in this plan. The complete threat taxonomy and executable release gates are in [08_data_leakage_threat_model.md](08_data_leakage_threat_model.md).

### 5.10 Other surveillance items

IPBind and several 2025–2026 preprints use bound/unbound complex geometry, ligand coordinates, or pose reconstruction despite “pose-free” language in adjacent tasks. They remain in the search log for periodic reassessment. Every new method must be classified from its actual data pipeline, not its abstract alone.

## 6. Sequence and ligand baselines

### DeepDTA

[DeepDTA](https://doi.org/10.1093/bioinformatics/bty593) learns from protein and SMILES character sequences. It is historically important and inexpensive, but weak splits and narrow benchmarks can exaggerate progress.

### GraphDTA and DGraphDTA

[GraphDTA](https://github.com/thinng/GraphDTA) replaces the ligand string with a molecular graph. [DGraphDTA](https://github.com/595693085/DGraphDTA) adds a predicted protein contact graph. These are appropriate controls for the incremental value of experimentally observed receptor geometry.

### PSICHIC

[PSICHIC](https://doi.org/10.1038/s42256-024-00847-1) uses protein sequence and ligand molecular information without a 3D structure. Its public [Apache-licensed source](https://github.com/huankoh/PSICHIC) makes it a strong modern sequence-only comparator and fallback.

### Ligand-only controls

A gMolAI-only regressor/classifier, ECFP with ridge/random forest/gradient boosting, and a matched-capacity ligand GNN are mandatory. Pocket value is not demonstrated if these controls are missing or evaluated on different splits.

## 7. Protein pocket encoders

### 7.1 Residue equivariant graphs

**GVP-GNN:** [ICLR 2021 paper](https://arxiv.org/abs/2009.01411) and [source](https://github.com/drorlab/gvp-pytorch). Geometric Vector Perceptrons maintain scalar and vector channels, naturally combining residue identity with directions and local frames. They are mature, understandable, and a strong first choice.

**GearNet:** [ICLR 2023 paper](https://arxiv.org/abs/2203.06125) and [source](https://github.com/DeepGraphLearning/GearNet). Relational edges and geometric pretraining provide a strong alternative with good ecosystem support.

**ProNet:** [paper](https://arxiv.org/abs/2207.12600) and [DIG source](https://github.com/divelab/DIG). It encodes hierarchical protein geometry efficiently using distances, angles, and torsions.

**ProteinWorkshop:** [ICLR 2024 project](https://github.com/a-r-j/ProteinWorkshop). It provides a configurable benchmark framework for protein geometric models and can reduce reimplementation risk, although its complete dependency/licence graph must be frozen.

**Recommendation:** begin with one GVP-style residue encoder and one simpler invariant distance-graph baseline. Do not compare many encoders until the data split is locked.

### 7.2 All-heavy-atom pocket graphs

Residue centroids and C-alpha graphs blur donor/acceptor geometry, metals, aromatic faces, and side-chain orientation. A local all-heavy-atom graph is therefore recommended in parallel. Node features should be receptor-derived chemical types; edges use invariant distances and equivariant directions within a fixed radius. Pool atom states to residues or let the ligand attend to both resolutions.

The model must distinguish missing/repaired atoms and must not infer atom types from the query ligand’s contacts.

### 7.3 Surface encoders

**MaSIF:** [Nature Methods paper](https://doi.org/10.1038/s41592-019-0666-6) and [Apache-licensed source](https://github.com/LPDI-EPFL/masif). It learns geodesic surface patch descriptors with chemical/geometric channels.

**dMaSIF:** [CVPR 2021 paper](https://openaccess.thecvf.com/content/CVPR2021/papers/Sverrisson_Fast_End-to-End_Learning_on_Protein_Surfaces_CVPR_2021_paper.pdf) and [source](https://github.com/FreyrS/dMaSIF). It avoids explicit mesh construction and is much faster, but the reviewed repository’s non-commercial/no-derivatives terms require legal review before reuse.

**ScanNet:** [Nature Methods paper](https://doi.org/10.1038/s41592-022-01490-7) and [source](https://github.com/jertubiana/ScanNet). It learns local atomic/residue neighbourhoods for binding-site prediction and supplies relevant geometric primitives.

Surface methods encode shape, curvature, electrostatics, and hydrophobic patches well, but add protonation, triangulation, solver, and licence fragility. The blueprint therefore defers them until residue+atom models establish a reproducible baseline.

### 7.4 Pharmacophore and field representations

GRIPHIN/GRAIL-style channels suggest receptor-only grids for hydrogen-bond donors/acceptors, hydrophobes, aromaticity, charge, exclusion volume, and metals. A receptor-centred field can be sampled by a sparse 3D CNN or equivariant point network. It is a useful ablation, but voxel resolution, rotational handling, protonation, and grid-centre leakage must be documented.

### 7.5 Pocket discovery

**P2Rank:** [Journal of Cheminformatics paper](https://doi.org/10.1186/s13321-018-0285-8); [MIT source](https://github.com/rdk/p2rank). Fast, receptor-only, reproducible, and suitable for top-k Mode B.

**fpocket:** [BMC Bioinformatics paper](https://doi.org/10.1186/1471-2105-10-168); [MIT source](https://github.com/Discngine/fpocket). A classical alpha-sphere baseline with low operational cost.

**DeepPocket:** [JCIM paper](https://doi.org/10.1021/acs.jcim.1c00799); [source](https://github.com/devalab/DeepPocket). A learned rescoring/refinement approach.

**GRaSP:** [source](https://github.com/charles-abreu/GRaSP). A graph-based site predictor relevant to predicted structures.

P2Rank and fpocket should be the initial top-k detectors. Pocket detection recall and affinity performance conditional on correct/incorrect pocket must be reported separately.

## 8. gMolAI v2.0 audit

**Source:** [gMolAI v2.0 repository](https://github.com/Jalil-Mahdizadeh/gMolAI-v2.0).

The reviewed model is a residual GINE-style 2D molecular graph encoder pretrained on approximately 223 million ZINC/PubChem molecular graphs. The public representation concatenates a 256-dimensional raw graph block and a 128-dimensional mean-node block, followed by released calibration/standardization; the mean-node block is weighted more strongly in its similarity workflow. The model implementation also exposes atom/node states, enabling both atom-level cross-attention and a stable global ligand vector.

This is unusually well matched to iScore3.0:

- atom tokens can participate in pose-free residue–atom compatibility;
- the calibrated global embedding supports strong ligand-only and late-fusion baselines;
- all production features can be derived from SMILES/2D connectivity; and
- one encoder can serve both retrieval/leakage diagnostics and prediction.

Before integration, the project must:

1. add or verify an explicit licence covering code, checkpoints, calibrators, and embeddings;
2. pin checkpoint, source revision, feature schema, standardization vectors, and atom ordering;
3. verify by test that no optional path constructs conformers;
4. document the training-corpus relationship to evaluation ligands;
5. freeze the encoder for initial experiments; and
6. compare frozen, adapter, and discriminative fine-tuning without contaminating locked tests.

The 384-dimensional released vector and internal atom embeddings serve different roles and should not be conflated.

The Gate-0/1 adapter audit completed these technical checks on the pinned revision and checkpoint. Across 80 unique pilot ligands and 1,971 atoms, the adapter verified canonical/input atom bijections, 48-dimensional atom and 15-dimensional bond inputs, 128-dimensional `node_z`, 256-dimensional raw graph states, and the released 384-dimensional molecule representation. Exact repeated and equivalent-SMILES CPU runs were bitwise identical; CPU/GPU maximum absolute deviation was at most `9.54e-06`, within the frozen tolerance. Labels were quarantined from the encoder path. The unresolved issues are upstream licensing and the absence of a pretraining identity ledger: the source documents 223,180,699 deduplicated ZINC/PubChem graphs, so exact pilot-entity exposure is unknown even though affinity-label pretraining was not detected. See [the Gate-0/1 audit](../reports/gate01/evidence/gmolai-adapter-audit-v1.json).

## 9. Source-code audit summary

| Project | Interface finding | Licence finding at review | Reuse recommendation |
|---|---|---|---|
| PSG-BAR | full-protein 3D graph + 2D ligand; raw xyz in standard GAT | no clear root licence found | clean-room comparator; fix invariance and k discrepancy |
| BANANA | pocket PDB + SMILES; strict ligand 2D | MIT | reproduce/pin |
| BigBind | assay-to-pocket mapping and 3D pocket splits | MIT | adapt curation concepts |
| CASTER-DTA | protein 3D + molecular graph, cross-attention | custom/restricted terms | clean-room baseline unless cleared |
| 3DProtDTA | separate protein graph and ligand graph/fingerprint | no clear root licence found | reimplement concepts only |
| HoloProt | hierarchical protein + molecular graph | MIT | reuse selectively |
| PLANET | pocket EGNN + ligand GAT + pair module; pose-derived targets | no clear root licence found | conceptual comparator only |
| PLMCA | 3D/PLM/physicochemical protein + ligand graph; implementation unavailable | not assessable | paper-level comparison; contact authors |
| AttentionMGT-DTA | AlphaFold/PDB pocket + ESM + 2D ligand; test-selected training | no clear root licence found | clean-room concepts only |
| BlendNet | 2D ligand/pocket sequence inference; PLIP complex teacher | no clear root licence found | privileged comparator only |
| AttentionSiteDTI | 3D receptor-site graph + 2D ligand; binary endpoint | CC BY 4.0 file | diagnostic/clean-room baseline |
| BindingSite-AugmentedDTA | site-augmented DTA tables; incomplete executable release | CC BY 4.0 file | mapping reference only |
| PGraphDTA | PLM + contact-map structure proxy | repository unavailable at cut-off | paper-level surveillance |
| MMPD-DTA | precomputed pose-dependent pocket–drug graph | no clear root licence found | non-strict comparator only |
| AlignNet | ligand MOL2 coordinates and cross-complex edges in released pipeline | MIT | non-strict comparator only |
| LigUnity | Uni-Mol ligand coordinates/conformers and ligand-defined sites | Apache-2.0 code; data terms separate | non-strict comparator only |
| CSCo-DTA | protein contact maps and interaction-network graph + 2D ligand | no clear root licence found | adjacent proxy baseline |
| GRIPHIN | ligand positions enter model/grid pathway | MIT | receptor-channel inspiration only |
| T-ALPHA | protein, ligand, and complex 3D channels | MIT | censoring/uncertainty inspiration |
| TankBind | ligand distance/conformer and pose optimization | MIT | pair-module inspiration |
| DrugCLIP | conformer-bearing molecule data | inspect transitive terms | non-core comparator |
| ProFSA | pseudo-complex fragment geometry | MIT | pocket-pretraining study only |
| PocketDTA | GraphMVP ligand pretraining includes 3D | MIT; detector terms separate | non-strict comparison |
| GEMS | complex geometry; CleanSplit metadata | repository licence to pin | use split/evaluation assets |
| gMolAI v2.0 | 2D GINE; atom and global states | no clear root licence found | resolve before incorporation |
| P2Rank | receptor-only pocket predictions | MIT | preferred detector |
| APObind | matched apo/holo construction; released scripts use co-crystal coordinates, hard-coded paths, BLAST/PyMOL/ProDy/RDKit | no clear root licence found | audit released pairs; reimplement strict mapping |
| CrossDocked2020 | >52m-file pose archive; pK mixes endpoints and some training labels are pose-conditioned | no clear root licence found for reviewed data/model repository | receptor/pocket-only extraction |
| PDBe SIFTS | residue-level structure–sequence mapping pipeline with MMseqs2/BLASTP, FASTA36, DuckDB, and mmCIF export | Apache-2.0 | pin as mapping tool plus construct checks |
| fpocket | receptor-only geometry | MIT | classical detector |
| MaSIF | mesh/surface pipeline | Apache-2.0 | optional later branch |
| dMaSIF | point/surface pipeline | restrictive terms observed | architecture inspiration only |
| PSICHIC | protein sequence + ligand 2D | Apache-2.0 | strong baseline |

A current licence inventory must be regenerated at each dependency lock. Publication access and source availability do not grant permission to copy code or redistribute trained weights.

## 10. Architecture implications

The review supports the following choices:

1. **Residue plus atom pocket encoding:** residue graphs are proven and efficient; local atoms recover chemical detail.
2. **No surface dependency in the first milestone:** surface methods are scientifically attractive but operationally expensive.
3. **gMolAI atom/global ligand branch:** this is a project-specific asset, provided its licence and provenance are resolved.
4. **Cross-attention after simple fusion baselines:** prior work already shows cross-modal attention, so value must be empirical.
5. **No pose-derived auxiliary task in the strict core:** use masked-pocket, residue-property, pocket-contrastive, or affinity multitask objectives instead.
6. **Pocket and conformation uncertainty:** current precedents rarely make this a first-class output.
7. **Hard evaluation as part of the contribution:** BigBind/BayesBind, GEMS/CleanSplit, LP-PDBbind, and CASP16 show why conventional scores are insufficient.

## 11. Open literature questions

- Obtain the exact Graph_RG implementation or sufficient author clarification for reproduction.
- Obtain PLMCA code/preprocessing or author clarification before treating its strict information boundary as verified.
- Reclassify PLANET v2.0 when final code and peer-reviewed details are available.
- Confirm transitive licences for all surface/electrostatics executables and pretrained weights.
- Audit whether any candidate ligand encoder was pretrained with coordinates, even if current inference accepts SMILES.
- Track post-August-2026 pose-free pocket/ligand models before freezing manuscript novelty language.
- Determine whether strict exclusion of pose-derived training targets should remain the headline track or whether a separately labelled privileged-information variant is scientifically worthwhile.

## 12. Review conclusion

The protein binding site should not be reduced to one fixed vector prematurely. The most defensible initial encoder is a multi-scale receptor-only representation: a residue-level equivariant graph plus a local all-heavy-atom graph, optionally augmented later by a receptor-only pharmacophore field. Its states can interact with gMolAI atom tokens through pose-free attention while a global pocket–ligand pathway captures assay-level priors.

The decisive experiment is not whether this architecture fits PDBbind. It is whether pocket geometry adds statistically reliable information beyond ligand and target priors on locked, dual-novel, conformation-stressed, and externally sourced targets.
