# Search and Review Log

## Review identity

- Project: iScore3.0 blueprint
- Review cut-off: 20 August 2026
- Review type: structured deep scoping review with source-code and dataset audit
- Primary question: methods that combine receptor 3D information with a coordinate-free ligand for affinity/activity prediction
- Secondary questions: protein pocket encoders, scalable datasets, leakage-resistant evaluation, licensing, and reproducibility
- Reviewer output: 02_literature_and_software_review.md, 03_dataset_strategy.md, 08_data_leakage_threat_model.md, and TSV evidence tables

This was not a preregistered systematic review and does not claim exhaustive bibliometric coverage. It was intentionally broad, iterative, and architecture-oriented.

## Sources searched

- DOI/Crossref landing pages and journal sites
- PubMed and PubMed Central
- arXiv and bioRxiv
- NeurIPS proceedings and OpenReview
- Nature, ACS, Oxford Academic, RSC, Springer/BMC, Wiley, and PLOS pages
- GitHub repositories, README files, licences, model/preprocessing modules, issues, and release pages
- official PDBbind+, BindingDB, ChEMBL, PLINDER, BioLiP2, DUD-E, LIT-PCBA, ATOM3D, TDC, CASP/Zenodo, Figshare, and Drug Design Data Resource pages
- Zenodo and Figshare deposits
- cited references and citing/related works for direct competitors

Search results from secondary pages were used for discovery only when a primary paper, repository, dataset page, or deposit could be identified.

## Core eligibility questions

For every method:

1. Does it predict affinity, activity, ranking, or virtual-screening outcome?
2. Does it use protein sequence, whole structure, a pocket, a surface, or a co-complex?
3. Is the query ligand represented by SMILES/2D graph, a conformer, a pose, or a complex?
4. Are ligand/complex coordinates used only at inference, only during training, or both?
5. Are pose-derived contacts/distances auxiliary targets?
6. Is the binding site supplied independently of the query ligand?
7. Are source, data, checkpoints, environment, and licence available?
8. Which datasets/splits were used, and do they control ligand and protein/pocket similarity?
9. What architecture or evaluation lesson transfers to iScore3.0?
10. Is it a strict comparator, privileged comparator, noncompliant inspiration, baseline, encoder, detector, or dataset resource?

## Search waves

### Wave 1 — exact concept and direct competitors

Example queries:

- "protein 3D ligand 2D graph binding affinity"
- "pose-free protein ligand affinity pocket SMILES"
- "docking-free scoring function receptor structure molecular graph"
- "equivariant protein graph standard molecular graph drug target affinity"
- "binding pocket PDB SMILES activity prediction"
- "cross attention protein residues ligand atoms affinity no docking"
- "separate protein pocket ligand graph affinity pose-free"
- "site:github.com protein pocket SMILES affinity"

Key inclusions:

- BANANA/BigBind
- CASTER-DTA
- 3DProtDTA
- HoloProt
- Graph_RG
- PLANET and PLANET v2.0

Key conclusion: the broad modality combination is already published; novelty must be narrower.

### Wave 2 — adjacent complex, pose, and screening models

Example queries:

- "protein pocket ligand graph virtual screening contrastive"
- "protein pharmacophore grid graph attention affinity"
- "protein ligand pair representation distance map pose prediction"
- "multimodal protein surface ligand affinity cross attention"
- "pocket DTA predicted binding site GVP ligand graph"
- "fragment surrounding pretraining protein pocket"
- "protein surface drug target binding affinity GitHub"

Key inclusions:

- GRIPHIN
- T-ALPHA
- TankBind
- DrugCLIP
- ProFSA
- PLANTAIN
- PocketDTA
- DTBind
- GEMS

Classification required inspection beyond abstracts because several apparently pose-free methods use ligand conformers, ligand-centred grids, co-complex channels, or pose-derived targets.

### Wave 3 — protein pocket encoders

Example queries:

- "protein geometric vector perceptron GitHub"
- "GearNet protein geometric pretraining"
- "residue atom protein graph equivariant encoder benchmark"
- "protein molecular surface deep learning MaSIF dMaSIF"
- "local atomic environment protein binding site ScanNet"
- "protein structure representation benchmark ProteinWorkshop"
- "protein pocket pharmacophore point cloud encoder"

Key inclusions:

- GVP-GNN
- GearNet
- ProNet
- ProteinWorkshop
- MaSIF
- dMaSIF
- ScanNet
- GRAIL/GRIPHIN receptor pharmacophore channels

Key conclusion: begin with residue equivariance plus local heavy atoms; defer surface preprocessing until it earns its cost.

### Wave 4 — pocket detection

Example queries:

- "receptor only ligand binding site prediction P2Rank"
- "fpocket open source pocket detection"
- "DeepPocket ligand binding site GitHub"
- "graph protein pocket detection AlphaFold"
- "binding site detection top k predicted protein structure"

Key inclusions:

- P2Rank
- fpocket
- DeepPocket
- GRaSP

Key conclusion: supplied-site evaluation must remain primary; pocket detection is an independently measured secondary task.

### Wave 5 — sequence and 2D baselines

Example queries:

- "drug target affinity protein sequence molecular graph baseline"
- "DeepDTA GraphDTA DGraphDTA source"
- "PSICHIC protein ligand interaction fingerprint source"
- "ligand only nearest neighbor binding affinity benchmark"

Key inclusions:

- DeepDTA
- GraphDTA
- DGraphDTA
- PSICHIC
- ECFP/gMolAI and nearest-neighbour controls

### Wave 6 — structural affinity datasets

Example queries:

- "PDBbind 2025 release number complexes"
- "PDBbind CleanSplit GEMS GitHub"
- "Leak Proof PDBbind LP-PDBBind split"
- "PDBbind Opt HiQBind BioLiP2 Opt dataset"
- "ATOM3D ligand binding affinity dataset"
- "MISATO protein ligand dataset"
- "Binding MOAD update binding data"

Verified current points:

- PDBbind v2025 release announcement and 29,001 protein–ligand complex count on PDBbind+
- PDBbind v2020R1 free-registration announcement
- official GEMS CleanSplit JSON/filter/search assets
- peer-reviewed LP-PDBbind DOI and official metadata repository
- BioLiP2-Opt Figshare deposit
- PLINDER current release/known issues

### Wave 7 — large activity databases

Example queries:

- "BindingDB current measurements compounds targets 2026 download"
- "BindingDB licence CC BY 4 publication date curation date"
- "ChEMBL current release 36 DOI licence"
- "Papyrus standardized bioactivity Zenodo"
- "BigBind counts pockets compounds activities"
- "target assay mapping PDB structure binding site"

Verified current points:

- BindingDB 202608 download counts and July 2026 update
- ChEMBL 36, July 2025, DOI 10.6019/CHEMBL.database.36
- ChEMBL CC BY-SA 3.0
- BigBind paper/repository mapping
- Papyrus parent-source overlap

### Wave 8 — external benchmarks

Example queries:

- "CASP16 ligand binding affinity pharmaceutical category data"
- "BayesBind benchmark EF Bayesian enrichment"
- "LIT-PCBA active inactive counts 15 targets"
- "DUD-E 102 targets 22886 actives decoy bias"
- "CASF-2016 285 complexes scoring ranking"
- "D3R SAMPL affinity challenge dataset"
- "Merck FEP benchmark 264 ligands 8 targets"

Key conclusion after the later leakage wave: CASP16 and BayesBind provide complementary external evidence; unmodified LIT-PCBA is diagnostic pending independent reproduction/rebuild, while CASF and DUD-E cannot be headline proof.

### Wave 9 — leakage and benchmark bias

Example queries:

- "protein ligand binding affinity data leakage PDBbind CASF"
- "PDBbind CleanSplit performance drop nearest neighbor"
- "LP-PDBbind leak proof protein ligand split"
- "Most ligand-based benchmarks memorization generalization"
- "DUD-E hidden bias deep learning"
- "In need of bias control virtual screening"
- "DataSAIL two dimensional interaction split"
- "protein pocket similarity low sequence identity leakage"
- "target mirroring protein ligand affinity"
- "Novelty-Tiered Affinity Benchmark"
- "HonestAffinity leak-aware pocket priors"
- "Data Leakage and Redundancy in the LIT-PCBA Benchmark"

Key inclusions:

- Wallach and Heifets 2018
- Sieg et al. 2019
- Chen et al. 2019
- Su et al. 2020
- Li et al. 2021 counterargument
- GEMS/CleanSplit 2025
- DataSAIL 2025
- LP-PDBbind 2026
- NTAB/systematic-leakage bioRxiv preprint 2026
- HonestAffinity preprint 2026
- APoc/pocket similarity evidence
- PLINDER similarity-aware splits
- BigBind/BayesBind
- LIT-PCBA duplicate/analogue audit preprint and public audit code 2025

Key conclusion: no single threshold or split earns the term leak-free; deployment regimes and continuous novelty audits are required.

## Source-code audit procedure

For each public repository, the review sought:

- root README and documented inputs;
- inference and preprocessing entry points;
- model classes and feature constructors;
- data/checkpoint download format;
- explicit ligand coordinate/conformer fields;
- contact/distance/pose targets;
- pocket extraction and centring logic;
- benchmark split files;
- environment/dependency definitions;
- root licence and transitive tool terms;
- commit/release recency; and
- known issues.

The most consequential inspections were:

- BANANA: pocket PDB + SMILES interface and separate residue/ligand message passing
- CASTER-DTA: GVP protein graph, 2D molecular GNN, residue–atom cross-attention
- PLANET: ProteinEGNN, LigandGAT, pair/contact and intraligand-distance objectives
- GRIPHIN: ligand Cartesian/Fourier input and ligand-centred construction
- T-ALPHA: explicit protein, ligand, and complex channels
- TankBind: ligand internal distances and pose optimization
- DrugCLIP: conformer-bearing molecule LMDB
- PocketDTA: GraphMVP ligand pretraining dependency
- GEMS: CleanSplit JSON, filtering code, similarity matrices, and search baseline
- PLINDER: disabled affinity query and issue #94
- gMolAI v2.0: residual GINE, per-atom states, global calibrated representation, and missing root licence
- HoloProt/MaSIF/dMaSIF: surface preprocessing and licensing burden

A repository was not cloned into blueprint/. No code or checkpoints were redistributed. Before implementation, every selected dependency must be pinned by commit and re-audited locally.

## Negative searches and unresolved items

### Graph_RG source

Targeted searches by method name, DOI, title fragments, authors, GitHub, and CASP16 references did not locate an official public source repository by 20 August 2026. This is recorded as “not found”, not “does not exist”. Contacting the authors is recommended if reproduction is required.

### PLANET v2.0

The January 2026 preprint and project site were located. Final peer-reviewed code/licence/input classification remains open and requires a later audit.

### Licences

No clear root licence was found during the reviewed pass for gMolAI v2.0, PLANET, and 3DProtDTA. Some other projects expose custom or restrictive terms. Absence of a visible licence is treated as no reuse permission, not implicit openness.

### Dataset rights

Web access, publication, or a free registration path does not establish redistribution rights. The blueprint records access observations; institutional/legal review remains required.

## Search stopping rule

A search wave was considered saturated when:

- repeated queries returned already classified direct methods;
- backward/forward citations produced encoder/dataset papers rather than new exact competitors;
- title/abstract and source-level coordinate classification agreed for all principal methods; and
- each architectural component and evaluation recommendation had at least one primary source.

Because the field is moving rapidly, this is a cut-off rather than a permanent stop. Surveillance is required before architecture freeze and manuscript submission.

## Evidence quality labels

- paper and source: primary publication plus implementation inspected
- paper and targeted source search: primary publication reviewed; code searched but not found
- preprint and source: non-peer-reviewed current evidence with implementation/deposit
- paper: primary publication only
- official site/deposit: current first-party dataset/project record
- secondary discovery only: not used as sole support for material conclusions

## Reproducibility notes

- Current web counts and licence statements are snapshots as of the access date.
- Repository default branches can change; implementation must pin commit hashes.
- Dataset counts can differ after filters and updates; reports must cite the exact manifest.
- Preprints are explicitly labelled and cannot silently replace peer-reviewed evidence.
- Direct quotations were not copied into the blueprint; findings are paraphrased.


## Final dataset and infrastructure expansion — 20 August 2026

A final category-completeness pass searched beyond the usual affinity benchmarks for receptor coordinates, predicted structures, apo/holo pairs, pocket archives, target-family resources, assay archives, mapping registries, conformational resources, and bias diagnostics.

Representative query groups included:

- CrossDocked2020 v1.3, APObind, scPDB, and PDBFlex official data/source pages
- RCSB PDB API/download/licence and AlphaFold DB release/licence
- KLIFS kinase pocket API and GPCRdb structures/refined-model documentation
- PubChem BioAssay downloads, Guide to PHARMACOLOGY 2026.2, and Drug Target Commons bulk/API access
- MUV, DEKOIS 2.0, and OpenFE/public FEP benchmark origins
- SIFTS PDB–UniProt mappings, UniProtKB downloads/licence, and CATH v4.4 classifications

Primary-page findings added to the registries:

- RCSB PDB archive coordinates and API data are CC0; entries still require dated, construct-aware processing.
- AlphaFold DB predictions are CC BY 4.0 and provide confidence/PAE, but omit ligand, metal, cofactor, water, and site context.
- CrossDocked2020 v1.3 documentation reports a 52,126,979-file full archive and exposes pose/RMSD plus mixed PDBbind-derived pK fields; only receptor/pocket views are eligible for the strict pipeline.
- APObind provides code and matched PDBbind-2019 apo/holo files, with unclear root licensing and derivative/mapping questions.
- KLIFS v3.2 was live and current at review, with standardized 85-residue kinase sites; interaction fingerprints are ligand-conditioned.
- GPCRdb supplies versioned structures, refined models, generic numbering, and activation-state information; post-cutoff refinements are a leakage channel.
- PubChem BioAssay and Guide to PHARMACOLOGY are valuable source-specific activity views, not automatic structure/site mappings.
- Drug Target Commons was verifiable through its papers and published API/download descriptions, but the live service was not reliably accessible; it remains conditional.
- MUV and DEKOIS 2.0 are bias diagnostics rather than quantitative-affinity evidence.
- PDBFlex is accessible but its displayed PDB snapshot ends in January 2017.
- OpenFE/public FEP series are scientifically useful but heavily reused and subject to the current target-mirroring debate.
- SIFTS, UniProtKB, and CATH are required mapping/OOD infrastructure; none alone proves pocket novelty.

No bulk resource was downloaded and no third-party coordinate, ligand, or checkpoint file was copied into the repository. The final evidence package contains 71 publication/resource assessments and 46 dataset/infrastructure records; these counts refer to reviewed registry rows, not a claim that future or private data sources do not exist.
