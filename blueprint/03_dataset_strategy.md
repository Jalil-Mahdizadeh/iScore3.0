# Dataset Strategy

**Snapshot date:** 20 August 2026  
**Principle:** separate structural examples, affinity measurements, pocket mappings, and evaluation benchmarks. They are related tables with provenance—not one merged CSV.

## 1. Dataset roles

iScore3.0 needs four kinds of evidence:

1. **Supervised structural affinity:** receptor structure, independently defined site, coordinate-free ligand identity, and quantitative affinity.
2. **Large-scale activity:** many target–compound measurements that can be mapped to a receptor and site with explicit uncertainty.
3. **Representation pretraining:** protein pockets and structures that need not carry reliable affinity labels.
4. **Locked evaluation:** hard, external, time-forward, or blinded targets not used for representation selection or hyperparameter tuning.

A resource can fill more than one role only through separately versioned views. Training and test records from the same source must be globally deduplicated before splitting.

## 2. Dataset inventory

Counts are release-dependent snapshots, not permanent properties. The machine-readable inventory in [evidence/datasets.tsv](evidence/datasets.tsv) records versions and access dates.

No finite list can literally contain every private collection, patent table, supplementary medicinal-chemistry series, or future release. This inventory therefore aims for **category completeness** across major accessible resources relevant to labels, receptor structures, sites, mappings, conformations, and evaluation. A resource must be added to the registry and pass the same provenance/licence audit before use; omission from this document is not implicit approval.

| Resource | Approximate reviewed scope | Affinity/activity | Protein 3D and site | Access/licence issue | Recommended role |
|---|---|---|---|---|---|
| PDBbind+ v2025 | 29,001 protein–ligand complexes on current site | curated Kd/Ki/IC50 and related values | co-complex structures | commercial/subscription; v2020R1 advertised free after registration | licensed structural core |
| PDBbind v2020/CleanSplit | about 19k general-set complexes before cleaning | quantitative | co-complex | structure rights from PDBbind; split metadata via GEMS | legacy core plus leakage audit |
| LP-PDBbind | 11,513/2,422/4,860 train/validation/test in published downstream use | PDBbind labels | co-complex | requires underlying PDBbind access | hard generalization evaluation |
| PDBbind-Opt / HiQBind | quality-optimized structural affinity collections | quantitative | co-complex | verify each release and derivative right | high-quality subset/external check |
| BioLiP2-Opt | more than 3k PDB entries and more than 5.6k complexes in deposited release | mapped/curated | co-complex | deposited CC BY 4.0; inspect source lineage | open structural complement |
| PLINDER 2024-06/v2 | >400k systems, >50k unique small molecules, extensive annotations | affinity field currently unsafe | holo plus linked apo/predicted structures | curated data Apache-2.0; upstream terms propagate; hundreds of GB | pocket/structure pretraining and conformation tests |
| BigBind | 582,957 activities, 399,090 compounds, 1,107 pockets in paper | activity/pChEMBL-like data | assays mapped to CrossDocked pockets | MIT code; source-data terms propagate | scalable pose-free activity training |
| BindingDB 202608 | 3,239,327 measurements, 1,433,271 compounds, 11,503 targets | Kd/Ki/IC50 and others, qualifiers | target sequences and some PDB links; sites require mapping | curator data CC BY 4.0; imported data retain source terms | large affinity source |
| ChEMBL 36 | large curated bioactivity release, July 2025 | normalized activities plus assay metadata | targets; structures/sites require mapping | CC BY-SA 3.0 | large activity source |
| Papyrus | standardized multi-source bioactivities | continuous/categorical endpoints | protein identifiers/sequences; no pocket structure mapping | version-specific open deposit | standardization and cross-source view |
| BioLiP2 | 385,160 protein chains and 50,064 interactions with affinity in paper snapshot | affinities imported from several sources | biologically filtered sites | weekly updates; overlapping provenance; BSD curation code | site annotations and open supplement |
| Binding MOAD | >32k curated complexes and 12,098 with binding data in 2019 update | quantitative/qualitative | co-complex | legacy service/access stability; verify terms | independent curation audit |
| ATOM3D LBA | PDBbind-derived identity-split task | PDBbind labels | co-complex | underlying PDBbind terms remain | preprocessing/baseline convenience |
| MISATO | roughly 20k PDBbind systems with QM/MD features | PDBbind-derived | dynamic complexes | derivative and source terms; ligand 3D | optional protein pretraining/non-core analysis |
| Davis | dense 68-ligand × 442-kinase matrix | 30,056 Kd values | target sequences; structures must be mapped | verify mirror/source rights | narrow baseline and within-target analysis |
| KIBA | 118,254 drug–target pairs, 2,111 drugs, 229 targets | aggregated KIBA score | structures require mapping | not a thermodynamic affinity | literature comparability only |
| Metz | kinase inhibitor affinity collection | quantitative | structures require mapping | verify original distribution | kinase baseline |
| CASF-2016 | 285 high-quality complexes in 57 target clusters | binding constants | native complex/pocket | PDBbind benchmark terms | legacy scoring/ranking comparability |
| CASP16 affinity category | 140 measurements across five systems | blinded pharmaceutical affinities | target structures/pockets | Zenodo benchmark terms | locked external ranking |
| BayesBind | hard targets dissimilar to BigBind training | binders/putative nonbinders with Bayesian enrichment | pockets | MIT BigBind repository; inspect component terms | hard virtual-screening evaluation |
| LIT-PCBA | 15 targets; paper reports 7,844 actives and 407,381 confirmed inactives | experimental binary activity | reference structures/sites | open site; pin archive; a 2025 preprint reports exact and analogue leakage | diagnostic or independently rebuilt split only |
| DUD-E | 102 targets, 22,886 clustered actives, 50 property-matched decoys per active | active/constructed decoy | reference pocket | dataset terms; known decoy bias | diagnostic only |
| D3R / SAMPL | small blinded challenge series | affinity/ranking and pose by challenge | target-specific structures | challenge-specific terms | external case studies |
| Merck FEP benchmark | 264 ligands over eight targets in common release | relative free energies | target structures | verify redistribution and experimental source | within-target ranking stress test |
| RCSB PDB / wwPDB | continuously updated archive | no uniform label | experimental structures and bound components | coordinate/API data CC0; freeze dates and obsolete entries | primary receptor source |
| AlphaFold DB | proteome-scale predicted structures | no label | predicted receptor with pLDDT/PAE; no ligand site | CC BY 4.0 | predicted-structure robustness |
| APObind | PDBbind-2019-derived matched collection | inherited PDBbind labels | matched apo/holo receptors | derivative rights and mapping need audit | conformation stress test |
| CrossDocked2020 v1.3 | very large redocked/cross-docked pose archive | RMSD and mixed PDBbind-derived pK fields | receptors and pocket groupings | component terms; >52 million files; contains ligand poses | receptor/pocket-only pretraining |
| scPDB | legacy curated druggable-site collection | no uniform primary label | PDB-derived sites | verify release/access/reuse terms | conditional site resource |
| KLIFS v3.2 | 326 kinases, 6,738 PDB structures on reviewed site | linked kinase bioactivities | standardized 85-residue pockets, states, IFPs | verify bulk and upstream terms | kinase site alignment/strata |
| GPCRdb | evolving GPCR structure/model resource | linked bioactivities | structures, refined models, generic numbering | source terms and post-release refinements | GPCR site/state robustness |
| PubChem BioAssay | very large contributor-deposited assay archive | outcomes, concentrations, heterogeneous readouts | targets; sites require mapping | contributor licences and assay heterogeneity | assay-specific pretraining/rebuilt screens |
| Guide to PHARMACOLOGY 2026.2 | curated ligand-target interactions | quantitative affinity/potency | targets; sites require mapping | ODbL plus CC BY-SA 4.0 | high-quality activity complement |
| Drug Target Commons 2.0 | ~6.5m protein bioactivities in last verifiable paper | quantitative, assay-annotated | targets; sites require mapping | CC BY-NC-SA 3.0; current service status must be rechecked | conditional audit source |
| MUV | 17 PubChem-derived target sets | active/inactive | receptor sites require mapping | verify original lineage | ligand-bias diagnostic |
| DEKOIS 2.0 | 81 target sets in 11 classes | actives/designed decoys | PDB targets plus pose files | mirror licence unclear; constructed decoys | diagnostic only |
| PDBFlex | legacy same-protein conformation clusters | no label | experimental conformational variation | archive uses PDB through January 2017 | optional flexibility annotation |
| OpenFE/public FEP benchmark | multi-target medicinal-chemistry series | relative/series affinity | prepared receptors | heavily reused; exact subset and lineage required | leakage-aware ranking stress test |
| SIFTS / UniProtKB | weekly/eight-week mapping registries | no label | residue mappings and versioned target sequences | public; pin exact releases | mandatory identity/mapping infrastructure |
| CATH / Gene3D v4.4 | 601,493 PDB domains and 6,573 superfamilies reported | no label | domain/fold classification | CC BY 4.0 | OOD strata, not a pocket-similarity substitute |


## 3. Structural supervised core

### 3.1 PDBbind+

The [PDBbind+ site](https://pdbbind-plus.org.cn/) reports that version 2025 was released on 8 February 2026 and contains 29,001 protein–ligand complexes. It is a commercial/subscription release. The site also advertises PDBbind v2020R1 for free download after registration. “Free download” is not equivalent to unrestricted redistribution.

PDBbind remains the natural structural-affinity core because every entry links a measured label to a co-crystal receptor and known site. Its disadvantages are equally important:

- crystallographic selection favours strong, tractable binders;
- holo structures encode ligand-induced conformation;
- protein, pocket, scaffold, series, and publication redundancy are extensive;
- Kd, Ki, IC50, and other labels are mixed;
- common benchmark cores overlap many training releases; and
- licences can restrict derived structure distribution.

**Decision:** use only after institutional acceptance of the exact release agreement. Store scripts, IDs, hashes, and derived non-coordinate metadata in version control; keep restricted structures in an access-controlled external data root.

### 3.2 CleanSplit/GEMS

The [GEMS project](https://github.com/camlab-ethz/GEMS) and [Zenodo artifacts](https://doi.org/10.5281/zenodo.14260170) expose PDBbind redundancy analysis and CleanSplit metadata. This is a required evaluation view, not an optional paper comparison. Use the authors’ exact split manifest and independently verify that no post-release preprocessing reconnects train and test through alternate PDB structures or standardized ligand duplicates.

Do not quote a sample count without naming the GEMS artifact and filter configuration; published downstream pipelines apply different quality/redundancy filters.

### 3.3 LP-PDBbind

[LP-PDBbind](https://github.com/THGLab/LP-PDBBind) was designed to reduce protein and ligand overlap. A commonly reported split contains 11,513 training, 2,422 validation, and 4,860 test examples. It is useful as a recognized hard benchmark, but the structures still originate from PDBbind. Run the official split unchanged and also run the project’s own union-component split; neither substitutes for the other.

### 3.4 PDBbind-Opt, HiQBind, and BioLiP2-Opt

[PDBbind-Opt](https://arxiv.org/abs/2411.01223) applies structure/label quality optimization and motivates HiQBind/BioLiP2-Opt. [BioLiP2-Opt on Figshare](https://figshare.com/articles/dataset/BioLiP2-Opt_Dataset/27430305/1) is an attractive open complement and is labelled CC BY 4.0 at the deposit. Before use, verify the final journal version, exact record-level source lineage, unit conversions, and whether all structures/ligands fall within the proposed chemical scope.

Quality-optimized sets can be cleaner but smaller and more selected. Compare performance and calibration with and without them; do not assume that curation monotonically improves external generalization.

## 4. Large-scale activity and affinity sources

### 4.1 BindingDB

The [BindingDB 202608 download](https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp) lists 3,239,327 measurements, 1,433,271 compounds, and 11,503 targets. The home page distinguishes approximately 1.6 million measurements curated directly by BindingDB and now provides publication and curation dates, which are valuable for temporal splitting.

Use the monthly 2D SDF/TSV and assay mapping; do not download the 3D ligand SDF for the strict pipeline. Preserve relation qualifiers and all endpoint-specific fields. Direct BindingDB-curated data are CC BY 4.0, while imported ChEMBL and other subsets retain their own licences.

PDB links in an activity database are not automatically assay-to-structure matches. A reported ligand may occur in a PDB with another target or construct. The PLINDER issue described below is a concrete warning.

### 4.2 ChEMBL

[ChEMBL release 36](https://chembl.gitbook.io/chembl-interface-documentation/downloads) was released in July 2025 with DOI [10.6019/CHEMBL.database.36](https://doi.org/10.6019/CHEMBL.database.36). ChEMBL provides target, assay, document, confidence, units, relations, and standardized activity fields under CC BY-SA 3.0.

Recommended inclusion begins with direct single-protein binding assays at high target-confidence scores and explicit molar endpoints. Retain lower-confidence or functional assays only as separately flagged sensitivity sets. Patent and journal records need an earliest-public-date policy.

### 4.3 Papyrus

[Papyrus](https://pmc.ncbi.nlm.nih.gov/articles/PMC9824924/) provides standardized, quality-stratified bioactivities assembled from multiple sources. It can accelerate a reproducible first pass and helps compare normalization decisions. Because it overlaps BindingDB, ChEMBL, and other databases, it must be treated as an alternate view, not additional independent observations. Pin a specific [Zenodo release](https://zenodo.org/records/7019874) and trace every record to its origin.

### 4.4 BigBind

The [BigBind paper](https://doi.org/10.1021/acs.jcim.3c01211) reports 582,957 activities for 399,090 compounds across 1,107 pockets. Its central innovation is mapping activity records to structurally clustered pockets without requiring an experimental structure for every compound. It is the most directly reusable data design for pose-free screening.

Use BigBind first for binary/weakly supervised pretraining and BANANA reproduction. For quantitative affinity, recover original endpoint, units, qualifier, and assay context rather than treating all pChEMBL-like values as identical. Recompute overlap against every locked evaluation ligand, pocket, target, and publication.

### 4.5 PubChem BioAssay, Guide to PHARMACOLOGY, and Drug Target Commons

[PubChem BioAssay](https://pubchem.ncbi.nlm.nih.gov/docs/bioassays) provides assay-level outcomes, target records, and downloadable SMILES at enormous scale. It is most useful for reconstructing individual confirmatory assays, counterscreen-aware pretraining, and independently rebuilt screening benchmarks. Pooling AIDs without assay lineage, target type, readout direction, and contributor licence would create a heterogeneous label soup.

The [IUPHAR/BPS Guide to PHARMACOLOGY 2026.2 download](https://www.guidetopharmacology.org/download.jsp) supplies manually curated ligand-target interactions, target/UniProt mappings, SMILES, and references under an ODbL/CC BY-SA framework. Its selectivity and documented provenance make it a high-quality complement, although affinity, potency, species, and assay contexts remain separate endpoints.

[Drug Target Commons 2.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC6146131/) reported about 6.5 million protein-target bioactivities with detailed assay annotation and a noncommercial share-alike licence. Because the live service was not reliably accessible during this review and the published release is old, it is conditional: verify availability, version, licence, and overlap before any pipeline dependency.

**Decision:** add these as source-specific views after global deduplication. None supplies a trustworthy pocket automatically; all require the S0-S4 mapping process.

## 5. Structure and pocket resources

### 5.1 PLINDER

[PLINDER](https://github.com/plinder-org/plinder) contains more than 400,000 protein–ligand interaction systems, more than 50,000 unique small molecules, hundreds of annotations, rich protein/pocket/ligand similarity information, and links from holo systems to apo and AlphaFold structures. The reviewed current full release is 2024-06/v2 and requires hundreds of gigabytes.

It is excellent for:

- self-supervised pocket pretraining;
- structure quality and site-definition experiments;
- holo/apo/predicted matched evaluation;
- pocket clustering and leakage analysis; and
- realistic pocket-detector benchmarks.

It is **not currently approved as an affinity-label source**. The repository disables ligand_binding_affinity queries, and [open issue #94](https://github.com/plinder-org/plinder/issues/94) explains that BindingDB PDB identifiers were matched without adequate sequence/assay validation, producing incorrect associations. The field may be used only after an upstream fixed release or an independently audited remapping.

PLINDER-curated data are advertised under Apache-2.0, its software repository uses GPL-2.0, and imported records retain source licences. A release-level licence manifest is required.

### 5.2 BioLiP2

The [BioLiP2 paper](https://doi.org/10.1093/nar/gkad630) reports 385,160 protein chains, 781,684 protein–ligand interactions, and 50,064 interactions with affinity data at its snapshot. It filters crystallization additives and exposes binding residues, SMILES, and non-redundant views. Its curation code is BSD-licensed.

The affinity total combines MOAD, PDBbind-CN, BindingDB, and manual records, so source counts overlap and are not independent labels. Use BioLiP2 primarily for biological-site annotation and source cross-checking. Pin its weekly release and retain the original affinity source.

### 5.3 Binding MOAD

The [2019 Binding MOAD update](https://pmc.ncbi.nlm.nih.gov/articles/PMC6589129/) described more than 32,000 curated complexes and 12,098 with binding data. It offers valuable independent manual curation and ligand validity information. Because current bulk access and service stability can change, confirm availability, rights, and identifiers before making it a dependency.

### 5.4 ATOM3D and MISATO

[ATOM3D](https://www.atom3d.ai/) supplies convenient PDBbind-derived ligand-binding-affinity tasks and protein-identity splits. It can validate geometric code, but it does not remove the underlying licence or leakage issues.

[MISATO](https://doi.org/10.1038/s43588-024-00627-2) adds quantum-mechanical and molecular-dynamics information to approximately 20,000 PDBbind complexes. Its ligand trajectories violate the strict ligand pathway and it should not provide affinity features. Receptor-only dynamic pretraining is an optional later study with a careful derivative-rights audit.

### 5.5 Primary experimental and predicted receptor sources

The [RCSB PDB/wwPDB archive](https://www.rcsb.org/docs/programmatic-access/file-download-services) is the canonical coordinate source. PDB archive coordinates and RCSB API data are CC0, but entries still require assembly selection, obsolete/superseding-entry tracking, construct parsing, alternate-location policy, and a freeze date. A PDB ligand or affinity annotation is provenance evidence, not automatic proof that a separately sourced assay measured the same construct/site.

[AlphaFold DB](https://alphafold.ebi.ac.uk/faq) provides predicted protein structures under CC BY 4.0 with pLDDT and PAE. It omits bound ligands, metals, cofactors, ions, and waters, making it a realistic receptor-only input but not a site oracle. Freeze the model version and relevant template/database dates; mask low-confidence geometry; run predicted-structure results as a distinct robustness stratum.

**Decision:** RCSB structures are the primary experimental receptor source. AlphaFold DB is an explicit fallback/evaluation mode, never a replacement label source.

### 5.6 Matched conformations and large pose archives

[APObind](https://github.com/devalab/Apobind) links apo structures to PDBbind-2019 holo complexes and is directly relevant to holo privilege. Source review found that its small construction toolkit is not turnkey: scripts contain hard-coded scratch paths and local BLAST/PDB database locations, invoke BLAST through the shell, depend on PyMOL/ProDy/RDKit, and use co-crystal ligand coordinates to select/measure the reference site. No clear repository licence was observed. Use only the released receptor pairs after independent sequence/construct/site remapping and PDBbind-rights review; do not import its ligand-conditioned pocket logic into receptor-only Mode B.

[CrossDocked2020 v1.3](https://github.com/gnina/models/blob/master/data/CrossDocked2020/README.md) contains redocked/cross-docked receptor-ligand poses, pocket groupings, RMSD labels, and a PDBbind-derived pK field that does not distinguish Kd, Ki, and IC50. The full archive contains more than 52 million files. For strict iScore3.0, only a receptor/pocket-derived view may be used; ligand coordinates, pose quality, and mixed pK fields cannot enter the affinity path.

[PDBFlex](https://pdbflex.org/) clusters highly similar PDB chains and exposes global/local conformational variation, but its reviewed archive is based on PDB through January 2017 and is not pocket-specific. It can supply a secondary flexibility annotation, not a current benchmark.

**Decision:** use APObind/PLINDER matched conformations first; consider PDBFlex for annotation. Treat CrossDocked2020 as receptor/site pretraining material only unless a separate pose task is deliberately approved.

### 5.7 Family-specific and legacy site resources

[KLIFS v3.2](https://klifs.net/) standardizes an 85-residue kinase binding site and exposes structures, conformations, subpockets, and interaction fingerprints. Its family alignment is valuable for kinase pocket novelty and conformation strata. Ligand interaction fingerprints and co-crystal-derived annotations are privileged, so strict receptor-only pocket detection must not consume them.

[GPCRdb](https://docs.gpcrdb.org/structures.html) supplies GPCR structures, generic residue numbering, activation-state annotations, refined structures, and homology models. Freeze the database/model version: later refinements can incorporate structures published after an evaluation cutoff. Use its generic positions and state labels for family-specific robustness, while isolating ligand-interaction annotations.

[scPDB](https://drugdesign.unistra.fr/scPDB/) is a useful legacy collection of annotated druggable PDB sites. Current bulk availability, processing version, licence, and overlap must be verified before it becomes a dependency.

**Decision:** these are specialist annotations and evaluation strata. They do not replace a general receptor-only site detector or confer independent affinity evidence.

### 5.8 Mapping and structural-novelty infrastructure

[SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/) releases weekly residue-level mappings among PDB, UniProt, and domain/function resources. [UniProtKB](https://www.uniprot.org/help/downloads) supplies versioned canonical/isoform sequences under CC BY 4.0. Both are mandatory inputs to construct-aware mapping, but neither proves that an assay target equals a PDB construct.

[CATH/Gene3D v4.4](https://www.cathdb.info/) supplies domain and superfamily classifications under CC BY 4.0. It supports fold-level OOD strata and diagnostics. CATH separation is coarser than local binding-site novelty and must accompany—not substitute for—continuous pocket similarity.

**Decision:** pin SIFTS, UniProt, and CATH snapshots in every data release and retain their identifiers in the provenance graph.

## 6. Benchmark resources

### 6.1 CASF-2016

[CASF-2016](https://doi.org/10.1021/acs.jcim.8b00545) contains 285 high-quality complexes across 57 target clusters and evaluates scoring, ranking, docking, and screening power. iScore3.0 can run scoring and ranking without a query pose by using the defined receptor pocket and SMILES; docking power is outside scope.

CASF is deeply represented in PDBbind training and hyperparameter history. It is retained for literature comparability only. A model chosen on CASF is not externally validated.

### 6.2 CASP16 pharmaceutical affinity

The [CASP16 affinity deposit](https://zenodo.org/records/16762332) provides blinded pharmaceutical affinity series used in the official assessment. The set’s 140 measurements over five systems are small but unusually valuable because methods could not tune to the labels before the challenge. Preserve the original challenge protocol, receptor/site inputs, relation qualifiers, and target-level ranking metrics.

### 6.3 BayesBind

[BayesBind](https://pmc.ncbi.nlm.nih.gov/articles/PMC10980085/) selects targets structurally dissimilar to BigBind training and introduces a Bayesian enrichment statistic that accounts for uncertain active fractions. The authors found evaluated methods did not materially exceed a nearest-neighbour baseline. It is a high-value test of whether iScore3.0 learns pocket-conditioned recognition beyond analogue retrieval.

### 6.4 LIT-PCBA

The [LIT-PCBA publication](https://pubmed.ncbi.nlm.nih.gov/32282202/) reports 15 targets, 7,844 true actives, and 407,381 experimentally confirmed inactives. However, the reproducible 2025 preprint [Data Leakage and Redundancy in the LIT-PCBA Benchmark](https://arxiv.org/abs/2507.21404) reports 2,491 unique inactives shared between its training and validation sets, three query compounds crossing into train/validation, and extensive analogue redundancy; its [audit code](https://github.com/sievestack/LIT-PCBA-audit) is public. Because this evidence is currently a preprint, iScore3.0 must reproduce the audit on a hash-pinned official archive. Unmodified LIT-PCBA is not confirmatory evidence. Use it as a shortcut diagnostic or rebuild a new split from the underlying assay records with exact-identity, analogue, provenance, and query-set isolation.

### 6.5 DUD-E

[DUD-E](https://doi.org/10.1021/jm300687e) contains 102 targets and 22,886 clustered actives, with 50 property-matched computational decoys per active. The original paper explicitly cautions that the construction is inappropriate for evaluating ligand-based 2D methods. A SMILES model can exploit decoy-selection artifacts.

Use DUD-E only to quantify bias: compare gMolAI-only, pocket-only, and full models; use asymmetric target/ligand resampling; never make it the headline benchmark.

### 6.6 D3R, SAMPL, and target-series data

D3R/SAMPL challenge series, the Merck FEP benchmark, and carefully curated medicinal-chemistry series probe within-target ranking and uncertainty. They are small and may emphasize relative rather than absolute affinity. Each needs its own endpoint, experimental covariance, receptor-conformation, and licence audit.

### 6.7 MUV, DEKOIS 2.0, and public FEP/OpenFE series

[MUV](https://doi.org/10.1021/ci8002649) comprises 17 PubChem-derived virtual-screening sets designed to reduce ligand-space benchmark bias. It is useful for testing whether a SMILES model still exploits obvious chemical shortcuts, but it has binary labels, no native general pocket package, and no claim of protein/pocket novelty.

[DEKOIS 2.0](https://doi.org/10.5281/zenodo.8131256) supplies 81 target sets with actives, designed decoys, PDB receptors, and pose files in the reviewed mirror. Strip ligand coordinates for strict experiments. Treat it as a construction-bias diagnostic; the mirror showed no clear licence during review, so original rights must be resolved.

The [OpenFE public benchmark](https://industrybenchmarks2024.readthedocs.io/en/latest/public/overview.html) and associated public binding-free-energy series offer realistic within-target ranking. They are now heavily reused, and a 2026 preprint specifically raises target-mirroring concerns for such series. Pin the exact subset, trace each experimental value, compare against the entire pretraining inventory, and never call it blind if its targets or labels influenced development.

**Decision:** these resources broaden diagnostics; none is primary proof of absolute, cross-target, dual-novel affinity prediction.

## 7. Measurement schema

The canonical observation table should include at least:

| Group | Required fields |
|---|---|
| identity | observation_id, source_record_id, source_version, source_url/DOI |
| ligand | original_smiles, standardized_smiles, InChIKey, connectivity_key, scaffold_id, standardization_version |
| target | source_target_id, UniProt accession/version, organism, construct sequence, mutations, complex/subunit context |
| assay | assay_id, assay_type, description, biochemical/cellular flag, temperature, pH, substrate/cofactor conditions |
| endpoint | Kd/Ki/IC50/KIBA/other, relation, raw value, raw unit, molar interval, transformed interval |
| structure | PDB/mmCIF or predicted-structure ID, chain mapping, structure date, method, resolution/confidence |
| site | site_id, definition method, residue set, detector/version/rank, mapping confidence |
| provenance | publication/patent, earliest public date, curation date, parent database, licence lineage |
| quality | parsing status, target confidence, unit confidence, replicate group, anomaly flags |
| split | protein cluster, pocket cluster, scaffold cluster, ligand-similarity component, temporal eligibility, final split hash |

Structures, ligands, measurements, assays, sites, and mappings should be normalized tables. A many-to-many association is expected; flattening early makes provenance and leakage errors difficult to detect.

## 8. Endpoint harmonization

### 8.1 Unit conversion

Convert only recognized concentration units to mol/L. For an exact molar value x, compute negative log10(x). Preserve the raw value and the converted interval. A qualifier maps to an interval before transformation; note that the inequality reverses under negative logarithm.

Reject dimensionless, ambiguous, zero, negative, or implausible concentrations to quarantine. Do not infer units from typical ranges.

### 8.2 Endpoint separation

- pKd and pKi are primary quantitative heads.
- pIC50 is separate and receives available assay-context features.
- KIBA remains a KIBA score.
- EC50 and cellular potency are separate exploratory tasks.
- binary activity has its own threshold/source definition.
- docking or computed energies are never experimental affinity labels.

A hierarchical multitask model may share representations while retaining endpoint-specific offsets/noise. Report each endpoint alone and test whether pooled training helps the primary endpoints.

### 8.3 Censoring and replicates

Represent “<”, “>”, and ranges as intervals and optimize an interval-censored likelihood. Keep individual replicate observations linked to a replicate group. Recommended derived views are:

1. all records with source/assay random effects;
2. robust median/consensus within a tightly defined assay group;
3. high-agreement subset; and
4. conflict set used to evaluate aleatoric uncertainty.

Never average Kd, Ki, and IC50 into one number.

## 9. Target-to-structure and target-to-site mapping

This is the highest-risk step when scaling beyond structural databases. Use confidence tiers:

| Tier | Evidence | Allowed use |
|---|---|---|
| S0 | affinity and co-complex are the same curated record/construct | structural supervised core and evaluation |
| S1 | exact construct/sequence and site mapped to another structure of the same target | training; external test if mapping frozen |
| S2 | exact UniProt target, compatible domain/construct, and a well-supported canonical site | large-scale pretraining/training with mapping token |
| S3 | homolog or predicted structure and inferred site | robustness/pretraining only |
| S4 | family/name match, ambiguous chain, multi-target assay, or multiple plausible sites | quarantine |

Checks include sequence alignment over the pocket domain, mutations, isoform, organism, biological assembly, cofactors, domain boundaries, and whether the assay mechanism corresponds to the selected site. Multi-site targets need explicit site assignment or a multiple-instance model; choosing the best pocket using the label is prohibited.

Each mapping produces a confidence tier and reason codes used as model features or stratification—not silently discarded metadata.

## 10. Deduplication

Deduplicate before splitting and across every source:

1. standardize molecules and link exact connectivity/stereochemical identities;
2. normalize target accessions and construct sequences;
3. group exact and near-identical receptor pockets across PDB entries;
4. link measurements by publication, assay, endpoint, target, and ligand;
5. identify records imported from one database into another;
6. preserve conflicting records while preventing one measurement from being counted repeatedly; and
7. compare all future external tests against the complete training/pretraining inventory.

A co-crystal and a BindingDB record can be the same experiment. Source multiplicity is provenance, not replication.

## 11. Leakage-resistant splitting

### 11.1 Unit of independence

Build a graph whose nodes are receptor-site/ligand examples. Add an edge when either side is too similar:

- ligand: exact connectivity, shared scaffold, or fingerprint similarity above a frozen threshold;
- protein: high sequence/domain identity;
- pocket: structural or residue-environment similarity above a frozen threshold;
- provenance: same assay series/publication when appropriate.

Connected components of this union graph are assigned wholly to train, validation, or test. This prevents transitive leakage—for example A similar to B through ligand chemistry and B similar to C through pocket structure.

### 11.2 Required split views

- random pair split, reported only as a diagnostic;
- ligand-scaffold novel;
- protein-family/pocket novel;
- dual novel in both ligand and pocket;
- temporal holdout by earliest public measurement/structure date;
- supplied-site versus receptor-only detected-site;
- holo versus matched apo versus predicted structure; and
- target-series holdout for ranking.

Thresholds are chosen using domain reasoning and sensitivity analysis before model comparisons. Report exact identities, maximum ligand similarity, protein identity, pocket similarity, and date gap to training for every test record.

### 11.3 Pretraining leakage

Unlabelled structural or molecular pretraining can contain evaluation entities. That is not automatically label leakage, but it can enable identity memorization. Report whether each test ligand, pocket, or target occurred in gMolAI, protein-language-model, PLINDER, or other pretraining corpora where knowable. Run exclusion or nearest-neighbour sensitivity experiments for headline sets.

## 12. Quality-control gates

A data release is eligible for model training only when automated reports show:

- 100% source/version/licence lineage for included records;
- 100% valid endpoint/unit/relation conversion;
- no query-ligand coordinates in strict ligand artifacts;
- no test component connected to train at frozen thresholds;
- no duplicated measurement assigned to different splits;
- audited target/construct/site mappings with tier distribution;
- pocket extraction reproducibility from structure plus site manifest;
- explicit counts for exclusions and failure reasons;
- visual review of a stratified sample including metals, alternate locations, missing loops, and multiple chains; and
- immutable SHA-256 manifests for raw references, processed records, features, and splits.

Any failed gate blocks the corresponding experiment; it does not justify silently dropping records after seeing model performance.

## 13. Recommended staged data build

### Stage D0 — open feasibility set

Pin RCSB PDB, UniProtKB, SIFTS, and CATH snapshots; use a rights-cleared structural subset such as BioLiP2-Opt, plus Davis/KIBA only for architecture smoke tests. Build the schema, strict feature tests, split code, and ligand/protein-only baselines. Results are not publication claims.

**Gate-0/1 result:** a provenance-first RCSB/BindingDB 202608 Kd pilot retained 61 exact uncensored measurements across five exact-construct components, plus five label-quarantined historical site references. BindingDB reconciliation removed 12 RCSB candidates whose omitted `>` qualifiers made censored lower bounds appear exact; a seven-row WDR5 group then failed the minimum group size. This incident makes source-level endpoint/relation reconciliation mandatory for every scaled ingestion. The strict-v2 pilot passed bounded data construction but is too small for external claims, lacks a validated local structure-similarity edge, and remains holo-privileged. See [the pilot data card](../reports/gate01/pilot_data_card.md).

### Stage D1 — licensed structural core

After approval, acquire and pin PDBbind. Generate supplied-site pockets, CleanSplit, LP-PDBbind, union-component splits, and PLINDER/APObind/AlphaFold-linked holo/apo/predicted views. CrossDocked2020 may contribute only receptor/pocket artifacts. Freeze test labels and manifests before tuning the full model.

### Stage D2 — scalable mapped activity

Ingest a pinned BindingDB or ChEMBL release, or use BigBind, with S0-S3 mapping tiers and endpoint-specific views. Add Guide to PHARMACOLOGY, PubChem BioAssay, Papyrus, or a revalidated Drug Target Commons release only as provenance-linked source-specific views. Pretrain or multitask on lower-confidence tiers; fine-tune on S0/S1 quantitative data.

### Stage D3 — external evaluation

Run CASP16, BayesBind, and approved target series exactly once for final model families. LIT-PCBA may join only if a pinned archive is independently audited and rebuilt into a valid split; otherwise it remains a red-team diagnostic. MUV, DEKOIS, DUD-E, CASF, and public FEP/OpenFE series are labelled diagnostics or secondary ranking tests. A separate evaluation owner or access-controlled script should protect labels during development.

## 14. Storage decision

Do not place raw third-party structures or bulk databases under blueprint/ or Git. The implementation repository should contain:

- download/request instructions;
- source URLs, access agreements, version IDs, and checksums;
- row-level provenance and exclusion manifests where permitted;
- deterministic processing code;
- immutable split identifiers; and
- small synthetic fixtures.

Restricted raw data live under a configurable external root with read permissions. Derived artifacts inherit source restrictions until legal review says otherwise.

## 15. Dataset decision

The recommended combination is:

- **PDBbind plus CleanSplit/LP-PDBbind** for supervised structure-linked affinity;
- **BindingDB/ChEMBL or BigBind** for scale;
- **RCSB PDB plus UniProtKB/SIFTS** for versioned receptor identity and residue mapping;
- **PLINDER, BioLiP2, APObind, AlphaFold DB, and specialist KLIFS/GPCRdb views** for pocket, conformation, family, and predicted-structure resources;
- **BioLiP2-Opt/HiQBind** for quality/open sensitivity;
- **CASP16 and BayesBind**, plus genuinely locked target series, for external evaluation; use LIT-PCBA only after an independently audited rebuild; and
- **MUV, DEKOIS 2.0, DUD-E, CASF, and public FEP/OpenFE series** only as bias, legacy, or secondary ranking diagnostics.

The most valuable dataset contribution may be the auditable target-to-pocket mapping and dual-novel split, not simply another aggregation of public measurements.
