# Technical Blueprint

## 1. Design objective

Build a scoring system in which all three-dimensional information comes from the receptor and all query-ligand information comes from SMILES-derived 1D/2D representations. The model should learn whether a ligand’s topology and chemistry are compatible with a pocket’s shape and chemical environment without inventing or inspecting a binding pose.

The first implementation must be deliberately simpler than the eventual research vision. It should make each source of improvement attributable.

## 2. System flow

1. Validate and standardize receptor structure.
2. resolve an independently supplied site or run a receptor-only pocket detector;
3. extract and cache residue- and atom-level pocket graphs;
4. standardize SMILES without conformer generation;
5. obtain gMolAI atom tokens and calibrated global ligand representation;
6. compute pose-free pocket–ligand compatibility;
7. predict an endpoint-specific affinity distribution;
8. calibrate uncertainty and attach applicability diagnostics; and
9. return provenance, warnings, and model/data versions.

Editable diagrams are in figures/system_architecture.mmd and figures/data_and_evaluation_flow.mmd.

## 3. Architectural invariants

The strict core must satisfy the following automated invariants:

- ligand objects have no coordinate field;
- receptor features are computable before the query ligand is known;
- no protein–ligand distance, angle, contact, or pose label is present;
- the pocket crop is independent of query-ligand coordinates;
- rotating/translating the receptor does not change scalar predictions beyond tolerance;
- permuting equivalent atom/residue order does not change predictions;
- receptor embeddings can be cached and reused across ligands;
- ligand embeddings can be cached and reused across receptors; and
- every feature tensor names its schema and producer version.

A separate privileged-information experiment may violate the third invariant only under a different model name, artifact tree, and manuscript label.

## 4. Receptor preprocessing

### 4.1 Input normalization

Accept mmCIF preferentially and PDB for compatibility. The preprocessing manifest records file checksum, assembly/model, chains, alternate locations, insertions, nonstandard residues, experimental method, resolution, B-factors, predicted confidence, repairs, and software versions.

Recommended initial policies:

- use one explicitly selected biological assembly/model;
- keep the highest-occupancy alternate location with deterministic tie-breaking;
- do not invent missing loops in the baseline;
- optionally repair ordinary side-chain atoms in a named sensitivity track;
- preserve metals and required cofactors with type masks;
- exclude bulk waters initially, then test high-occupancy conserved waters separately;
- reject or flag covalently bound query ligands;
- retain coordinate-origin masks for observed, repaired, and predicted atoms; and
- preserve chain/residue identifiers in every derived node.

Protein protonation can materially alter donor/acceptor and electrostatic features. The baseline should use robust heavy-atom features that do not pretend hydrogens are experimentally resolved. A named protonation pipeline may be added as an ablation.

### 4.2 Site extraction

A site manifest supplies residue IDs or a centre/radius. The extractor creates:

- a core set exactly identified by the manifest;
- one or more fixed context shells;
- an atom graph over selected residues/cofactors;
- a mapping between atom and residue nodes; and
- boundary flags marking truncated residues/edges.

Recommended starting radii are selected globally on validation data and frozen. Report sensitivity to pocket size. When using a holo-derived annotation, remove all query-ligand records before feature construction.

For Mode B, P2Rank and fpocket produce top-k site manifests through the same interface. Detector score/rank is metadata; oracle choice after affinity prediction is not allowed.

## 5. Multi-scale pocket representation

### 5.1 Residue graph

Let residue i carry scalar state s_i and vector state v_i.

**Scalar node features**

- amino-acid identity and modified/unknown token;
- chain, terminus, insertion, missing-atom, repair, and pocket-boundary flags;
- optional frozen protein-language-model embedding projected to a small dimension;
- DSSP-like secondary structure and backbone torsions where defined;
- solvent-accessible surface area and relative depth;
- experimental B-factor or predicted-confidence features with source indicator;
- coarse side-chain chemistry: charge, donor/acceptor, aromatic, hydrophobic;
- proximity to metal/cofactor and residue conservation when legitimately available.

**Vector/geometric node features**

- backbone local frame from N, C-alpha, and C when observed;
- normalized backbone directions;
- side-chain direction such as C-alpha to C-beta, with an explicit missing mask; and
- optional surface normal derived solely from receptor geometry.

**Edges**

Use both sequential and spatial relations. Candidate spatial edges are radius or k-nearest neighbours over representative coordinates, capped deterministically. Edge features include:

- radial-basis expansion of distance;
- sequence separation and same-chain flag;
- edge type: covalent sequence, spatial, disulfide, or cofactor relation;
- relative direction represented in local frames;
- relative frame orientation; and
- masks for uncertain coordinates.

A GVP-GNN is the recommended first equivariant encoder because it exposes scalar/vector reasoning clearly. GearNet or an invariant message-passing graph is the predeclared alternative. Hyperparameter searches compare no more than these families in milestone 1.

### 5.2 All-heavy-atom graph

Residue-level geometry loses pharmacophoric detail. The local atomic branch uses receptor heavy atoms only.

**Atom features**

- element and atom/residue identity;
- formal charge when reliably assigned;
- aromaticity, hybridization, ring membership;
- donor, acceptor, hydrophobe, cation, anion, and halogen-bond capabilities;
- metal/cofactor/water category;
- backbone/side-chain indicator;
- observed/repaired/predicted mask; and
- partial charge only in a separately versioned and licence-cleared feature set.

**Edges**

- covalent bonds within residues/cofactors;
- spatial radius neighbours with radial-basis distance;
- same-residue, same-chain, and atom-pair type; and
- equivariant direction vectors.

Use a shallow E(3)-equivariant or invariant radial message-passing network. Pool atom states to their parent residue through gated attention, while retaining a bounded number of atomic tokens for direct ligand cross-attention.

### 5.3 Global pocket descriptors

A transparent low-dimensional vector provides a useful residual pathway and baseline:

- pocket residue/atom count and volume proxy;
- surface area/depth and enclosure;
- hydrophobic/polar/charged fractions;
- donor/acceptor/aromatic counts;
- net formal charge proxy;
- metal/cofactor indicators;
- missing-coordinate fraction;
- detector confidence/rank; and
- structure quality/confidence summaries.

These features prevent the deep encoder from having to relearn basic statistics and make pocket-only leakage easy to audit.

### 5.4 Optional receptor-only pharmacophore field

After milestone 1, construct a receptor-centred sparse point/grid field with donor, acceptor, hydrophobe, aromatic, positive, negative, metal, and exclusion-volume channels. The coordinate frame and crop must be defined before the query ligand is loaded. A sparse 3D CNN or equivariant point model can pool it into pocket tokens.

This branch is inspired by receptor grids such as GRAIL/GRIPHIN but intentionally removes ligand coordinates and ligand-centred grids.

### 5.5 Deferred surface branch

A MaSIF/dMaSIF-style surface encoder is deferred because mesh/protonation/electrostatic pipelines introduce substantial operational and licensing risk. It enters only if:

1. residue+atom baselines are stable;
2. the surface toolchain is reproducible on a held-out structure audit;
3. licences permit the intended distribution; and
4. an ablation on dual-novel/apo tests justifies compute and failure cost.

## 6. Ligand representation

### 6.1 Strict ligand graph

From standardized SMILES build only:

- atom identity, formal charge, aromaticity, hybridization, valence, chirality, ring and hydrogen-count features;
- bond type, aromaticity, conjugation, stereochemistry, ring flag;
- graph topology and shortest-path/separation encodings;
- transparent global descriptors; and
- validity/domain masks.

No ETKDG, force field, distance geometry, conformer embedding, 3D descriptor, or coordinate-bearing SDF may be called.

### 6.2 gMolAI v2.0 integration

The reviewed gMolAI implementation provides internal per-atom states and a released calibrated 384-dimensional global representation assembled from graph- and node-pooled blocks. Define a stable adapter:

- AtomTokens: n_atoms × d_atom, projected from the internal node representation;
- GlobalLigand: 384, produced by the released calibration path;
- RawChem: transparent descriptors and counts; and
- LigandMask: unsupported chemistry and standardization flags.

Pin atom ordering so atom tokens align with explicit chemical features. Unit tests compare adapter output to the upstream checkpoint on reference SMILES.

Training schedule:

1. freeze gMolAI and train only projections/fusion/heads;
2. add small adapters or unfreeze the final graph blocks at a lower learning rate;
3. consider full discriminative fine-tuning only after locked validation shows benefit; and
4. always retain a frozen gMolAI-only baseline.

Because gMolAI’s pretraining corpus is vast, disclose overlap or nearest-neighbour exposure for benchmark ligands where knowable. This is a molecular prior, not evidence that pocket information helped.

### 6.3 Fallback ligand encoders

Maintain interchangeable adapters for:

- ECFP plus a shallow learner;
- a randomly initialized or task-trained GINE;
- a SMILES transformer; and
- a simple descriptor vector.

The scientific result should not depend on one proprietary or unclearly licensed checkpoint. Adapter equivalence also tests whether improvement comes from the fusion/pocket design or the ligand prior.

## 7. Pose-free fusion

### 7.1 Inputs

Let P_r be residue pocket tokens, P_a selected atomic pocket tokens, l_j ligand atom tokens, p_g global pocket vector, and l_g global ligand vector. Protein tokens contain receptor geometry; ligand tokens contain only topology.

### 7.2 Baseline fusion ladder

Run in order:

1. **additive baseline:** independent ligand and pocket predictions added;
2. **late concatenation:** concatenate pooled pocket and ligand vectors into an MLP;
3. **bilinear fusion:** low-rank outer product of pooled vectors;
4. **bidirectional cross-attention:** ligand atoms query pocket tokens and vice versa;
5. **pair-state model:** explicit residue/atom × ligand-atom compatibility tensor.

No later model proceeds if it cannot reliably improve the previous step on frozen validation components.

### 7.3 Cross-attention design

Project protein and ligand tokens into a shared latent space. Attention logits may include learned chemistry terms derived from token types, but not receptor–ligand distances:

\[
a_{ij} = \frac{(W_q p_i)^\top(W_k l_j)}{\sqrt d} + b_{\text{chem}}(p_i,l_j).
\]

The chemistry bias can use donor/acceptor, charge, hydrophobe, aromatic, metal, and size compatibility encoded from each side independently. Protein self-attention receives 3D geometric bias; ligand self-attention receives graph-distance/bond bias. Cross-attention receives neither relative coordinates nor a fabricated pose.

Several blocks alternate:

- geometry-aware protein self-update;
- topology-aware ligand self-update;
- protein-to-ligand and ligand-to-protein cross-update; and
- gated residual/global updates.

Pool with masked attention and combine the interaction summary with p_g and l_g. Regularize cross-attention entropy cautiously; sparse attention is not automatically chemically correct.

### 7.4 Pair-state extension

For bounded pocket and ligand sizes, initialize pair state z_ij from projected token products, differences, and independently derived chemistry compatibility. Update z through:

- rows sharing one pocket token;
- columns sharing one ligand atom;
- protein–protein geometric edges; and
- ligand–ligand bond/topological edges.

This allows coherent compatibility patterns without placing ligand atoms in space. Complexity is O(N_p × N_l × d), so cap tokens and profile memory. This is version 1.5, not the starting architecture.

### 7.5 Multi-scale integration

Two reasonable implementations should be predeclared:

- **hierarchical:** atom branch pools to residues; ligand cross-attends mainly to residues, with top atomic detail pooled into each residue;
- **dual-token:** ligand attends to residue and atom token sets with separate heads, followed by gated fusion.

Begin hierarchical for stability. Adopt dual-token only if atomic detail improves dual-novel and metal-free strata rather than merely fitting holo structures.

## 8. Prediction heads and losses

### 8.1 Endpoint-conditioned distribution

Use endpoint token e and optional assay-context vector a. For each endpoint, predict mean mu and positive scale sigma, for example with softplus plus a floor:

\[
(\mu,\sigma) = h(z_{\text{pair}}, p_g, l_g, e, a).
\]

Separate output calibrators are fitted for pKd, pKi, and pIC50. KIBA uses an independent native-scale head.

### 8.2 Censored regression

For exact y under a Gaussian observation model, use negative log likelihood. For interval [low, high], optimize the probability mass:

\[
-\log\left[\Phi\left(\frac{high-\mu}{\sigma}\right)-
\Phi\left(\frac{low-\mu}{\sigma}\right)\right].
\]

One-sided qualifiers use the corresponding CDF or survival function. Numerically stable log-CDF implementations and unit tests are mandatory.

### 8.3 Auxiliary objectives allowed in the strict core

Allowed objectives use no ligand coordinates:

- binary binder/activity classification;
- within-target pairwise ranking;
- masked residue/atom property prediction;
- receptor structure denoising or contrastive pocket learning;
- ligand masked-graph or gMolAI distillation;
- pocket–ligand contrastive learning from labelled or defensible negative pairs;
- endpoint/assay consistency; and
- replicate/noise modelling.

Not allowed in the strict core:

- native contact maps;
- protein–ligand distances;
- ligand internal 3D distances;
- pose RMSD or coordinate reconstruction; and
- ligand-centred pocket labels.

Loss weights should be fixed from validation-independent scale estimates or learned with explicit uncertainty weighting, then ablated.

### 8.4 Ranking and classification

Affinity regression alone does not optimize prospective triage. Add a within-target ranking loss only for measurements whose endpoint and assay context are comparable. For screening data, use class-balanced or prevalence-aware classification and report calibration. Do not create random negatives across targets and call them inactive.

## 9. Pretraining strategy

### 9.1 Pocket-only pretraining

Use PLINDER/BioLiP2 pockets without affinity labels for:

- masked residue/atom recovery;
- coordinate-noise denoising while preserving E(3) symmetry;
- contrastive matching of holo/apo/predicted views of the same site;
- pocket-cluster discrimination with debiasing;
- binding-site versus non-site discrimination; and
- residue chemical-environment prediction.

No ligand coordinates or ligand-defined contact map enters the strict pretraining track. If a site was originally located by a ligand, the ligand is removed before features are built and this provenance is retained.

### 9.2 Pair pretraining

BigBind or curated BindingDB/ChEMBL mappings can train binder classification and cross-modal alignment. Negatives should combine:

- experimentally confirmed inactives where available;
- same-assay thresholded records with qualifiers;
- cross-pocket negatives sampled to minimize trivial target/family shortcuts; and
- hard chemistry neighbours with uncertain status down-weighted or treated as positive-unlabelled data.

A positive-unlabelled or Bayesian treatment is preferable to asserting that every unmeasured pair is inactive.

### 9.3 Fine-tuning

Fine-tune on high-confidence S0/S1 quantitative records with endpoint-aware censored loss. Lower-confidence mappings can remain in pretraining or receive confidence-dependent weights. Hyperparameters are chosen on validation components, never CASP16/BayesBind/LIT-PCBA.

## 10. Uncertainty and applicability

Use three layers:

1. **Aleatoric:** endpoint-conditioned sigma learned from data and replicate noise.
2. **Epistemic:** a small deep ensemble across seeds/data resamples or a validated approximation.
3. **Input uncertainty:** variation across plausible receptor conformations, site radii, and detector top-k pockets.

Calibrate on held-out target/pocket components using grouped conformal prediction or endpoint-specific scaling. Calibration sets must exchange at the target-cluster level; ordinary row-wise conformal assumptions are inappropriate.

Applicability diagnostics include:

- maximum gMolAI/ECFP similarity to training ligands;
- protein sequence/domain identity;
- pocket structural and learned-embedding similarity;
- mapping confidence tier;
- structure quality/confidence and missingness;
- ensemble disagreement; and
- distance in the joint representation.

These are warnings and stratification variables, not guarantees of correctness. Report performance versus each diagnostic and define an abstention policy before external testing.

## 11. Holo, apo, and predicted structures

For matched systems, use the same ligand SMILES and transferred biological site on:

- a holo receptor with all ligands removed;
- a matched apo receptor aligned to the holo receptor; and
- an AlphaFold or other predicted structure with confidence features.

Measure absolute prediction change, rank stability, calibration, and failure rate. Train with conformation augmentation only from training targets. A model that succeeds only on holo geometry may be recognizing crystallographic induced fit rather than general receptor compatibility.

Coordinate perturbation tests should include small Gaussian noise, side-chain masking, missing loops, alternative protonation, and pocket-boundary changes within realistic ranges.

## 12. Model versions

### iScore3.0-v0 — diagnostic baselines

- gMolAI-only regression/classification;
- ECFP shallow models;
- pocket-only descriptors and residue graph;
- additive, concatenation, and low-rank bilinear fusion;
- BANANA-style pocket/ligand model; and
- sequence-only PSICHIC/GraphDTA-style control.

**Exit criterion:** data/splits stable; all leakage and ligand-coordinate tests pass; confidence intervals generated automatically.

### iScore3.0-v1 — reference strict model

- residue GVP pocket encoder;
- local heavy-atom pocket branch pooled hierarchically;
- frozen gMolAI atom/global ligand adapter;
- bidirectional pose-free cross-attention;
- endpoint-specific censored heteroscedastic heads;
- ensemble uncertainty and applicability diagnostics; and
- supplied-site primary evaluation with apo/predicted stress tests.

**Exit criterion:** statistically reliable improvement over v0 on predeclared dual-novel validation and at least one external set without worse calibration or prohibitive failure rate.

### iScore3.0-v1.5 — research extensions

- explicit pair-state tensor;
- receptor-only pharmacophore field;
- pocket-only pretraining at PLINDER scale;
- top-k pocket marginalization;
- adapter/full ligand fine-tuning; and
- optional surface branch.

Each extension is accepted independently through ablation. The final model is not necessarily the largest one.

## 13. Inference and caching

Precompute one PocketArtifact per receptor/site/conformation containing processed graphs, encoder states, global features, checksums, and warnings. Precompute one LigandArtifact per standardized SMILES containing graph, gMolAI states, global features, and schema version.

Scoring then runs only fusion and heads. For one pocket against a large library:

- cache pocket keys/values for cross-attention;
- bucket ligands by atom count;
- batch ligand encodings and fusion;
- retain deterministic input order and failure records;
- stream outputs rather than holding the full matrix; and
- report wall time, energy/memory, and compounds per second on named hardware.

Do not set performance targets before profiling v0/v1.

## 14. Proposed research API

Conceptual interfaces:

- prepare_receptor(structure, structure_policy) → ReceptorArtifact
- define_site(receptor_artifact, site_manifest) → PocketArtifact
- detect_sites(receptor_artifact, detector_config) → list of PocketArtifact
- prepare_ligand(smiles, ligand_policy) → LigandArtifact
- score(pocket_artifact, ligand_artifacts, endpoint, assay_context) → PredictionBatch
- explain(prediction_id) → compatibility hypotheses and applicability record

Each artifact is immutable and content-addressed. API validation rejects ligand coordinate fields rather than ignoring them.

## 15. Failure handling

Return an explicit failure or warning code for:

- invalid/ambiguous SMILES, mixtures, unsupported element or excessive size;
- missing or implausible site;
- too many missing pocket atoms/residues;
- covalent attachment;
- metal/cofactor or membrane context outside validated strata;
- receptor and site identifier mismatch;
- model/feature schema mismatch; and
- prediction outside calibrated applicability thresholds.

A numerical score without a validated input is more dangerous than an abstention.

## 16. Technical decision

The recommended first pocket encoding is a GVP-style residue graph plus a receptor heavy-atom graph, hierarchically fused and cross-attended by gMolAI atom tokens. This directly answers the user’s protein-site encoding question while keeping the implementation auditable.

A surface mesh, pharmacophore grid, pair tensor, and large pretrained pocket model are valuable follow-on experiments, not prerequisites. The project should first prove that receptor geometry contributes information on hard and conformation-stressed splits.
