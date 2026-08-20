# Scope, Definitions, and Requirements

## 1. Purpose

This document freezes the scientific information boundary for iScore3.0 before data curation or model selection. Without this boundary, a nominally “docking-free” system can still receive a ligand pose, ligand conformer, co-complex interaction graph, or ligand-conditioned pocket and thereby answer an easier and different question.

## 2. Primary prediction task

For protein structure P, binding-site definition S, ligand molecular graph L derived from SMILES, endpoint metadata E, and optional assay context A, learn

\[
f(P,S,L,E,A) \rightarrow (\hat y,\hat \sigma,\text{applicability})
\]

where y is a specified affinity endpoint on a negative base-10 molar scale, sigma is predictive uncertainty, and applicability summarizes whether the protein pocket and ligand lie within the training domain.

The core model may use the receptor’s 3D coordinates. It may use ligand atom/bond connectivity, stereochemical SMILES annotations, formal charge, calculated 2D descriptors, fingerprints, and learned graph/SMILES embeddings. It may not use ligand Cartesian coordinates or any feature derived from them.

## 3. Terminology

### 3.1 Docking-free

No docking engine, pose search, pose optimization, or docking score is invoked as an input to or component of the affinity predictor. Merely skipping an external docking executable is insufficient if a crystallographic ligand pose is supplied.

### 3.2 Pose-free

The predictor receives no protein–ligand relative coordinates and constructs no protein–ligand distance graph. The model can learn chemical compatibility but cannot inspect a proposed complex.

### 3.3 Strict ligand-3D-free

No ligand conformer, coordinates, 3D descriptors, distance matrix, pose-derived contact label, ligand-centred receptor crop, or pretrained ligand embedding whose production requires coordinates is used in training, validation, test, or inference for the core model.

This strict definition is stronger than “no ligand pose at inference.” It deliberately excludes training-only pose supervision from the principal claim.

### 3.4 Binding site

A set of receptor atoms or residues supplied independently of the query ligand. A site may originate from a curated site annotation, a receptor-only pocket detector, a historical reference ligand that is not the query, or an expert-defined region. Its provenance must be recorded.

### 3.5 Affinity

A measurement with an explicit endpoint, relation qualifier, numeric value, and unit. Kd, Ki, IC50, and KIBA are different outcomes. They must not be pooled as if they were exchangeable measurements of one latent quantity without an explicit model and sensitivity analysis.

## 4. Operating modes

| Mode | Inputs | Pocket source | Primary use | Report separately? |
|---|---|---|---|---|
| A: supplied site | receptor 3D, site definition, SMILES | curated or user supplied | primary development and scientific comparison | yes |
| B: discovered site | complete receptor 3D, SMILES | receptor-only top-k detector | virtual screening where site is unknown | yes |
| C: conformation stress test | holo, apo, or predicted receptor; site; SMILES | transferred or detected | robustness and translational validity | yes |
| D: sequence fallback | protein sequence, SMILES | none | baseline or missing-structure fallback | yes; not called the main model |

The main performance claim applies to Mode A. A system that selects the best result after seeing several pockets receives an additional multiple-instance advantage; top-1, top-3, and oracle-pocket results must therefore be distinguished.

## 5. Information-boundary matrix

| Information | Core model | Optional ablation | Prohibited for strict claim |
|---|---:|---:|---:|
| receptor atomic coordinates | allowed | — | — |
| receptor residue sequence | allowed | — | — |
| backbone frames, dihedrals, distances, orientations | allowed | — | — |
| receptor confidence/B-factor and alternate-location metadata | allowed | — | — |
| receptor-only surface or electrostatic features | allowed if reproducible | yes | — |
| ions, metals, waters, cofactors | policy-controlled | yes | — |
| query-ligand SMILES and 2D graph | allowed | — | — |
| gMolAI atom/global embeddings generated from the 2D graph | allowed after licence verification | yes | — |
| 1D/2D descriptors and fingerprints | allowed | baseline | — |
| ligand conformer or 3D descriptor | no | non-core comparator only | prohibited |
| crystallographic or docked query-ligand coordinates | no | non-core comparator only | prohibited |
| protein–query-ligand contact/distance graph | no | non-core comparator only | prohibited |
| training target computed from ligand coordinates | no | privileged-information study only | prohibited |
| ligand-centred grid/crop based on query coordinates | no | non-core comparator only | prohibited |
| label or analogue from the locked test target | no | — | prohibited |

## 6. Pocket-definition policy

A supplied pocket should be represented by one of the following stable forms:

1. receptor chain identifiers plus residue identifiers;
2. a receptor-coordinate centre and fixed radius;
3. a receptor-only pocket detector identifier and rank; or
4. a reference-site mapping with documented source and sequence/structure alignment.

A crystallographic reference ligand may be used once to annotate the biological site when that ligand is not the query and when the evaluation scenario genuinely assumes a known site. The query ligand must never determine its own crop. For unbiased apo/predicted evaluation, the holo-defined site must be transferred by receptor alignment without transferring ligand coordinates or ligand-conditioned features.

Default pocket shells should be evaluated at several radii rather than tuned per test target. Residues or atoms at unresolved coordinates are flagged, not silently imputed from the bound ligand.

## 7. Protein structure policy

The structure pipeline must:

- preserve biological chain, residue, insertion-code, alternate-location, and model provenance;
- standardize nonstandard residues conservatively;
- define explicit policies for protonation, missing atoms/loops, waters, metals, covalent ligands, and cofactors;
- retain a mask distinguishing observed, repaired, predicted, and low-confidence coordinates;
- prevent any feature generator from reading the query-ligand file; and
- hash both the input structure and processed pocket artifact.

Holo receptors are permitted for development but are a potential privileged signal because binding can induce a ligand-specific conformation. Matched apo and predicted-structure tests are therefore mandatory.

## 8. Ligand policy

SMILES standardization must be versioned and deterministic. The pipeline should retain the original record and produce:

- parsed molecular graph and sanitization status;
- largest-fragment or mixture decision;
- normalized charge and tautomer policy;
- stereochemistry and isotopic information;
- canonical, non-isomeric, and standardized SMILES where applicable;
- InChIKey and connectivity block;
- Bemis–Murcko scaffold;
- fingerprint and similarity metadata for splitting; and
- atom mapping used to align raw features with learned gMolAI tokens.

No conformer-generation call may occur in the strict feature pipeline, including inside an ostensibly pretrained model. Tests should fail if a coordinate-bearing ligand object reaches the model interface.

## 9. Endpoints and labels

Preferred primary endpoints are pKd and pKi:

\[
pK_x=-\log_{10}(K_x[\mathrm{mol/L}]).
\]

pIC50 is a secondary endpoint with an endpoint token or separate head because IC50 depends on assay conditions and mechanism. KIBA is an aggregated score and remains in its native scale. EC50, percent inhibition, docking scores, and qualitative activity labels are excluded from regression unless assigned to explicitly separate tasks.

Required label fields include endpoint, relation operator, numeric value, unit, converted value, assay identifier/type, temperature where known, pH where known, source, publication/patent, measurement date or earliest public date, and curation status.

Censored values such as “>10 µM” are intervals, not point observations. Replicates should remain linked; robust consensus labels and their dispersion may be derived without erasing the underlying records.

## 10. Outputs

For each receptor-site/SMILES pair, the research API should return:

- predicted endpoint and scale;
- predictive mean and calibrated interval;
- assay/endpoint head used;
- pocket source and rank;
- ligand, pocket, and pair applicability indicators;
- nearest-training ligand and pocket similarities;
- model, feature schema, and data-manifest versions; and
- structured warnings for unsupported chemistry, poor structures, metals, covalency, or out-of-domain inputs.

Optional atom–residue compatibility maps may be exposed as hypotheses. They must not be called interaction poses or mechanistic explanations without independent validation.

## 11. Intended use

Initial intended use is prospective prioritization within medicinal-chemistry campaigns when a plausible receptor structure and binding site are available but high-throughput docking is undesirable. The model supports ranking and triage; it does not replace experimental affinity measurements.

Initial chemical domain is drug-like, noncovalent, single-component organic molecules representable by sanitizable SMILES. The receptor domain is structured protein pockets with ordinary amino acids and documented heteroatom handling.

## 12. Exclusions for version 1

Unless a dedicated work package is approved, version 1 excludes:

- covalent inhibitors and reaction-aware scoring;
- peptides, proteins, nucleic acids, polymers, mixtures, and ill-defined salts;
- membrane/lipid free-energy effects not represented by the input;
- induced-fit pose or binding-mode prediction;
- kinetics, residence time, permeability, selectivity, toxicity, or efficacy claims;
- mutation effects without adequate matched data;
- metalloprotein chemistry requiring quantum treatment; and
- absolute free-energy interpretation beyond the calibrated empirical endpoint.

These cases should be detected and reported rather than silently scored.

## 13. Functional requirements

The implementation shall:

- accept PDB/mmCIF receptor structures, an explicit site definition, and SMILES;
- run without ligand coordinates and without a docking executable;
- batch many ligands against a cached pocket representation;
- expose deterministic preprocessing and schema validation;
- support Kd, Ki, and IC50 heads with censor-aware training;
- export uncertainty and applicability diagnostics;
- reproduce predictions from immutable configuration and artifact hashes; and
- make the strict information boundary machine-testable.

## 14. Non-functional requirements

The project should target:

- receptor preprocessing once per site;
- ligand preprocessing independent of receptor;
- cached scoring suitable for library-scale triage;
- mixed-precision GPU inference with a documented CPU fallback;
- no test-set network calls or mutable external lookups;
- seed, hardware, package, checkpoint, and data provenance capture;
- unit, integration, scientific-regression, and leakage tests; and
- licence-compatible, independently replaceable components.

Throughput targets should be set only after the reference baselines are profiled on named hardware.

## 15. Minimum claim language

Permitted after successful validation:

> iScore3.0 is a pose-free affinity model that combines a three-dimensional receptor-pocket representation with a coordinate-free ligand representation derived from SMILES.

Not permitted without further evidence:

- “the first docking-free scoring function”;
- “predicts binding free energy”;
- “generalizes to unseen proteins and chemotypes” without dual-novel tests;
- “learns interactions” solely from an attention map; or
- “outperforms docking” using biased decoys or overlapping structures.

## 16. Scope-freeze decisions required from the committee

Before implementation begins, record approval or modification of:

1. strict prohibition of ligand 3D during both training and inference;
2. Mode A as the primary task and Modes B/C as secondary evaluations;
3. pKd/pKi as primary labels and pIC50 as a separate endpoint;
4. version-1 exclusions, especially covalent ligands, metals, waters, and cofactors;
5. commercial versus open-only data strategy;
6. intended licensing and release model for code, weights, and processed metadata; and
7. whether gMolAI v2.0 may be reused and relicensed in the intended distribution.

The decision record belongs in 07_risks_decisions_and_governance.md. Any later change to the information boundary creates a new experimental track and must not overwrite strict-core results.
