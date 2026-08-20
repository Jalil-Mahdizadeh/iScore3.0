# iScore3.0 — Executive Blueprint

**Committee review draft:** 20 August 2026  
**Proposed system:** pose-free prediction of protein–small-molecule affinity from a receptor binding site in 3D and a ligand represented only by its SMILES-derived 1D/2D information.

**Gate-0/1 addendum:** the bounded feasibility study is complete and did not meet the predeclared geometry-value criterion. Full architecture training is not recommended yet; the evidence supports only a narrow remediation phase. See [`reports/gate01/GATE_0_1_REPORT.md`](../reports/gate01/GATE_0_1_REPORT.md).

## Decision requested

Do not approve full architecture training from the current evidence. Consider approval of a bounded remediation phase to expand the number of independent high-confidence components, add validated structural-similarity and apo/predicted views, resolve gMolAI provenance/rights, and reproduce direct and modern sequence baselines. Progress remains conditional on the six gates below.

## Scientific conclusion from the review

The broad idea is feasible, useful, and no longer unprecedented. PSG-BAR, AttentionMGT-DTA, PLMCA, BANANA/BigBind, CASTER-DTA, 3DProtDTA, parts of HoloProt, and Graph_RG already combine receptor structure with a ligand graph without requiring a docked query complex. PLANET and BlendNet are highly relevant but use pose-derived or complex-interaction supervision during training. GRIPHIN, T-ALPHA, TankBind, DrugCLIP, PocketDTA, MMPD-DTA, AlignNet, LigUnity, and related methods inform the design but do not meet the strict “no ligand 3D anywhere in the core model” rule.

Accordingly, iScore3.0 should not be presented as the first docking-free or first protein-3D/ligand-2D scoring function. A defensible contribution is the combination of:

1. a strict ligand-3D-free information boundary at training, validation, and inference;
2. a pocket-focused, multi-scale E(3)-aware receptor encoder;
3. atom-level and global ligand representations from gMolAI v2.0;
4. pose-free cross-modal compatibility learning without protein–ligand distance edges;
5. explicit handling of pocket, receptor-conformation, endpoint, censoring, and prediction uncertainty; and
6. leakage-resistant, structure-aware and chemistry-aware evaluation including apo and predicted receptors.

Whether that package is novel enough for publication must be reassessed immediately before manuscript submission because this area is moving quickly.

## Proposed product definition

The primary task is:

> Given a protein receptor structure, an externally supplied binding-site definition, and a canonicalized ligand SMILES, estimate binding affinity without docking, generating a ligand conformer, or using any ligand coordinates.

The committee should distinguish three operating modes:

- **Mode A — supplied site:** a residue list or receptor-coordinate centre is provided. This is the primary scientific task.
- **Mode B — receptor-only site discovery:** P2Rank or another receptor-only method proposes top-k pockets before scoring.
- **Mode C — known target with uncertain conformation:** holo, matched apo, and AlphaFold/other predicted structures are scored to quantify conformation sensitivity.

Mode A avoids conflating affinity prediction with pocket detection. Modes B and C are essential translational stress tests but should be reported separately.

## Recommended model

The minimum credible architecture has four parts:

- **Pocket geometry:** a residue-level GVP/GearNet-style graph using backbone frames, distances, orientations, residue chemistry, solvent exposure, and optional protein-language-model features.
- **Pocket chemistry:** a local all-heavy-atom graph encoding element, formal charge, donor/acceptor/aromatic/hydrophobic character, metals, cofactors, and 3D neighbourhoods.
- **Ligand:** gMolAI atom embeddings plus its calibrated global embedding, augmented by transparent 2D chemical features. The encoder is frozen first and selectively fine-tuned only after stable baselines.
- **Pose-free fusion:** bidirectional cross-attention or a learned residue/atom compatibility tensor constrained by protein geometry and ligand topology, but never by receptor–ligand distances.

Separate calibrated heads should model pKd, pKi, and pIC50. Censored observations require interval-aware likelihoods. KIBA scores and other non-thermodynamic endpoints must not be silently converted into molar affinity.

## Data recommendation

No single dataset is sufficient.

- **Structural supervised core:** a licensed, version-pinned PDBbind release, complemented by LP-PDBbind and the official PDBbind CleanSplit/GEMS metadata.
- **Scale:** BindingDB, ChEMBL, Papyrus, or BigBind after rigorous target-to-structure and target-to-pocket mapping.
- **Structure and pocket pretraining:** PLINDER and BioLiP2/BioLiP2-Opt. PLINDER’s affected binding-affinity field must not be used until its documented parsing issue is fixed or independently remapped.
- **External evaluation:** CASP16 pharmaceutical affinity data, BayesBind, and carefully separated target series. Use LIT-PCBA only as a red-team diagnostic or after an independently rebuilt split passes audit; a 2025 reproducible preprint reports exact and analogue leakage. CASF-2016 remains legacy comparability, not proof of generalization.
- **Diagnostics only:** DUD-E, because its decoy construction can be exploited by ligand-only models.

Every observation must retain source, publication, assay type, units, qualifiers, construct/target mapping, structure identifier, site identifier, and licence lineage.

## Evaluation that can support a scientific claim

Random pair splits are unacceptable as headline evidence. The locked test sets should be connected components of a union graph that links examples by protein/pocket similarity or ligand scaffold/similarity, with a temporal holdout where dates permit. Report the nearest training similarity for every test point.

Required controls are:

- ligand-only gMolAI and fingerprint models;
- protein/pocket-only models;
- additive and concatenation baselines;
- BANANA reproduction or architecture-matched reimplementation;
- sequence-only PSICHIC/GraphDTA-style controls;
- holo versus matched apo versus predicted receptor structures;
- supplied pocket versus receptor-only top-1 and top-k pockets; and
- label-permutation and pocket-permutation negative controls.

Metrics must include within-target rank correlation, MAE/RMSE with target-cluster bootstrap confidence intervals, prospective classification/enrichment where appropriate, calibration, interval coverage, and failure rates. Pair-level bootstrap intervals overstate certainty and are not sufficient.

## Six approval gates

### Gate 1 — scope and rights

Approve only after the exact ligand-3D-free rule, intended-use statement, dataset licences, gMolAI reuse terms, and redistribution policy are signed off.

### Gate 2 — data integrity

Proceed to model development only when duplicate measurements, endpoint semantics, censoring, target/construct mapping, site assignment, and cross-source overlap have auditable tests.

### Gate 3 — baseline value

The 3D-pocket model must outperform ligand-only, protein-only, and simple late-fusion baselines on locked dual-novel splits with uncertainty intervals.

### Gate 4 — geometry value

Pocket geometry must add value on apo or predicted structures, not only on ligand-conditioned holo pockets. Performance must be robust to realistic coordinate perturbations and alternative pocket definitions.

### Gate 5 — external validity

At least one untouched external or time-forward benchmark must show useful ranking or calibrated affinity performance. CASF alone cannot satisfy this gate.

### Gate 6 — reproducibility

A fresh environment must regenerate manifests, splits, features, training, and evaluation from versioned configurations. Test-set labels remain inaccessible to model selection.

## Expected research outputs

The staged programme should yield:

1. a rights-cleared, provenance-rich benchmark and immutable split manifests;
2. a suite of leakage-diagnostic baselines;
3. a strict pose-free iScore3.0 reference implementation;
4. calibrated predictions with applicability-domain indicators;
5. a model card, data card, and negative-results record; and
6. a manuscript whose novelty claim is supported by contemporaneous competitor surveillance.

## Bottom line

The research question remains worthwhile, but the current pilot does not justify the full architecture. Its strongest opportunity is not the generic absence of docking; it is demonstrating, under hard leakage controls, when and why receptor-pocket geometry improves affinity prediction over molecular and sequence priors while the ligand remains entirely coordinate-free. That demonstration has not yet been achieved.
