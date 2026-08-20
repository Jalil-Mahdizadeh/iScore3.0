# Risks, Decisions, and Governance

## 1. Governance objective

iScore3.0 combines restricted structural resources, noisy biochemical measurements, large pretrained encoders, and fast-moving prior art. Governance is therefore part of scientific validity. This document identifies conditions under which the project should proceed, revise scope, or stop.

Statuses are proposals until the scientific committee records a decision.

## 2. Committee decision register

| ID | Decision requested | Recommended decision | Status | Committee record |
|---|---|---|---|---|
| D-001 | Correct project directory name | use blueprint/ | implemented | user confirmed the corrected spelling |
| D-002 | Primary task | receptor 3D + independently supplied site + ligand SMILES | proposed | pending |
| D-003 | Ligand information boundary | prohibit ligand coordinates and pose-derived labels in strict training and inference | proposed | pending |
| D-004 | Secondary modes | evaluate receptor-only site detection and apo/predicted structures separately | proposed | pending |
| D-005 | Primary endpoints | pKd and pKi; pIC50 separate; KIBA native scale | proposed | pending |
| D-006 | Version-1 chemistry | noncovalent, drug-like, single-component organic ligands | proposed | pending |
| D-007 | Protein heteroatoms | explicit metal/cofactor policy; waters as later ablation | proposed | pending |
| D-008 | Structural data | licensed PDBbind core plus open complements if agreement accepted | proposed | pending |
| D-009 | Evaluation | dual-novel union components plus temporal/external tests; random split diagnostic only | proposed | pending |
| D-010 | Pocket encoder | residue GVP + local heavy-atom graph; surface branch deferred | proposed | pending |
| D-011 | Ligand encoder | pinned gMolAI atom/global adapter after licence resolution | proposed | pending |
| D-012 | Fusion | baseline ladder ending in pose-free cross-attention for v1 | proposed | pending |
| D-013 | Uncertainty | endpoint-specific interval model, ensemble, applicability, abstention | proposed | pending |
| D-014 | Reuse/release | no third-party code, weights, or derived data released without rights review | proposed | pending |
| D-015 | External labels | access-controlled/frozen one-pass evaluation | proposed | pending |
| D-016 | Project success | structural benefit must exceed predeclared practical effect on hard/external tests | proposed | pending |
| D-017 | Gate-0/1 progression | do not start full architecture; authorize only bounded remediation | recommended | evidence complete; committee decision pending |

For each pending item, the committee should record approve, modify, defer, or reject; date; responsible person; justification; and consequences.

## 3. Risk register

Probability and impact are qualitative before WP1 pilots.

| ID | Risk | Probability | Impact | Early indicator | Primary mitigation |
|---|---|---:|---:|---|---|
| R-01 | broad novelty already published | high | high | direct method with same modalities and rigorous results | narrow claim; continual surveillance; emphasize strict boundary/evaluation |
| R-02 | ligand 3D leaks through preprocessing/pretraining | medium | critical | coordinate fields, conformer calls, pose-derived targets | machine-enforced strict interfaces and lineage tests |
| R-03 | train/test chemical or pocket leakage | high | critical | close neighbours dominate test predictions | union-component splits, independent similarity audit, per-test neighbour report |
| R-04 | holo receptor encodes query-induced conformation | high | high | performance collapses on apo/predicted structures | matched conformation tests and training augmentation |
| R-05 | wrong activity-to-structure/site mapping | high at scale | critical | label assigned to incompatible construct/site | S0–S4 tiers, sequence/domain/site audit, quarantine ambiguity |
| R-06 | endpoints are mixed or units/qualifiers mishandled | high | critical | implausible pK, contradictory replicates | strict schema, interval conversion tests, separate heads |
| R-07 | putative negatives/decoys create ligand bias | high | high | ligand-only model matches full model | confirmed inactives/PU learning, BayesBind, audited/rebuilt LIT-PCBA, negative controls |
| R-08 | gMolAI terms do not permit intended reuse | medium | high | no explicit licence for code/checkpoint/calibrator | add/obtain licence or use replaceable ligand adapter |
| R-09 | competitor code has restrictive/unclear licence | medium | high | custom, missing, or non-commercial licence | clean-room concepts; legal inventory; no copying |
| R-10 | PDBbind/derivative redistribution is restricted | high | high | agreement limits structures/features | controlled storage; distribute IDs/scripts only; rights review |
| R-11 | PLINDER affinity labels are incorrect | known | critical if used | open issue #94/disabled query field | prohibit field until fixed or independently remapped |
| R-12 | surface/electrostatics pipeline is fragile | high | medium | preprocessing failures/version drift | defer branch; containerize; residue+atom baseline |
| R-13 | protein 3D adds no value beyond ligand/target priors | medium-high | high | H1/H2 intervals include no useful gain | early v0 gate; publish negative result or revise task |
| R-14 | target-family imbalance drives global metrics | high | high | protein-only/target-mean performance is strong | target-macro metrics, cluster resampling, balanced sampling |
| R-15 | uncertainty is overconfident OOD | high | high | poor interval coverage on novel targets | grouped calibration, ensembles, risk–coverage, abstention |
| R-16 | site detector error obscures affinity model | medium | medium-high | correct-site performance high, Mode B low | separate detector recall and conditional scoring |
| R-17 | multiple binding sites/mechanisms are unresolved | medium | high | same target maps to inconsistent sites | site-specific records, multi-instance sensitivity, quarantine |
| R-18 | metals/cofactors/waters/covalency exceed model scope | medium | high | clustered large errors/failures | scope detection, strata, explicit exclusions/branches |
| R-19 | compute/storage expands before feasibility | medium | medium | large PLINDER/PDBbind feature build before v0 | pilot 1–5%, gated scale, cache/deduplicate |
| R-20 | mutable databases make results irreproducible | high | high | counts/records change between runs | release snapshots, checksums, immutable manifests |
| R-21 | external benchmarks become indirectly tuned | medium | critical | repeated runs or result-driven preprocessing | evaluation owner, freeze, immutable first predictions |
| R-22 | attention maps are overinterpreted | high | medium | claimed contacts without pose evidence | label as compatibility hypotheses; independent validation |
| R-23 | experimental label noise caps ranking | high | medium-high | replicate/series disagreement | interval/noise model, report estimated ceiling, avoid overclaim |
| R-24 | receptor preprocessing encodes dataset/source | medium | high | model predicts source/resolution/label from metadata | source baselines, harmonization, source-held-out tests |
| R-25 | pretrained corpus contains benchmark identities | high | medium-high | exact test ligands/pockets in pretraining | exposure disclosure, exclusion/similarity sensitivity |
| R-26 | irreproducible dependency or checkpoint | medium | high | upstream moves/disappears | pin revision, hash assets, container/SBOM, permitted archive |
| R-27 | performance metric does not match use | medium | high | global RMSE improves but target ranking does not | intended-use endpoints, macro ranking, screening calibration |
| R-28 | committee scope changes after results | medium | high | information boundary or primary metric moves | ADR/amendment process; preserve original confirmatory track |

## 4. Detailed critical risks

### 4.1 Novelty and prior art

The strongest direct precedents include PSG-BAR, PLMCA, AttentionMGT-DTA, BANANA, CASTER-DTA, 3DProtDTA, HoloProt, Graph_RG, and PLANET. The project must not claim to invent the modality combination. The literature review should be refreshed:

- before WP3 architecture freeze;
- before external evaluation;
- before abstract/manuscript submission; and
- during peer review if the field changes.

A surveillance record includes query, source, date, candidate, classification, code, and effect on claims. PLANET v2.0 and post-2026 systems deserve particular attention.

**Trigger for scope review:** a public method matches the strict information boundary, multi-scale pocket design, gMolAI-like atom/global fusion, and hard evaluation before iScore3.0 is submitted.

**Response:** emphasize independent reproduction/data/uncertainty findings, compare directly, or change the research question. Do not manufacture novelty through terminology.

### 4.2 Privileged ligand information

Leakage can occur through:

- a ligand SDF that silently carries coordinates;
- RDKit conformer generation inside a pretrained encoder;
- ligand 3D descriptors;
- contact/distance auxiliary labels;
- a ligand-centred pocket crop/grid;
- holo side chains conditioned on the exact query ligand;
- docking scores used as features; or
- complex-pretrained embeddings.

The first four violate the strict core outright. Holo conformation is permitted receptor information but must be disclosed and stress-tested. Feature-lineage tests are release blockers.

### 4.3 Activity-to-pocket mapping

Target databases often identify a protein, not a physical site or exact construct. PDB occurrence does not prove that an assay was performed on that structure. Incorrect mapping produces highly confident but scientifically meaningless pairs.

All scaled records need mapping tiers and reason codes. S3/S4 cannot enter locked structural evaluation. The proportion of training data by mapping tier is reported with an ablation. A large performance gain from low-confidence mapping is treated skeptically.

### 4.4 Data licences

At least four layers can differ:

1. database software licence;
2. curated metadata licence;
3. imported upstream data licence;
4. structure/checkpoint/derived-artifact terms.

For example, PLINDER code, PLINDER-curated data, BindingDB imports, ChEMBL imports, PDB coordinates, and PDBbind curation do not share one licence. “Open on GitHub” is not a licence.

The data steward maintains a rights matrix for internal use, publication of aggregate results, release of IDs/metadata, release of coordinates/features, model weights, and commercial use. Unknown means not approved.

### 4.5 Scientific signal failure

It is plausible that pose-free receptor geometry contributes little once ligand and target priors are controlled. The project is designed to detect this early.

Possible outcomes:

- geometry helps dual-novel/apo targets: proceed;
- geometry helps only random/holo data: narrow to ligand-conditioned structural use or stop headline claim;
- no reliable geometry benefit: publish benchmark/negative result, improve data mapping, or redirect to screening/site-conditioned classification;
- benefit exists only in certain families/cofactor classes: define a restricted applicability domain.

Do not add pose-derived labels merely to rescue the strict claim; that creates a different model track.

## 5. Stop, revise, and proceed rules

### Stop or pause

- required data or gMolAI rights cannot support intended use/release;
- ligand-coordinate leakage cannot be eliminated from the chosen encoder;
- target-to-site mappings fail audit at an unusable rate;
- dual-novel splits leave insufficient independent targets for inference;
- external-test governance has been compromised;
- reproducibility cannot be achieved from pinned artifacts; or
- full model fails to improve on strong ligand-only/pocket-prior baselines beyond the approved practical threshold after predeclared v1 work.

### Revise scope

- useful results occur only for supplied holo sites;
- only ranking, not absolute affinity, is reliable;
- performance is family-specific;
- pIC50 assay variability overwhelms shared learning;
- top-k pocket detection dominates total error;
- metals/covalency form a distinct learnable subproblem; or
- a new direct competitor changes the novelty position.

### Proceed

- data rights and lineage are clear;
- strict feature and leakage tests pass;
- geometry adds replicated value on hard validation;
- conformation sensitivity is understood and reflected in uncertainty;
- external results beat non-structural/nearest-neighbour baselines; and
- outputs are calibrated enough for the declared prioritization use.

## 6. Data and model release policy

Three release levels should be considered independently:

### Level A — fully open implementation

Release source, synthetic fixtures, configuration schemas, split-generation logic, evaluation code, documentation, and permitted small metadata. This should be the default target.

### Level B — restricted-data reproducibility

Release acquisition instructions, source IDs, checksums, exclusion lists, and scripts that registered users can run against their own licensed data. Do not release restricted coordinates or reversible derived features.

### Level C — weights/predictions

Release only if training-data and upstream-checkpoint terms permit it and membership/privacy analysis finds no unacceptable risk. Publish external predictions and aggregate metrics when benchmark terms allow, even if weights cannot be distributed.

The repository licence covers original code only. THIRD_PARTY_NOTICES and the data rights matrix remain separate.

## 7. Scientific integrity

The project will:

- report all preregistered models and seeds;
- preserve the first external prediction;
- disclose failed inputs and exclusions;
- keep random-split and biased-decoy findings clearly labelled;
- distinguish statistical from practical significance;
- avoid interpreting attention as a pose;
- avoid equating predicted pK with a physical free-energy calculation;
- cite datasets, software, and pretrained models;
- archive negative results and abandoned branches; and
- state all known test/pretraining overlap.

Corrections create new artifact versions and a visible erratum record.

## 8. Security and sensitive data

The planned public biochemical datasets are not personal data, but licences, unpublished partner assays, credentials, and proprietary structures may be sensitive.

- credentials live in approved secret storage, never code/config/logs;
- unpublished projects use separate access-controlled data roots and experiment tracking;
- external services are disabled unless approved;
- logs avoid raw restricted structures, full compound libraries, and credentials;
- releases undergo secret/data/licence scanning; and
- collaboration exports are reviewed by the data steward.

Any future clinical or patient-derived data requires a new governance plan.

## 9. Reproducibility review checklist

Before each gate:

- [ ] source versions and checksums resolve;
- [ ] licences and citations are complete;
- [ ] environment/container digest is available;
- [ ] workflow rebuilds a clean small fixture;
- [ ] strict ligand-coordinate tests pass;
- [ ] split/leakage audit passes independently;
- [ ] model/config/checkpoint hashes are recorded;
- [ ] seed-level results and failures are complete;
- [ ] statistical analysis uses correct independent units;
- [ ] external labels were not used for selection;
- [ ] model/data cards reflect current scope;
- [ ] decisions, deviations, and incidents are recorded.

## 10. Committee approval checklist

The committee is asked to approve or amend:

- [ ] the primary task and three operating modes;
- [ ] the strict no-ligand-3D rule at training and inference;
- [ ] endpoint hierarchy and censoring policy;
- [ ] chemical/protein exclusions;
- [ ] PDBbind/open-data acquisition strategy;
- [ ] gMolAI licensing and integration path;
- [ ] initial residue+atom pocket encoder;
- [ ] baseline and ablation suite;
- [ ] dual-novel/temporal/external evaluation protocol;
- [ ] minimal important effects and statistical plan;
- [ ] compute/storage pilot and work-package gates;
- [ ] code/data/weight release intent;
- [ ] named owners and independent evaluation control; and
- [ ] novelty surveillance schedule.

Approval should be recorded in a dated committee minute and mirrored as ADRs after repository initialization.

## 11. Open questions for the committee

1. Is the intended downstream use academic/open research, internal institutional work, or eventual commercial deployment?
2. Is a commercial PDBbind release authorized, or must the first project remain entirely open-data?
3. Should strict exclusion cover pose-derived training supervision, as recommended, or should a separately labelled privileged-information model also be funded?
4. Which assay endpoints and target classes are highest priority for the intended use?
5. Are metals, cofactors, conserved waters, covalent ligands, membrane proteins, and multi-site targets in version 1?
6. Can gMolAI code/checkpoints/calibrators receive an explicit licence compatible with the intended release?
7. Who will control external labels and approve confirmatory runs?
8. What practical effect size and library-screening use case justify continuation?
9. What compute, storage, staffing, and calendar constraints apply?
10. Which artifacts may be shared with the committee, collaborators, and public?

These questions do not prevent blueprint review. Their answers are Gate 1 inputs and should not be guessed during implementation.

## 12. Governance conclusion

The highest project risks are not neural-network selection. They are overbroad novelty, hidden ligand/complex information, wrong assay-to-site mappings, train–test redundancy, endpoint semantics, and data rights. The proposed gates make those risks observable early and keep the project valuable even if the final scientific result is negative.
