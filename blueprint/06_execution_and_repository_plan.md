# Execution and Repository Plan

## 1. Delivery strategy

Implementation should proceed through gated work packages. Data engineering and baselines come before architectural expansion. The blueprint directory remains a read-only decision package after approval; implementation decisions are recorded rather than silently rewriting the approved scope.

Durations below are planning ranges for a small research team and should be revised after compute, staffing, and data-access decisions.

## 2. Work packages

### WP0 — governance, scope, and rights

**Duration:** 2–3 weeks  
**Depends on:** committee review

Tasks:

- approve the information boundary and intended use;
- decide open-only versus licensed PDBbind track;
- resolve gMolAI code/checkpoint/calibrator licence;
- audit all planned dataset and dependency terms;
- name data steward, evaluation owner, and model owner;
- approve compute/storage envelope;
- freeze endpoint and version-1 exclusion policies; and
- initialize Git, issue templates, code review, and decision log.

Deliverables:

- signed scope and rights matrix;
- approved data-source registry;
- architecture decision records 0001–0004;
- threat model for leakage and restricted data; and
- go/no-go Gate 1 record.

### WP1 — reproducible data foundation

**Duration:** 5–8 weeks  
**Depends on:** WP0 data access

Tasks:

- implement molecule, target, assay, structure, site, and observation schemas;
- build version-pinned ingestion and checksum manifests;
- implement strict ligand standardization with coordinate-prohibition tests;
- implement receptor/site normalization and pocket extraction;
- curate endpoint intervals and replicates;
- assign target-to-site confidence tiers;
- deduplicate sources and construct similarity graphs;
- produce random, ligand, pocket, dual-novel, and temporal splits; and
- create data cards and QC dashboards.

Deliverables:

- D0 open feasibility manifest;
- D1 licensed structural manifest if approved;
- immutable split files;
- synthetic and redistributable test fixtures;
- exclusion/provenance/licence reports; and
- Gate 2 data-integrity review.

### WP2 — baseline suite

**Duration:** 4–6 weeks  
**Depends on:** stable WP1 development/validation manifests

Tasks:

- implement ECFP and gMolAI-only models;
- implement target/pocket-only and nearest-neighbour models;
- implement transparent concatenation/bilinear models;
- reproduce BANANA or a protocol-equivalent baseline;
- implement sequence-only and 3DProtDTA-style controls;
- build common training/evaluation interfaces;
- profile caching, throughput, and failure rates; and
- publish seed-level baseline report.

Deliverables:

- iScore3.0-v0 model registry;
- reproducible baseline leaderboard;
- leakage and negative-control report;
- initial compute/storage forecast; and
- Gate 3 decision on whether pocket information is promising.

### WP3 — strict reference model

**Duration:** 6–10 weeks  
**Depends on:** WP2 and Gate 3

Tasks:

- implement residue GVP and invariant graph alternatives;
- implement local receptor heavy-atom branch;
- build the pinned gMolAI atom/global adapter;
- implement bidirectional pose-free cross-attention;
- implement endpoint-conditioned censored likelihood;
- implement ensemble uncertainty and applicability features;
- train only on approved train/validation components; and
- run predeclared ablations.

Deliverables:

- iScore3.0-v1 checkpoints and model card;
- geometry, atomic-resolution, and fusion ablation report;
- invariant/information-boundary test evidence;
- learning curves and calibration report; and
- Gate 4 geometry/conformation decision.

### WP4 — conformation, pocket, and robustness

**Duration:** 4–7 weeks  
**Depends on:** frozen v1 candidates

Tasks:

- assemble matched holo/apo/predicted structures;
- run transferred-site and receptor-only P2Rank/fpocket evaluation;
- quantify coordinate, missingness, protonation, cofactor, and shell sensitivity;
- calibrate top-k pocket/conformation uncertainty;
- define abstention thresholds; and
- freeze the confirmatory candidate and calibrator.

Deliverables:

- conformation/pocket robustness report;
- detector recall and conditional-affinity report;
- frozen external-evaluation bundle; and
- Gate 5 readiness approval.

### WP5 — scale and optional research extensions

**Duration:** 6–12 weeks, parallel after v1 stability  
**Depends on:** audited BindingDB/ChEMBL/BigBind or PLINDER pipelines

Candidate tasks:

- large-scale mapped activity or positive-unlabelled pretraining;
- PLINDER/BioLiP2 pocket-only pretraining;
- pair-state compatibility tensor;
- receptor-only pharmacophore field;
- gMolAI adapters/fine-tuning;
- surface branch after licence/toolchain approval; and
- efficiency/distillation studies.

Each extension has its own hypothesis, branch, configuration family, and ablation. None can change the v1 confirmatory result.

### WP6 — locked external evaluation and dissemination

**Duration:** 4–6 weeks  
**Depends on:** Gate 5 frozen bundle

Tasks:

- execute CASP16, BayesBind, approved target-series protocols, and either diagnostic-only or independently rebuilt LIT-PCBA analyses;
- generate target-cluster statistical analyses;
- audit all failed inputs and post hoc diagnostics;
- prepare reproducibility archive, code/data/model cards, and committee report;
- re-run literature/source/licence surveillance;
- define a precise novelty claim; and
- decide code, weights, and metadata release scope.

Deliverables:

- immutable external predictions;
- final validation report;
- manuscript-ready tables and figures;
- reproducibility bundle;
- negative-results appendix; and
- Gate 6 release/publication decision.

## 3. Proposed repository layout

After committee approval, scaffold the following structure:

    iScore3.0/
    ├── README.md
    ├── LICENSE
    ├── CITATION.cff
    ├── CONTRIBUTING.md
    ├── SECURITY.md
    ├── blueprint/                  # approved committee package
    │   ├── 00_executive_summary.md
    │   ├── ...
    │   ├── evidence/
    │   ├── figures/
    │   └── references/
    ├── docs/
    │   ├── architecture/
    │   │   └── decisions/          # numbered ADRs
    │   ├── data_cards/
    │   ├── model_cards/
    │   ├── protocols/
    │   └── operations/
    ├── src/
    │   └── iscore3/
    │       ├── chemistry/
    │       ├── protein/
    │       ├── pockets/
    │       ├── fusion/
    │       ├── models/
    │       ├── objectives/
    │       ├── uncertainty/
    │       ├── data/
    │       ├── evaluation/
    │       ├── inference/
    │       └── cli/
    ├── configs/
    │   ├── data/
    │   ├── features/
    │   ├── model/
    │   ├── experiment/
    │   ├── evaluation/
    │   └── cluster/
    ├── scripts/
    │   ├── acquire/
    │   ├── curate/
    │   ├── train/
    │   ├── evaluate/
    │   └── release/
    ├── workflows/
    │   ├── rules/
    │   ├── profiles/
    │   └── Snakefile
    ├── environments/
    │   ├── lock/
    │   └── containers/
    ├── tests/
    │   ├── unit/
    │   ├── integration/
    │   ├── scientific/
    │   ├── leakage/
    │   ├── invariance/
    │   └── fixtures/
    ├── data/
    │   ├── README.md
    │   ├── registry/
    │   ├── manifests/
    │   ├── raw/
    │   ├── external/
    │   ├── interim/
    │   ├── processed/
    │   ├── features/
    │   └── splits/
    ├── artifacts/
    │   ├── README.md
    │   ├── checkpoints/
    │   ├── predictions/
    │   ├── metrics/
    │   ├── logs/
    │   └── profiles/
    ├── reports/
    │   ├── generated/
    │   ├── figures/
    │   └── tables/
    ├── notebooks/
    │   ├── exploratory/
    │   └── archived/
    ├── third_party/
    │   ├── README.md
    │   └── notices/
    └── cluster/
        ├── slurm/
        └── logs/

Only blueprint/ exists during committee review. The remaining tree should be generated after Gate 1 so it reflects approved tooling and licences.

## 4. Directory contracts

### blueprint/

Committee-facing, human-readable decision package plus machine-readable evidence. After approval, changes require a dated amendment and decision record.

### data/

Never a miscellaneous drop folder. Each subdirectory has a README and registry entry.

- raw/: immutable source snapshots; restricted; never modified.
- external/: immutable locked benchmark inputs/labels with stronger access controls.
- interim/: restartable normalized tables and structure repairs.
- processed/: analysis-ready canonical observations and pockets.
- features/: versioned ligand/pocket artifacts, rebuildable from processed data.
- splits/: immutable identifiers and component/similarity audit files.
- manifests/: checksums, lineage, exclusions, schemas, licences, and build reports.

Raw, external, interim, processed, and feature payloads are ignored by Git. Manifests and permitted small metadata are versioned.

### artifacts/

Machine-generated experiment outputs. An experiment directory is immutable after completion and named by semantic label plus content/config hash, not “final” or “latest”. Large checkpoints live in approved artifact storage, referenced by checksum.

### reports/

Generated tables and figures only. Every item records the command/config and metric artifact that produced it. Manual figure editing occurs in a separate source file and is disclosed.

### notebooks/

Exploration only. A result used in a report must move to src/, scripts/, or workflows/ with tests. Notebooks use small samples and may not embed restricted records or secrets.

### third_party/

Not a place to copy arbitrary repositories. Record pinned source URL, revision, licence, local modifications, and citation. Prefer declared package dependencies or isolated submodules only after review.

## 5. Data registry

Each source receives a registry file containing:

- canonical name and role;
- release/version/date;
- official URL and citation;
- access method and credential owner;
- licence/terms and redistribution classification;
- expected files, sizes, and checksums;
- raw and processed schema versions;
- acquisition/build command;
- lineage dependencies;
- known issues;
- allowed train/validation/test uses; and
- deprecation/replacement history.

No acquisition script may write outside a source-specific validated directory. Download scripts stop on checksum mismatch and never overwrite an existing raw snapshot.

## 6. Content-addressed artifacts

Recommended identity:

    artifact_id = sha256(
        input_manifest_hash
        + feature_schema_hash
        + source_revision
        + resolved_configuration
        + checkpoint_hash_if_applicable
    )

Every model checkpoint records:

- code revision and dirty-worktree diff hash;
- environment/container digest;
- train/validation split hashes;
- feature and label schema hashes;
- seed and deterministic settings;
- resolved hyperparameters;
- hardware;
- parent/pretrained checkpoint hashes; and
- exact metrics/prediction artifact IDs.

A human alias may point to an artifact; it must not replace the immutable identity.

## 7. Configuration and workflow

Use declarative, composed configuration with schema validation. Separate:

- data identity and filters;
- feature schemas;
- model architecture;
- optimization;
- compute/cluster resources; and
- evaluation protocol.

Do not encode dataset paths or split logic inside model classes. Resolve absolute paths only at runtime from one project data-root setting.

A workflow engine such as Snakemake is recommended for acquisition-to-report dependency tracking and HPC profiles. Each rule writes to a temporary target and atomically promotes it only after validation. Failed partial outputs are never treated as complete.

## 8. Environments and dependencies

Maintain:

- a lightweight developer lock;
- a CUDA/GPU training lock;
- a CPU inference lock;
- a documented container digest for reproducibility; and
- a software bill of materials/licence report.

Pin Python, PyTorch, PyTorch Geometric/equivariant libraries, RDKit, Biopython/gemmi, pocket tools, and gMolAI revision. External executables such as DSSP, P2Rank, fpocket, MSMS, APBS, or protonation software have explicit wrappers, version checks, and licences.

Never let package import trigger a network download. Weights and databases are acquired by explicit, checksummed commands.

## 9. Testing strategy

### Unit tests

Parsing, unit conversion, qualifiers, canonicalization, graph construction, masks, losses, metrics, and artifact hashing.

### Integration tests

A tiny redistributable fixture runs receptor/site/SMILES through preprocessing, training smoke test, cached inference, and report generation.

### Scientific regression tests

Freeze small expected predictions/metrics with tolerances. These detect silent changes from RDKit, graph libraries, precision modes, or preprocessing.

### Leakage tests

Check duplicate identities, connected components, publication overlap, nearest similarities, feature lineage, and external-label access.

### Invariance tests

Rotation, translation, permutation, batching, padding, and deterministic caching.

### Property/fuzz tests

Malformed SMILES, unusual elements, mmCIF edge cases, alternate locations, insertion codes, missing atoms, multiple chains, qualifiers, and extreme values.

CI runs small CPU tests on every change; scheduled GPU and data-integration tests run in a controlled environment without exposing restricted data in logs.

## 10. Experiment tracking

A local/self-hosted tracker is preferred until data-governance review approves an external service. The tracker stores metadata and permitted scalar summaries; never raw structures, SMILES, unpublished labels, credentials, or restricted filenames by default.

Model selection uses a registry state machine:

- candidate;
- validation-frozen;
- confirmatory-frozen;
- externally-evaluated;
- rejected; or
- released.

State transitions require an owner, timestamp, evidence artifact, and reason.

## 11. Storage and compute planning

### Storage tiers

- Tier 0: Git — code, configs, manifests, small fixtures, documents.
- Tier 1: project filesystem — restricted raw/processed data and active features.
- Tier 2: artifact/object store — checkpoints, predictions, logs, archives.
- Tier 3: scratch — restartable caches with expiry and no sole-copy artifacts.

PLINDER alone is hundreds of gigabytes; full derived graphs and embeddings can exceed raw size. Estimate storage from a 1–5% pilot before full ingestion. Use columnar tables, compressed coordinate formats where permitted, and deduplicated per-pocket/per-ligand artifacts.

### Compute sequence

1. CPU preprocessing pilot;
2. one-GPU v0/v1 profiling;
3. learning-curve and batch-size estimate;
4. fixed baseline compute budget;
5. v1 ablation budget; and
6. extensions only after gates.

Record requested and consumed GPU-hours. Early stopping and pruning criteria use validation only. Do not launch broad sweeps before profiling data throughput.

## 12. Naming conventions

- source releases: source__release__rawhash
- processed views: source__view__schemahash
- splits: dataset__split-policy__splithash
- experiments: model__task__split__confighash__seed
- predictions: experiment-id__benchmark__predictionhash
- reports: protocol__dataset-version__model-family

Use lowercase ASCII, hyphens within human labels, and machine identifiers without spaces. Dates use ISO 8601. Do not use names such as new, final, final2, test, or best without an immutable pointer.

## 13. Branching and review

After Git initialization:

- protect the main branch;
- use short feature branches and reviewed pull requests;
- require tests, formatting, type/schema checks, and licence scan;
- prohibit committing raw/restricted data and secrets through hooks/CI;
- tag approved dataset schemas, model candidates, and releases;
- retain a CHANGELOG and migration notes; and
- never rewrite a released tag.

At least one reviewer independent of the author approves split/evaluation changes. External-test code changes after freeze require an incident record.

## 14. Roles

Minimum role separation:

| Role | Responsibility |
|---|---|
| scientific lead | scope, hypotheses, biological interpretation |
| data steward | rights, provenance, schema, mapping, immutable splits |
| model lead | architecture, training, unit/scientific tests |
| evaluation owner | locked labels, metrics, statistical analysis |
| infrastructure owner | environments, workflows, storage, compute |
| independent reviewer | leakage audit and gate evidence |

One person may hold multiple roles in a small team, but the evaluation-owner review should remain independent where possible.

## 15. Definition of done

A work package is done only when:

- code and configuration are reviewed;
- tests pass in the locked environment;
- inputs/outputs have immutable manifests;
- licences and citations are recorded;
- expected and failed records are counted;
- resource use is reported;
- results are regenerated by workflow, not notebook state;
- documentation and decision records are updated; and
- the next gate has an explicit approve/revise/stop decision.

## 16. Immediate actions after committee approval

1. initialize the approved repository scaffold and licence;
2. create ADR-0001 for the strict ligand-3D-free boundary;
3. resolve gMolAI licensing and pin its reference checkpoint;
4. approve/acquire the structural data track;
5. implement canonical schemas and synthetic fixtures;
6. build strict ligand and receptor/site preprocessing tests;
7. freeze first-pass union-component split logic; and
8. train only the v0 diagnostic baselines.

This order preserves the option to stop cheaply if the data rights, mappings, or structural signal are inadequate.
