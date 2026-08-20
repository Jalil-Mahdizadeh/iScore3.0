# iScore3.0 Research Blueprint

**Status:** Gate-0/1 completed; full-model progression not recommended
**Evidence cut-off:** 20 August 2026

This package defines a staged research programme for iScore3.0: a pose-free protein–small-molecule affinity model that accepts a three-dimensional receptor binding site and a ligand represented only by SMILES-derived 1D/2D information.

## Central conclusion

The broad concept is feasible but is not unprecedented. PSG-BAR, PLMCA, AttentionMGT-DTA, BANANA, CASTER-DTA, 3DProtDTA, HoloProt, Graph_RG, and PLANET are important prior art. The defensible research opportunity is a strictly ligand-3D-free and provenance-first system that can *demonstrate* incremental geometry value beyond ligand, sequence, metadata, and nearest-neighbour controls under conformation and leakage stress tests.

The bounded Gate-0/1 pilot did not meet that progression criterion. The current recommendation is to remediate data independence, structural-similarity auditing, conformation coverage, licensing, and direct baselines before training the proposed full architecture. See the [Gate-0/1 report](../reports/gate01/GATE_0_1_REPORT.md).

Data leakage is treated as a release-blocking scientific threat. The package does not label any one split “leak-free”; it requires exact-contamination prevention, union-component and two-dimensional splitting, continuous nearest-neighbour audits, deployment-specific novelty tiers, ligand-/protein-/metadata-only red teams, and locked external evaluation.

## Committee reading order

1. [00_executive_summary.md](00_executive_summary.md) — recommendation, novelty boundary, and approval gates
2. [01_scope_and_requirements.md](01_scope_and_requirements.md) — exact task, information boundary, intended use, and exclusions
3. [08_data_leakage_threat_model.md](08_data_leakage_threat_model.md) — dedicated leakage debate, threat model, audits, and release gates
4. [02_literature_and_software_review.md](02_literature_and_software_review.md) — publication and source-code synthesis
5. [03_dataset_strategy.md](03_dataset_strategy.md) — dataset inventory, licences, mappings, curation, and splits
6. [04_technical_blueprint.md](04_technical_blueprint.md) — pocket encoding, ligand adapter, fusion, objectives, and inference
7. [05_evaluation_and_validation.md](05_evaluation_and_validation.md) — hypotheses, controls, metrics, statistics, and acceptance criteria
8. [06_execution_and_repository_plan.md](06_execution_and_repository_plan.md) — gated work plan and scalable repository layout
9. [07_risks_decisions_and_governance.md](07_risks_decisions_and_governance.md) — risk register, stop rules, and committee decisions

## Evidence package

- [evidence/publications.tsv](evidence/publications.tsv) — 82 evaluated publications/methods/resources with coordinate-use classification
- [evidence/software.tsv](evidence/software.tsv) — 47 source-code, input-interface, reproducibility, and licence audits
- [evidence/datasets.tsv](evidence/datasets.tsv) — 47 dataset, benchmark, structure, and mapping records with versions, roles, access, and caveats
- [evidence/leakage_threats.tsv](evidence/leakage_threats.tsv) — 26 leakage/shortcut channels and release rules
- [evidence/search_log.md](evidence/search_log.md) — search waves, example queries, review protocol, negative searches, and limitations
- [references/references.bib](references/references.bib) — 77 core bibliographic records
- [figures/system_architecture.mmd](figures/system_architecture.mmd) — editable system architecture
- [figures/data_and_evaluation_flow.mmd](figures/data_and_evaluation_flow.mmd) — editable data/evaluation flow
- [figures/leakage_audit_flow.mmd](figures/leakage_audit_flow.mmd) — editable leakage-audit flow

## Proposed first model

The recommended v1 protein encoder is a residue-level GVP-style graph plus a local receptor heavy-atom graph. gMolAI supplies aligned atom tokens and its calibrated global molecular representation. Bidirectional cross-attention learns pose-free compatibility; protein self-updates may use receptor geometry and ligand self-updates may use graph topology, but cross-modal layers receive no receptor–ligand distances.

Simple ligand-only, protein-only, nearest-neighbour, late-fusion, 3DProtDTA-style, sequence-only, and BANANA baselines are completed first. A surface branch, pharmacophore field, and pair tensor remain gated extensions.

## Data position

The recommended programme combines a rights-cleared PDBbind structural core with official LP-PDBbind and CleanSplit protocols, BindingDB/ChEMBL or BigBind for scale, and PLINDER/BioLiP2 for pockets and receptor conformations. PLINDER’s currently disabled/incorrect affinity mapping is prohibited. Unmodified LIT-PCBA is downgraded to a red-team diagnostic because a reproducible 2025 preprint reports exact and analogue leakage; any confirmatory derivative must be independently rebuilt and audited.

No third-party dataset, paper, checkpoint, or repository is copied into this package. Future acquisition uses source-specific registries, rights records, checksums, immutable manifests, and controlled external storage.

## Decision needed

The committee should record approve/modify/reject decisions for the strict ligand-3D-free rule, primary supplied-site task, endpoint policy, version-1 domain, data licences, gMolAI licence, pocket encoder, leakage protocol, practical success threshold, external-label control, and release intent. The checklist is in 07_risks_decisions_and_governance.md.
