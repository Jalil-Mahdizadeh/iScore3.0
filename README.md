# iScore3.0

iScore3.0 is a research programme for predicting protein–small-molecule binding affinity from a receptor structure and a ligand SMILES string, without using or generating query-ligand 3D coordinates or a docking pose.

The repository contains completed **Gate-0/1 and bounded Gate-2 feasibility
phases**. Gate-2 expanded the strict provenance-first dataset, added structural
leakage edges and predicted/unoccupied receptor views, and tested a low-capacity
ligand–pocket interaction control. It did not train the proposed full architecture.

## Start here

- [`blueprint/`](blueprint/) — committee blueprint and literature/data/leakage evidence.
- [`reports/gate01/`](reports/gate01/) — Gate-0/1 report and reproducibility record.
- [`reports/gate02/`](reports/gate02/) — Gate-2 no-go decision, component/OOD results, and reproducibility record.
- [`configs/gate01/`](configs/gate01/) — frozen pilot and baseline configurations.
- [`configs/gate02/`](configs/gate02/) — frozen interaction-feasibility protocol and amendments.
- [`src/iscore3/`](src/iscore3/) — tested implementation.
- [`data/README.md`](data/README.md) — data contracts and local rebuild instructions.

Large source data, external repositories, container images, features, and checkpoints are intentionally not committed. Their immutable metadata and checksums are versioned instead.

## Scientific boundary

The strict track permits receptor coordinates and ligand connectivity/features derived directly from SMILES. It forbids query-ligand conformers, docked/co-crystal poses, query-derived contacts/distances, and pretrained ligand features whose provenance cannot rule out unintended 3D or label exposure. See [`blueprint/01_scope_and_requirements.md`](blueprint/01_scope_and_requirements.md) and [`blueprint/08_data_leakage_threat_model.md`](blueprint/08_data_leakage_threat_model.md).

## Reproducing Gate-0/1

Use [`reports/gate01/reproducibility.md`](reports/gate01/reproducibility.md). The supplied Apptainer image is treated as a local runtime artifact and identified by SHA-256; it is not stored in Git.

## Status

Gate-0/1 and Gate-2 are complete. Gate-2 is **no-go for full architecture
training**: the predeclared low-capacity interaction term did not add reproducible
signal beyond matched additive fusion under strict component/OOD evaluation. See
the [`Gate-2 report`](reports/gate02/GATE_2_REPORT.md). No model is suitable for
clinical, regulatory, or prospective drug-discovery decisions.
