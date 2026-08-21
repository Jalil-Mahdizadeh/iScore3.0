# Gate-4A representation preregistration

This document contains the initial plan. The outcome-generating freeze is
[`DELTA3D_LIGAND_PROTOCOL.md`](DELTA3D_LIGAND_PROTOCOL.md). Before feature
extraction, the tentative Uni-Mol2 choice below was replaced by the mature,
explicit-coordinate Uni-Mol v1 molecule-all-H checkpoint because Uni-Mol2 had
unresolved interface/batching reproducibility concerns.

## gMolAI 2D

Use the existing audited gMolAI adapter and its canonical atom mapping. Global
and atom-level outputs remain frozen. Any new checkpoint must receive a separate
hash/provenance audit before use.

## Deterministic free conformers

Input is canonical isomeric SMILES only. The initial configuration is:

- RDKit ETKDGv3, fixed seed `20260821`;
- 32 attempted conformers;
- 0.5 angstrom RMS pruning;
- MMFF94s optimization, single thread, at most 500 iterations;
- retain converged conformers within 10 kcal/mol of the minimum;
- Boltzmann ensemble summaries at 298.15 K;
- no crystallographic or docked ligand coordinates for generation, filtering,
  selection, scoring, supervision, or quality control.

Descriptors are rotation/translation invariant shape summaries, intramolecular
pharmacophore-pair distance summaries, conformer energy statistics, and internal
conformer diversity. Feature presence/count information is labelled separately
because it is largely 2D chemistry.

The initial aggregator is deterministic weighted mean/standard deviation. No
learned conformer attention is permitted.

## Frozen pretrained ligand 3D

The executed model is pretraining-only Uni-Mol v1 molecule-all-H from
`unimol_tools` v0.1.6 at commit `4596596a`, checkpoint
`mol_pre_all_h_220816.pt` with SHA-256 `7f5f14bb...3446a7`. It consumes the exact
free-conformer atom/coordinate arrays, is never affinity-fine-tuned, and aggregates
per-conformer 512-vectors by frozen Boltzmann mean and standard deviation.

Actual-coordinate embeddings require an identical-encoder geometry-ablation
control. A gain unique to the pretrained embedding, without deterministic or
coordinate-ablation corroboration, is not evidence for ligand geometry.

## Matched interaction projection and budget

Each frozen ligand-2D, ligand-3D, and pocket vector is standardized and projected
to 32 dimensions by training-fold-only, label-free PCA. PCA is unwhitened; its sign
is fixed by making each component's largest-absolute loading positive. A fold with
rank below 32 is invalid rather than silently changing capacity.

The reference interaction has one 2D-by-pocket rank-8 bilinear factor:
`2 × 32 × 8 = 512` trainable factor parameters. The augmented model has separate
rank-4 2D-by-pocket and rank-4 3D-by-pocket factors:
`2 × 32 × 4 + 2 × 32 × 4 = 512`. Factor biases are forbidden, projections are
non-trainable, additive heads are identical, and tuning budgets are paired.

The final five-seed build generated all 69 ligands per seed, retained 1–32
conformers (median 16–18), and produced 69 unique geometry hashes per seed. Both
deterministic and pretrained 3D branches failed the preregistered progression
criteria; no trainable representation is authorized.
