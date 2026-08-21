# Gate-4A representation preregistration

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

The planned model family is a pretraining-only Uni-Mol2 molecular encoder. Model
size, exact checkpoint URL, source commit, licence, file hash, input convention,
coordinate normalization, and pooling are **not yet frozen**. No embeddings may
be generated until those fields are completed and it is verified that the
checkpoint is not docking- or affinity-fine-tuned.

Actual-coordinate embeddings require an identical-encoder geometry-ablation
control. A gain unique to the pretrained embedding, without deterministic or
coordinate-ablation corroboration, is not evidence for ligand geometry.
