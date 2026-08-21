# Gate-4A reproducibility

## Completed isolated Delta3D-ligand experiment

The outcome-generating protocol and code were sealed in pre-outcome commit
`55c176f`. Runtime provenance is in
`reports/gate4a/evidence/environment-delta3d-v2.json`; the feature, raw-ensemble,
OOF and result hashes are in the three tracked `delta3d-*-v1.json` evidence files.

Use a disposable replay copy with ignored source/checkpoint/raw snapshots restored
at their manifest hashes. Move the three tracked output evidence JSONs aside in
that disposable copy before running: feature extraction, evaluation and validation
intentionally refuse to overwrite them. Install the small Uni-Mol overlay and run:

```bash
apptainer exec containers/gmolai-pyg-25.09-arm64.sif \
  python -m pip install --target environments/local/gate4a-delta3d/vendor \
  -r environments/gate4a-delta3d-requirements.txt

DELTA3D_PY=(apptainer exec --nv containers/gmolai-pyg-25.09-arm64.sif env \
  PYTHONPATH=src:environments/local/gate4a/vendor:environments/local/gate4a-delta3d/vendor:third_party/source_cache/unimol_tools-v0.1.6 \
  CUBLAS_WORKSPACE_CONFIG=:4096:8 python)

"${DELTA3D_PY[@]}" scripts/gate4a/extract_delta3d_ligand_features.py --device cuda

apptainer exec containers/gmolai-pyg-25.09-arm64.sif env \
  PYTHONPATH=src:environments/local/gate4a/vendor \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python scripts/gate4a/evaluate_delta3d_ligand.py

apptainer exec containers/gmolai-pyg-25.09-arm64.sif env PYTHONPATH=src \
  python scripts/gate4a/validate_delta3d_results.py
```

Feature extraction is label/receptor blind. Evaluation verifies the feature hash
before reading Davis labels and never fits receptor or interaction models. Compare
per-array feature hashes (the NPZ container timestamp itself need not be byte-
identical), model metrics, component-bootstrap summaries and reconstructed
likelihoods with the tracked evidence.

The complete test suite is:

```bash
apptainer exec --nv containers/gmolai-pyg-25.09-arm64.sif env \
  PYTHONPATH=src:environments/local/gate4a/vendor:environments/local/gate4a-delta3d/vendor:third_party/source_cache/unimol_tools-v0.1.6 \
  ISCORE3_TEST_DEVICE=cuda pytest -q
```

It passed 50/50 tests on 2026-08-21.

## Admission and provenance-closure history

This phase performs data admission and label-blind input audits only. It trains no
predictive model. The recorded environment is
`reports/gate4a/evidence/environment-v1.json`; the SIF itself is intentionally
Git-ignored, but its size and SHA-256 are frozen there. Raw publisher/API payloads
are also ignored and are admitted only through tracked manifests containing their
official URLs, acquisition times, byte counts, and hashes.

Define the common Bash command array:

```bash
GATE4A_PY=(apptainer exec containers/gmolai-pyg-25.09-arm64.sif \
  env PYTHONPATH=src:environments/local/gate4a/vendor python)
```

Then run the admission scripts into an isolated temporary directory:

```bash
"${GATE4A_PY[@]}" scripts/gate4a/adjudicate_davis_compounds.py \
  --output /tmp/iscore3-g4a-compounds-replay.tsv \
  --audit-output /tmp/iscore3-g4a-compounds-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_receptor_admission.py \
  --output /tmp/iscore3-g4a-receptors-replay.tsv \
  --audit-output /tmp/iscore3-g4a-receptors-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_apo_candidates.py \
  --output /tmp/iscore3-g4a-apo-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_okl_confirmation.py \
  --output /tmp/iscore3-g4a-okl-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_kirhub_confirmation.py \
  --output /tmp/iscore3-g4a-kirhub-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_admitted_design.py \
  --components-output /tmp/iscore3-g4a-components-replay.tsv \
  --audit-output /tmp/iscore3-g4a-design-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_admitted_3d_informativeness.py \
  --feature-output /tmp/iscore3-g4a-admitted-conformer-replay-v1.json \
  --logical-feature-path data/features/gate4a/davis-admitted-conformer-audit-v1.json \
  --audit-output /tmp/iscore3-g4a-admitted-conformer-replay-audit-v1.json
```

Compare each temporary output with its tracked counterpart using `cmp`; every
listed audit was byte-identical on 2026-08-21. Frozen hashes
are in `reports/gate4a/evidence/admission-reproducibility-v1.json` and
`reports/gate4a/evidence/conformer-reproducibility-v1.json`.

Acquisition scripts are deliberately separate because they make network calls and
refuse to overwrite immutable raw evidence. Their tracked manifests cover Davis,
PubChem, ChEMBL, KLIFS/UniProt/AlphaFold, RCSB, OKL, and KIRHub sources. A replay
must first verify those hashes; a changed upstream response is a new source
version, never an in-place update.

## Provenance-closure replay

The closure added Biopython 1.85 and Gemmi 0.7.3 to the hash-locked overlay. Install
the lock with `--require-hashes`; do not use an unpinned host parser. Raw structure,
KLIFS, RCSB and IDG-DREAM payloads are immutable and ignored by Git. Verify them
against the tracked coordinate/candidate/noise manifests before processing.

The label-blind processing commands are:

```bash
"${GATE4A_PY[@]}" scripts/gate4a/freeze_ligand_identity_qa.py \
  --qa-packet /tmp/iscore3-g4a-ligand-qa-replay.tsv \
  --freeze /tmp/iscore3-g4a-ligand-freeze-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/close_alphafold_receptors.py \
  --derived-root /tmp/iscore3-g4a-af-pockets-replay \
  --ledger /tmp/iscore3-g4a-af-ledger-replay.tsv \
  --manifest /tmp/iscore3-g4a-af-manifest-replay.json \
  --audit /tmp/iscore3-g4a-af-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/qualify_apo_tiers.py \
  --ledger /tmp/iscore3-g4a-apo-ledger-replay.tsv \
  --manifest /tmp/iscore3-g4a-apo-manifest-replay.json \
  --audit /tmp/iscore3-g4a-apo-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/freeze_confirmation_ledgers.py \
  --output-dir /tmp/iscore3-g4a-confirmation-replay \
  --audit /tmp/iscore3-g4a-confirmation-audit-replay.json
"${GATE4A_PY[@]}" scripts/gate4a/audit_practical_equivalence_evidence.py \
  --output /tmp/iscore3-g4a-equivalence-replay.json
```

The 41,328-pair US-align replay is deterministic but computationally longer. Run it
only after the coordinate ledger replay passes:

```bash
"${GATE4A_PY[@]}" scripts/gate4a/finalize_structural_leakage.py \
  --coordinates /tmp/iscore3-g4a-af-ledger-replay.tsv \
  --all-pairs /tmp/iscore3-g4a-structure-pairs-replay.tsv \
  --edges /tmp/iscore3-g4a-structure-edges-replay.tsv \
  --components /tmp/iscore3-g4a-components-final-replay.tsv \
  --audit /tmp/iscore3-g4a-structure-audit-replay.json
```

US-align is pinned as version 20260819 with binary SHA-256
`b5b44c885e61bba1352a20fd8726a4d73f1bd9798e13dc926e315c72b550e4d3`.
Tracked output hashes are recorded in each evidence JSON and collectively in
`reports/gate4a/evidence/provenance-closure-reproducibility-v1.json`.

The clean replay reproduced all parsed US-align metrics and both the edge and
component ledgers byte-for-byte. The per-row `stdout_sha256` alone is mount-path
dependent because US-align prints its absolute input filenames (`/workspace` in the
SIF versus the host path). It is raw-execution provenance, not a scientific field.
This expected distinction and every other replay comparison are frozen in
`reports/gate4a/evidence/provenance-closure-replay-v1.json`.

Network acquisition is intentionally separate:

```bash
"${GATE4A_PY[@]}" scripts/gate4a/acquire_site_unoccupied_candidates.py
"${GATE4A_PY[@]}" scripts/gate4a/acquire_idg_dream_noise_evidence.py
```

These commands refuse to overwrite the existing snapshot. Reacquisition belongs in
a new versioned raw root/manifest. No command above trains or evaluates a predictive
model.
