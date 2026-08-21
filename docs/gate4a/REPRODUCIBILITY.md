# Gate-4A admission reproducibility

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
