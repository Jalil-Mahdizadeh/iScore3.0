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
