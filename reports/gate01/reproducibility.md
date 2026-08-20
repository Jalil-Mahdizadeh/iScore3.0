# Gate-0/1 Reproducibility Record

## Reproduction level

The tracked repository contains code, configs, small processed tables, splits, predictions, metrics, audit reports, and immutable manifests. Large/mutable artifacts remain local and ignored:

- 11.1 GB Apptainer SIF;
- RCSB API snapshots and 86 mmCIF files;
- BindingDB 202608 ZIP and 8.98 GB extracted TSV;
- gMolAI checkpoint/source mirror and generated arrays; and
- competitor source mirrors.

Their official URLs, revisions, sizes, and SHA-256 hashes are tracked. Exact historical raw-byte reproduction still requires access to the preserved local files or an archive serving the same bytes; a future live API response is not assumed identical.

## Runtime

The run was executed directly on Arrhenius node `n40`, an ARM64 NVIDIA GH200 node, with Apptainer 1.5.3. The canonical environment record is [`evidence/environment-v1.json`](evidence/environment-v1.json).

Container:

```text
containers/gmolai-pyg-25.09-arm64.sif
SHA-256 69cfad3e38de6397a94184a132a84d44bea8d1f8f258a1b80ba6a5dca5c714e2
```

No extra SIF was needed. Gemmi 0.7.5 and Biopython 1.85 were installed as a local, ignored layer from [`environments/lock/gate01-arm64-requirements.txt`](../../environments/lock/gate01-arm64-requirements.txt). The exact container invocation used throughout was:

```bash
apptainer exec --nv \
  --bind /nobackup/proj/disk/theo-storage/personal/jalil/iScore3.0:/work \
  --pwd /work \
  containers/gmolai-pyg-25.09-arm64.sif \
  bash -lc 'PYTHONPATH=/work/environments/local/gate01:/work/src <COMMAND>'
```

Replace `<COMMAND>` below while retaining the prefix.

## Third-party source reconstruction

The source trees live under ignored `third_party/source_cache/`. For gMolAI:

```bash
git clone https://github.com/Jalil-Mahdizadeh/gMolAI-v2.0.git third_party/source_cache/gmolai-v2.0
git -C third_party/source_cache/gmolai-v2.0 checkout 29ced8352886692c00a792d24a183f05e5de0059
sha256sum -c <independently generated list from configs/gate01/gmolai_adapter.yaml>
```

The adapter refuses a revision, source file, checkpoint, calibrator, schema, or resolved-config hash mismatch. Competitor revisions and inspected-file hashes are in [`third_party/reviewed_sources.tsv`](../../third_party/reviewed_sources.tsv).

## From-source data workflow

### 1. Acquire and select the RCSB candidate pilot

These commands query a live API and therefore belong in a fresh versioned raw directory; compare every returned hash to the 2026-08-20 acquisition manifest.

```bash
python scripts/gate01/rcsb_pilot.py acquire-metadata \
  --raw-root data/raw/rcsb/2026-08-20 \
  --manifest data/manifests/rcsb-kd-2026-08-20-acquisition.json \
  --endpoint Kd --batch-size 100

python scripts/gate01/rcsb_pilot.py select \
  --raw-root data/raw/rcsb/2026-08-20 \
  --output data/processed/gate01/rcsb-kd-pilot-v1.tsv \
  --manifest data/manifests/rcsb-kd-2026-08-20-selection-v1.json \
  --min-supervised-per-construct 8 \
  --max-supervised-per-construct 20 \
  --max-constructs 12 \
  --replicate-tolerance-pkd 0.30 \
  --selection-seed 20260820

python scripts/gate01/rcsb_pilot.py acquire-structures \
  --selection data/processed/gate01/rcsb-kd-pilot-v1.tsv \
  --coordinate-root data/raw/rcsb/2026-08-20/structures \
  --manifest data/manifests/rcsb-kd-2026-08-20-structures-v1.json
```

Expected candidate-table SHA-256: `c60d291f538ea9b59df5e0b3c41ba854f6ea667ccfc7611972164e69edea49ae`.

### 2. Build S0/S1 candidate pocket views

```bash
python scripts/gate01/build_pockets.py \
  --pilot data/processed/gate01/rcsb-kd-pilot-v1.tsv \
  --coordinate-root data/raw/rcsb/2026-08-20/structures \
  --output data/processed/gate01/rcsb-kd-pocket-views-v1.tsv \
  --site-manifest data/manifests/rcsb-kd-2026-08-20-sites-v1.json \
  --build-manifest data/manifests/rcsb-kd-2026-08-20-pockets-v1.json \
  --cutoff-angstrom 6.0 --minimum-coverage 0.80
```

Expected candidate pocket-table SHA-256: `8fc078d03294971c060c77ea4dcba911cc862fa6b88147e5c77e259623e2727e`.

### 3. Audit and encode gMolAI

```bash
python scripts/gate01/audit_gmolai.py \
  --pilot data/processed/gate01/rcsb-kd-pilot-v1.tsv \
  --source-root third_party/source_cache/gmolai-v2.0 \
  --adapter-config configs/gate01/gmolai_adapter.yaml \
  --container containers/gmolai-pyg-25.09-arm64.sif \
  --feature-root data/features/gate01/gmolai-v2-pilot-v1 \
  --manifest data/manifests/gmolai-v2-pilot-v1.json \
  --audit-report reports/gate01/evidence/gmolai-adapter-audit-v1.json \
  --device cuda
```

Expected manifest SHA-256: `8f74da3fe5a800887e893f222a023de7c086c5d70ddafbcbdc776e3fc4e2200a`. The audit must report all checks true. CPU/GPU outputs are tolerance-equal, not bitwise identical.

### 4. Acquire and reconcile BindingDB

Official archive used:

```text
https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202608_tsv.zip
SHA-256 8dc30541d668e5e403ba1fe9a682b76c42a50e1711aac44329c62ba24ccd57db
bytes 592986963
```

After downloading that exact archive:

```bash
python scripts/gate01/audit_bindingdb.py \
  --pilot data/processed/gate01/rcsb-kd-pilot-v1.tsv \
  --archive data/raw/bindingdb/202608/BindingDB_All_202608_tsv.zip \
  --output data/processed/gate01/bindingdb-202608-provenance-v2.tsv \
  --report reports/gate01/evidence/bindingdb-provenance-audit-v2.json \
  --source-url https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_202608_tsv.zip \
  --expected-archive-sha256 8dc30541d668e5e403ba1fe9a682b76c42a50e1711aac44329c62ba24ccd57db \
  --pockets data/processed/gate01/rcsb-kd-pocket-views-v1.tsv \
  --strict-pilot data/processed/gate01/rcsb-kd-pilot-strict-v2.tsv \
  --strict-pockets data/processed/gate01/rcsb-kd-pocket-views-strict-v2.tsv \
  --min-supervised-per-construct 8
```

Expected canonical hashes:

```text
4be3ec4bb841987249c633e2547b5bfb323d41982f2f7572031baaa11d3dca40  bindingdb-202608-provenance-v2.tsv
a51d2174d9d733a5f58d32e80908134d32da8d5e8a52ca7d22a26c5b3e22d532  rcsb-kd-pilot-strict-v2.tsv
1421234dbe9a6cb337081fa7fd74aba96a8ab381b05049fb65417475dfcff996  rcsb-kd-pocket-views-strict-v2.tsv
```

The audit’s top-level status is expected to be `FAIL` because it documents rejection of the candidate table; `strict_dataset.status` must be `PASS`.

### 5. Run the strict shallow controls

```bash
python scripts/gate01/run_baselines.py \
  --pilot data/processed/gate01/rcsb-kd-pilot-strict-v2.tsv \
  --pockets data/processed/gate01/rcsb-kd-pocket-views-strict-v2.tsv \
  --sites data/manifests/rcsb-kd-2026-08-20-sites-v1.json \
  --gmolai-manifest data/manifests/gmolai-v2-pilot-v1.json \
  --gmolai-feature-root data/features/gate01/gmolai-v2-pilot-v1 \
  --config configs/gate01/baselines.yaml \
  --split-output data/splits/gate01/pilot-splits-strict-v2.tsv \
  --leakage-output reports/gate01/evidence/leakage-diagnostics-strict-v2.json \
  --prediction-output reports/gate01/evidence/baseline-predictions-strict-v2.tsv \
  --hyperparameter-output reports/gate01/evidence/baseline-hyperparameters-strict-v2.tsv \
  --metric-output reports/gate01/evidence/baseline-metrics-strict-v2.json \
  --manifest-output data/manifests/gate01-baselines-strict-v2.json
```

The frozen config predates BindingDB reconciliation and its narrative `only_six_independent_components` string is a pre-run estimate. It was not altered after fitting because doing so would invalidate its recorded SHA-256. The executable result and manifest correctly contain five components; the statistical criterion itself is unchanged.

Expected result hashes are recorded inside [`data/manifests/gate01-baselines-strict-v2.json`](../../data/manifests/gate01-baselines-strict-v2.json). The progression check must return `FAIL`, point delta `-0.0399876464`, and interval `[-0.1517929649, 0.0796584652]`.

## Verification suite

Run from the same container:

```bash
python -m pytest -q
python -m compileall -q src scripts tests
git diff --check
```

The repository includes unit, integration, and information-boundary tests for atom mapping, adapter rejection behaviour, deterministic embeddings, pocket geometry invariance, historical-label quarantine, and BindingDB exact/qualified matching. The upstream gMolAI suite is also run separately at its pinned revision.

Final observed result: **18/18 iScore3.0 tests passed** and **67/67 upstream gMolAI tests passed** in the pinned SIF. `compileall` passed; Ruff was not installed in the image and is therefore reported as not run rather than silently assumed.

## Canonical versus superseded artifacts

Only filenames containing `strict-v2` and `bindingdb-provenance-audit-v2.json` support the Gate decision. `*-v1` baseline and provenance outputs document the pre-BindingDB-reconciliation run and are retained to make the correction visible. They must not be pooled with strict-v2 results.
