# Gate-2 reproducibility record

The repository tracks processed tables, frozen configuration, splits, code,
predictions, metrics, audits, and manifests. Raw BindingDB/RCSB/AlphaFold data,
feature arrays, checkpoints, US-align source/binary, and the SIF are local and
ignored. All final inputs and outputs are hash-recorded in
`data/manifests/gate02-baseline-experiment-v1.json`.

The final run used the supplied ARM64 container on Arrhenius `n40`, plus the
existing pinned Biopython 1.85 local layer in
`environments/lock/gate01-arm64-requirements.txt`.

```text
containers/gmolai-pyg-25.09-arm64.sif
SHA-256 69cfad3e38de6397a94184a132a84d44bea8d1f8f258a1b80ba6a5dca5c714e2
```

Run from the repository root with the project mounted at `/work`:

```bash
apptainer exec --cleanenv \
  --env PYTHONPATH=/work/src:/work/environments/local/gate01 \
  --env PYTHONHASHSEED=0 --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind /nobackup/proj/disk/theo-storage/personal/jalil/iScore3.0:/work \
  --pwd /work containers/gmolai-pyg-25.09-arm64.sif \
  python -u scripts/gate02/run_baselines.py \
    --pilot data/processed/gate02/rcsb-kd-strict-v3.tsv \
    --pockets-s01 data/processed/gate02/rcsb-kd-pocket-views-strict-v3.tsv \
    --pockets-s2 data/processed/gate02/rcsb-kd-pocket-views-s2-v3.tsv \
    --pockets-s3 data/processed/gate02/rcsb-kd-pocket-views-s3-v1.tsv \
    --site-manifest data/manifests/gate02-rcsb-kd-sites-v3.json \
    --prefit-split data/splits/gate02/prefit-union-components-v3.tsv \
    --gmolai-manifest data/manifests/gate02-gmolai-v2-v3.json \
    --gmolai-feature-root data/features/gate02/gmolai-v2-v3 \
    --esm2-manifest data/manifests/gate02-esm2-v3.json \
    --esm2-feature-root data/features/gate02/esm2-v3 \
    --structural-allpairs data/processed/gate02/structural-similarity-allpairs-v3.tsv \
    --config configs/gate02/feasibility-effective-v3.yaml \
    --required-audit bindingdb_provenance=reports/gate02/evidence/bindingdb-provenance-audit-v3.json \
    --required-audit bindingdb_ties=reports/gate02/evidence/bindingdb-tie-equivalence-audit-v1.json \
    --required-audit prefit_components=reports/gate02/evidence/prefit-union-component-audit-v3.json \
    --required-audit structural_similarity=reports/gate02/evidence/structural-similarity-audit-v3.json \
    --required-audit gmolai_adapter=reports/gate02/evidence/gmolai-adapter-audit-v3.json \
    --required-audit esm2_adapter=reports/gate02/evidence/esm2-adapter-audit-v3.json \
    --required-audit apo_views=data/manifests/gate02-apo-structure-views-v1.json \
    --split-output data/splits/gate02/baseline-splits-v1.tsv \
    --leakage-output reports/gate02/evidence/leakage-diagnostics-final-v1.json \
    --prediction-output reports/gate02/evidence/baseline-predictions-v1.tsv \
    --hyperparameter-output reports/gate02/evidence/baseline-hyperparameters-v1.tsv \
    --metric-output reports/gate02/evidence/baseline-metrics-v1.json \
    --manifest-output data/manifests/gate02-baseline-experiment-v1.json
```

`immutable_write` refuses a differing overwrite. Use fresh output names for an
independent rerun. Expected result hashes are:

```text
a28bcabbed0f48bcc8443784c14a6a645df326ef7d53a4494ac6f0ff51c38e0a  data/splits/gate02/baseline-splits-v1.tsv
c1d14128ae4207c3af86c043f22706fac36b1d44e1f0c8de484f9f0d6ae8de39  reports/gate02/evidence/leakage-diagnostics-final-v1.json
9a2bfa8a4a8ae7f8d21dbd2f1de14b7d65c6bedb0f0491926e6bfdef5c3b1b7b  reports/gate02/evidence/baseline-predictions-v1.tsv
2db09262dbd595701cf35934f194c12ad1a143620aecbe877128a93b0c404e5c  reports/gate02/evidence/baseline-hyperparameters-v1.tsv
4d9447b0e9d188530c1cfbec8adc0f05d38460f6e6561facef9e03da9eb9c24e  reports/gate02/evidence/baseline-metrics-v1.json
```

The progression check must be `FAIL` /
`NO_GO_FOR_FULL_ARCHITECTURE_UNDER_FROZEN_GATE`; the primary tensor-versus-
additive delta is −0.002042 pKd with interval approximately
[−0.007448, +0.004358]. No full cross-attention model is part of this record.

Verification:

```bash
apptainer exec --cleanenv \
  --env PYTHONPATH=/work/src:/work/environments/local/gate01 \
  --bind /nobackup/proj/disk/theo-storage/personal/jalil/iScore3.0:/work \
  --pwd /work containers/gmolai-pyg-25.09-arm64.sif pytest -q

python -m compileall -q src scripts tests
git diff --check
```
