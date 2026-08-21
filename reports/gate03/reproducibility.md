# Gate-3 reproducibility record

The final Gate-3 run used the supplied ARM64 Apptainer image on Arrhenius `n40`:

```text
containers/gmolai-pyg-25.09-arm64.sif
SHA-256 69cfad3e38de6397a94184a132a84d44bea8d1f8f258a1b80ba6a5dca5c714e2
```

The evaluator records prefit contract commit `e919fe3` and final outcome-blind
evaluator/amendment commit `9b16414`. It refuses tracked worktree changes at run
start, verifies every frozen input/audit hash, uses one numerical thread, and
enforces deterministic Torch algorithms. Result writes are immutable.

Run from the repository root:

```bash
apptainer exec --cleanenv \
  --env PYTHONPATH=/work/src:/work/environments/local/gate01 \
  --env PYTHONHASHSEED=0 --env OMP_NUM_THREADS=1 \
  --env OPENBLAS_NUM_THREADS=1 --env MKL_NUM_THREADS=1 \
  --bind /nobackup/proj/disk/theo-storage/personal/jalil/iScore3.0:/work \
  --pwd /work containers/gmolai-pyg-25.09-arm64.sif \
  python -u scripts/gate03/run_evaluation.py \
    --dataset data/processed/gate03/gate03-strict-kd-v3.tsv \
    --splits data/splits/gate03/gate03-component-scaffold-splits-v3.tsv \
    --pocket data/processed/gate03/receptor-pocket-v2-features-v1.tsv \
    --pocket data/processed/gate03/receptor-pocket-v2-features-s3-v1.tsv \
    --esm2-manifest data/manifests/gate03-esm2-v1.json \
    --esm2-root data/features/gate03/esm2-v1 \
    --esm-if1-manifest data/manifests/gate03-esm-if1-v2.json \
    --esm-if1-root data/features/gate03/esm-if1-v2 \
    --gmolai-manifest data/manifests/gate03-gmolai-v2-v1.json \
    --gmolai-root data/features/gate03/gmolai-v2-v1 \
    --structural-allpairs data/processed/gate03/structural-similarity-allpairs-v2.tsv \
    --config configs/gate03/evaluation-effective-v1.yaml \
    --amendment configs/gate03/protocol-amendment-04.yaml \
    --amendment configs/gate03/protocol-amendment-05.yaml \
    --required-audit reports/gate03/evidence/dataset-gate-audit-v3.json \
    --required-audit reports/gate03/evidence/primary-source-verification-v4.json \
    --required-audit reports/gate03/evidence/receptor-view-audit-v1.json \
    --required-audit reports/gate03/evidence/apo-view-audit-v1.json \
    --required-audit reports/gate03/evidence/esm2-adapter-audit-v1.json \
    --required-audit reports/gate03/evidence/esm-if1-adapter-audit-v2.json \
    --required-audit reports/gate03/evidence/gmolai-adapter-audit-v1.json \
    --required-audit reports/gate03/evidence/structural-similarity-audit-v2.json \
    --split-output data/splits/gate03/gate03-evaluation-folds-v1.tsv \
    --prediction-output reports/gate03/evidence/evaluation-predictions-v1.tsv \
    --hyperparameter-output reports/gate03/evidence/evaluation-hyperparameters-v1.tsv \
    --leakage-output reports/gate03/evidence/leakage-diagnostics-v1.json \
    --metric-output reports/gate03/evidence/evaluation-metrics-v1.json \
    --manifest-output data/manifests/gate03-evaluation-v1.json
```

Canonical output hashes:

```text
3e1bc568fc6e3854c1325825072536c4170575e8c8393145e3b52454a5c6458a  data/splits/gate03/gate03-evaluation-folds-v1.tsv
17444d41e69982192535bc96239c6e682113adbcd01fbbeec21387077fa77180  reports/gate03/evidence/evaluation-predictions-v1.tsv
34ad668c102041b7fe35a1e9cbe802a483198a3b19da8be989740aa270ea6593  reports/gate03/evidence/evaluation-hyperparameters-v1.tsv
5b61f39c1186aa291cf8c9a57f2576178c022caee2b4e18da3a9a6aada78a224  reports/gate03/evidence/leakage-diagnostics-v1.json
9700a3e4cb0b1d85058ee92739c929bc5aa48fcee688735bcd48f5c6858001a5  reports/gate03/evidence/evaluation-metrics-v1.json
38c9aa27cc1f9970d5b1f08228b8f2665d456c590e531d4a6c15cbf4207b8cfe  data/manifests/gate03-evaluation-v1.json
```

The canonical run must report `HYPOTHESIS-NO-GO`, 50,870 predictions, 175 metric
summaries, 46 paired comparisons, and leakage status `PASS`. The manifest records
3,550 fitted fold/model records, 661 observations, 34 series, runtime package
versions, and the information-boundary assertions.

Verification without any new efficacy fit:

```bash
apptainer exec --cleanenv \
  --env PYTHONPATH=/work/src:/work/environments/local/gate01 \
  --bind /nobackup/proj/disk/theo-storage/personal/jalil/iScore3.0:/work \
  --pwd /work containers/gmolai-pyg-25.09-arm64.sif pytest -q

python -m compileall -q src scripts tests
git diff --check
sha256sum -c reports/gate03/output-sha256.txt
```

No full cross-attention model was trained. Do not delete, overwrite, or reinterpret
these outputs as a tuning set for another model under the same hypothesis.
