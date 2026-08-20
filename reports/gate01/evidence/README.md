# Gate-0/1 evidence map

## Canonical decision evidence

- `gmolai-adapter-audit-v1.json` — adapter identity, atom mapping, deterministic and cross-device checks on the original 80 candidate ligands
- `bindingdb-provenance-audit-v2.json` — candidate failure and strict-v2 reconciliation
- `leakage-diagnostics-strict-v2.json` — union-edge thresholds, component composition, similarity maxima, and unresolved structure edge
- `baseline-metrics-strict-v2.json` — primary/diagnostic metrics, paired component bootstrap, and progression decision
- `baseline-predictions-strict-v2.tsv` — all 3,050 out-of-fold predictions
- `baseline-hyperparameters-strict-v2.tsv` — fold-local model selection record
- `environment-v1.json` — HPC/GPU/container/runtime identity

Canonical data/split manifests are `data/manifests/gate01-pilot-strict-v2.json` and `data/manifests/gate01-baselines-strict-v2.json`.

## Deliberately retained superseded evidence

Files ending in `-v1` for BindingDB provenance, leakage, metrics, predictions, hyperparameters, and splits were generated before the BindingDB qualifier reconciliation. They document why the analysis changed but must not support scientific conclusions or be combined with strict-v2 outputs.

The top-level status in `bindingdb-provenance-audit-v2.json` is `FAIL` by design: the input candidate table failed. The nested `strict_dataset.status` is `PASS` and identifies the corrected canonical outputs. This distinction prevents a failed audit from being cosmetically rewritten as a pass.

Large feature arrays, raw databases, source mirrors, and the SIF are ignored; tracked manifests preserve their hashes.
