# Gate-2 report bundle

Gate-2 is **no-go for full iScore3.0 architecture training** under the frozen,
strict component/OOD protocol. The result is specific: the shallow ligand–pocket
interaction control did not show reproducible incremental value beyond matched
additive fusion. No cross-attention or other full iScore3.0 model was trained.

- [GATE_2_REPORT.md](GATE_2_REPORT.md) — decision, methods, results, risks, and recommendation
- [reproducibility.md](reproducibility.md) — pinned inputs, runtime, commands, and expected hashes
- [evidence/baseline-metrics-v1.json](evidence/baseline-metrics-v1.json) — all out-of-fold metrics and paired intervals
- [evidence/leakage-diagnostics-final-v1.json](evidence/leakage-diagnostics-final-v1.json) — split-overlap and residual-risk audit

Canonical inputs are Gate-2 v3 and baseline outputs are v1. Earlier v1/v2
candidate/manifests remain an immutable mapping and protocol-amendment trail.
