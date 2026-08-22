# Gate-4A closure and recovery

Gate-4A is terminated as a failed scientific hypothesis. Under its frozen
ligand-component-OOD protocol, neither deterministic free-conformer descriptors
nor frozen Uni-Mol v1 supplied reproducible predictive information beyond frozen
gMolAI 2D. The geometry controls did not support a geometry-specific explanation:
actual Uni-Mol coordinates were worse than coordinate-destroyed controls, and
single-versus-ensemble and energy-permutation contrasts were inconclusive.

The primary pooled NLL gains over gMolAI were `-0.00257` for deterministic 3D
(95% component-bootstrap interval `[-0.01780, 0.01359]`) and `-0.01356` for
Uni-Mol (95% interval `[-0.03111, 0.00410]`), where a positive value would favour
3D. Only one of five deterministic seeds and none of five Uni-Mol seeds improved
over gMolAI. The preregistered joint progression rule therefore failed.

The Davis receptor-additive and ligand-pocket interaction branches were never
trained. They remained blocked because the frozen structural-leakage union placed
323 of 338 standardized receptors in one component, making the intended broad
receptor-structural-novelty estimand non-identifiable. No leakage rule was relaxed
after observing that collapse.

The complete terminal repository—including the final report, machine-readable
evidence, frozen protocols, manifests, ledgers, ignored-artifact hashes, source
code, and tests—is recoverable at:

- Git commit: `56be8ebba3a5e8a58cdcf2931c610a34172223f2`
- Git tag: `gate4a-terminal-2026-08-22`
- Terminal report at that tag:
  `reports/gate4a/GATE4A_DELTA3D_LIGAND_REPORT.md`
- Machine result at that tag:
  `reports/gate4a/evidence/delta3d-ligand-results-v1.json`

Gate-4A protocols, Davis outcomes, control results, and receptor-component
diagnostics are spent for architecture selection. They must not be used as an
apparently fresh development or confirmation set for a successor hypothesis.
Any retrospective reuse must be explicitly labelled, outcome-unblinded, and
incapable of determining tuning, stopping, or progression.

Gate-4A-specific tracked artifacts were removed from the active tree after the
terminal tag was verified. Ignored local raw structures and rebuildable features
were not erased; they remain inactive local recovery material rather than an
admitted dataset.
