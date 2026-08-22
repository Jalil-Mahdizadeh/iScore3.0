# iScore3.0

iScore3.0 is currently between hypotheses. Gate-4A is terminated with a
**NO-GO** decision: independently generated free-ligand conformers did not add
reproducible predictive information beyond frozen gMolAI 2D under the frozen,
ligand-component-OOD experiment. Receptor-additive and ligand-pocket interaction
models were not trained because receptor structural novelty was not identifiable
under the preregistered Davis leakage graph.

## Historical recovery

- [`reports/history/GATE4A_CLOSURE.md`](reports/history/GATE4A_CLOSURE.md) records
  the terminal decision and recovery instructions.
- [`reports/history/GATE3_CLOSURE.md`](reports/history/GATE3_CLOSURE.md) records
  the preceding ligand-2D plus pocket-3D interaction no-go.

The complete Gate-4A repository—including protocols, reports, evidence,
manifests, ledgers, source code, and tests—is frozen at Git tag
`gate4a-terminal-2026-08-22`. Gate-4A-specific artifacts have been removed from
the active tree, following the Gate-3 archive pattern. Reusable provenance,
ligand-encoding, pocket, apo-view, and structure-view utilities remain.

No successor hypothesis is active or authorized for implementation. Historical
gate data are spent for architecture selection and must not be silently reused
for tuning, stopping, or progression decisions in a future project.

Large source data, container images, external repositories, generated features,
and checkpoints are ignored by Git. Their historical provenance is recoverable
from the terminal tags; local ignored files are not treated as an active dataset.

No model in this repository is suitable for prospective, clinical, or regulatory
use.
