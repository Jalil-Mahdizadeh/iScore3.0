# Design only: prospective kinase-selectivity interaction gate

Status: `DESIGN_ONLY / NO TRAINING AUTHORIZED` (2026-08-21).

## Choice and rationale

The preferred next interaction experiment is a newly preregistered **kinase-specific selectivity** estimand, not a post-hoc weakening of the Davis broad receptor-structural-novelty estimand. A genuinely heterogeneous multi-fold crossed affinity matrix would be cleaner for broad transfer, but no already-audited public panel currently satisfies dense crossing, comparable quantitative affinity semantics, ligand identity, construct provenance and non-holo receptor requirements. Claiming one exists would substitute dataset convenience for estimand clarity.

The proposed question is narrower: can ligand free-3D and a ligand-independent kinase pocket predict selectivity across prospectively held-out kinase groups/families and ligand scaffolds? This does not claim transfer to unrelated protein folds.

## Fresh-data governance

- Davis and all Gate-3 outcomes are permanently spent for architecture, threshold and control selection. They may be cited as prior evidence but never enter fitting, early stopping or hyperparameter choice.
- A new outcome matrix remains sealed until ligand identity, assay semantics, construct mapping, receptor views, folds, models, budgets and decision rules are frozen. OKL and KIRHub remain label-locked candidates, not automatically admitted datasets.
- Development, confirmation and final test panels must be publication/vendor-batch separated. Molecules connected at exact parent, Bemis-Murcko or Morgan2 Tanimoto >= 0.60 are one leakage component across all panels.
- Target units are prospectively defined by wild-type canonical kinase-domain estimands. Mutants, phosphorylated-state duplicates, fusion constructs and unresolved constructs are excluded from primary analysis or isolated as separate sensitivity strata.

## Dataset admission

Admission requires a crossed quantitative binding matrix (exact Kd preferred) with explicit censoring, at least 100 admitted ligand parents, at least 120 kinase-domain targets, all major kinase groups represented, at least 8 independent kinase family holdout units, and at least 20 ligand scaffold components per held-out target unit. Each admitted target needs a hashed predicted wild-type view and either a strict-apo or audited binding-site-unoccupied view for replication. No minimum is waived after outcome access.

Panels based only on single-concentration percent inhibition, mixed IC50/Ki/Kd without strata, imputed inactives, unresolved assay direction, or target-dependent ligand sampling fail admission. Raw repeats or a separately acquired repeat subset must support metric-specific practical-equivalence margins; otherwise only evidence against zero gain can be reported, not practical equivalence.

## Frozen estimands and splits

Primary targets are held out by recognized kinase group/family units defined before coordinate similarity is computed. Within development groups, leave-one-family-out folds are used; final testing seals entire kinase groups when sample size permits. Ligand scaffold/similarity components are held out simultaneously, producing double-cold cells. Sequence identity, KLIFS85 identity and pocket-structure similarity are reported as continuous diagnostics and forbidden from being threshold-tuned against outcomes.

Primary response is censor-aware affinity likelihood. Selectivity endpoints are within-ligand target ranking and pairwise target preference; within-target ligand ranking and exact-cell RMSE/MAE are co-primary supportive endpoints. Absolute pKd is secondary because batch and target offsets can mimic receptor information.

## Matched models and falsification controls

The first gate remains low capacity. Training a full cross-attention architecture is forbidden.

1. ligand-only gMolAI and gMolAI+free-3D main effects;
2. pocket-only and sequence-only controls;
3. additive ligand+pocket controls;
4. capacity-matched bilinear/FiLM interactions comparing 2D×pocket with 2D×pocket + 3D×pocket under the already frozen projection widths and factor budgets;
5. ligand-coordinate destruction, topology-fake 3D, energy permutation and conformer-seed repeats;
6. pocket permutation within kinase group, receptor-coordinate destruction, and target-label permutation;
7. nearest-neighbour baselines on ligand, sequence and pocket structure.

All pocket definitions use conserved kinase-alignment positions on predicted/apo structures and are independent of query ligands. Holo coordinates, contact labels, docking poses/scores, complex-derived supervision and ligand↔pocket correspondence remain forbidden.

## Progression rule

Progression requires an interaction-specific gain over both the matched additive model and the capacity-matched 2D×pocket model, with component-bootstrap uncertainty excluding zero on the primary censor-aware endpoint; concordant within-ligand selectivity ranking and double-cold scaffold/family performance; failure of pocket and ligand geometry permutation controls to reproduce the gain; positive replication in at least one non-holo receptor view and a sealed independent panel; and practical magnitude relative to repeat-derived metric margins.

Failure to admit a sufficiently crossed panel, collapse of the prospectively chosen family units, or failure of the interaction-specific controls is a stop decision. Model complexity must not be increased to compensate.
