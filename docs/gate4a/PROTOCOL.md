# Gate-4A protocol: free-conformer interaction identifiability

## Objective

Determine whether independently generated free-ligand 3D information provides
reproducible target-dependent signal beyond frozen ligand 2D, additive receptor
effects, and a capacity-matched ligand-2D-by-pocket interaction control.

Gate-4A does not test docking, pose prediction, contact prediction, or a full
cross-attention architecture.

## Predictive estimands

Let `x` be frozen gMolAI 2D features, `z` free-conformer features, and `p` a
frozen receptor representation. Lower held-out risk is better.

| Model | Terms |
|---|---|
| M2D | f(x) |
| M3D | f(x) + g(z) |
| MA | f(x) + g(z) + h(p) |
| MI2 | MA + I2(x,p) |
| MI23 | MA + I2(x,p) + I3(z,p) |

Primary contrasts are:

- `delta_3d_ligand = risk(M2D) - risk(M3D)`;
- `delta_pocket_additive = risk(M3D) - risk(MA)`;
- `delta_3d_x_pocket = risk(MI2) - risk(MI23)`.

The full additive 2-by-2 factorial is mandatory to reveal order dependence. The
primary interaction comparison assigns the same total low-rank budget to MI2
and MI23. A nested extra-capacity comparison is sensitivity-only.

`delta_3d_x_pocket` is statistical non-additivity. It does not by itself prove
an atomistic binding mechanism.

## Phase ordering

1. **4A-0 source qualification:** recover source-cell semantics, audit censoring,
   map constructs and compounds, quantify matrix density and 3D informativeness.
2. **4A-1 deterministic representations:** fixed conformer generator, invariant
   shape/pharmacophore-distance descriptors, geometry-destruction controls.
3. **4A-2 frozen pretrained representation:** one pretraining-only Uni-Mol2
   molecular checkpoint; no docking/PDBbind-fine-tuned weights.
4. **4A-3 low-capacity effects:** additive and low-rank bilinear/gated models only.
5. **4A-4 locked replication:** one authorized evaluation on a distinct panel.

No trainable E(3) ligand encoder is allowed in Gate-4A.

## Outcome handling

The primary absolute endpoint is censor-aware negative log likelihood. Numeric
Kd cells additionally receive exact-value MAE/RMSE, explicitly conditional on
quantification. Supportive endpoints are classification at 100 nM and 1
micromolar, within-target ligand ranking, within-ligand target/selectivity
ranking, and two-way residual diagnostics.

Censored cells contribute inequality constraints and are never ranked among
themselves. A gain confined to exact-value regression is not progression-grade.

## Evaluation

- Double-cold ligand-scaffold and target structural/family components.
- Multiway paired bootstrap over scaffold and receptor components.
- Canonical AlphaFold structures for the admitted standardized wild-type reference
  domain are the primary view after residue-level coordinate qualification.
- Strict global-zero-nonpolymer X-ray structures are the non-holo replication view;
  candidate availability alone does not qualify a coordinate view.
- Same folds, preprocessing, tuning budget, and seeds for every paired contrast.
- No test-label access until dataset, representations, model ranks, metrics,
  equivalence regions, and exclusions are frozen.

Mandatory controls include ligand 2D only, ligand free-3D only, pocket only,
additive fusion, 2D-by-pocket interaction, KNN, sequence baseline, one-conformer
versus ensemble, topology-only distances, destroyed geometry, permuted energies,
pocket permutation, and non-holo receptor replication.

## Progression and termination

Practical-equivalence regions must be tied to assay repeatability and frozen
before fitting; no single arbitrary pKd threshold controls the decision.

Progress to a separate Gate-4B trainable E(3) study only if:

1. `delta_3d_ligand` is directionally consistent across deterministic and frozen
   pretrained representations and survives geometry-specific ablations;
2. `delta_3d_x_pocket` beats a capacity-matched 2D-by-pocket control;
3. the interaction improves both within-target and within-ligand ranking;
4. gains survive double-cold splits, censor-aware analysis, pocket permutation,
   apo/predicted views, and multiway uncertainty;
5. the result replicates on the locked external panel.

Terminate the free-conformer pivot if `delta_3d_ligand` is absent. If ligand 3D
helps but `delta_3d_x_pocket` does not, retain the result only as a ligand
representation finding and terminate the no-correspondence complementarity
hypothesis. Pocket-additive gain alone is target calibration, not interaction.

If the available data lack assay semantics, crossed coverage, conformational
diversity, or statistical power, stop as non-identifiable rather than increasing
model capacity or shopping for favourable datasets.

## Provenance-closure result

The closure phase does not authorize fitting. Davis still has 69 provisionally
admitted parents and 338 standardized wild-type reference-domain receptor rows,
while exact KINOMEscan constructs remain unavailable. The 69-row secondary chemical
review packet is technically complete but lacks a named independent signature.

Fixed KLIFS columns 1–85 define the site without a ligand. Canonical AlphaFold v6
coordinates pass exact mapping and coordinate checks for 288/338 estimands; 50 fail
closed. The corrected contract explicitly represents alignment-gap columns rather
than requiring 85 physical residues. A strict global-zero-nonpolymer apo view passes
for 34 targets, and the audited binding-site-unoccupied tier passes for 39. The
predeclared structural-edge union is now frozen, yielding only 10 components and a
dominant component of 323/338 estimands. The threshold will not be loosened after
observing this collapse.

OKL 2026 and KIRHub 2026 outcome-blind eligibility ledgers are frozen without using
outcomes for tuning. Both release zero strict pairs: OKL lacks double-cold receptor
support under the frozen components and KIRHub lacks machine-readable compound
structures plus completed construct validation. Public data still do not contain a
representative same-estimand raw paired-Kd repeat panel, so numeric metric-specific
equivalence regions remain blocked.

`delta_3d_ligand` remains blocked until an external chemical reviewer signs the
frozen packet. `delta_pocket_additive` and `delta_3d_x_pocket` remain blocked. The
dominant structural component is now a design-level identifiability problem, not
merely unfinished provenance.
