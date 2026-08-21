# Gate-4A initiation report

## Decision

**Continue source qualification only; do not train.** The free-conformer
pipeline is technically viable, but the scientific dataset is not yet admissible.
Compound identities, assay constructs, receptor views, the locked external panel,
and practical-equivalence regions remain unresolved.

This report contains no affinity-model result and makes no claim that ligand 3D
or ligand-3D-by-pocket interaction is predictive.

## Completed boundary

- Removed obsolete Gate-0--3 working data, generated features, environments,
  scripts, reports, and stale graph products after preserving the complete
  tracked terminal state at Git tag `gate3-terminal-2026-08-21`.
- Froze the three estimands (`delta_3d_ligand`, `delta_pocket_additive`, and
  `delta_3d_x_pocket`), capacity-matched interaction controls, censor-aware
  endpoint, double-cold evaluation, and termination rules.
- Acquired and hash-locked the official Davis and Karaman supplementary sources.
- Implemented strict manifest verification, affinity-cell semantics, Davis
  source/metadata audits, and deterministic invariant conformer descriptors.
- Acquired a frozen, candidate-only PubChem name-resolution snapshot and created
  a 72-row manual identity-review ledger.

## Evidence

| Check | Result | Status |
|---|---:|---|
| Davis matrix dimensions | 442 targets x 72 compounds = 31,824 tested pairs | pass |
| Exact numeric Davis Kd | 9,424 (29.6%); range 0.016--9,900 nM | pass |
| Right-censored Davis Kd | 22,400 (70.4%) | pass |
| Exact numeric 10,000 nM values | 0 | pass; common imputation is invalid |
| Invalid/missing source cells | 0 / 0 under the publisher semantics | pass |
| Cross-table target order | exact match | pass |
| Cross-table compound names | 1 mismatch: `INCB018424` vs `INCB18424` | adjudicate |
| Assay target labels / unique RefSeq accessions | 442 / 384 | descriptive |
| Mutant assay rows | 54 | requires state-aware mapping |
| Exact construct sequence/boundaries | absent | **fail/blocker** |
| PubChem name-resolution candidates | 71 unique CIDs/InChIKeys | candidate-only |
| Unspecified source derivative | 1 (`BIBF-1120 (derivative)`) | quarantined |
| Manually accepted compound mappings | 0 / 72 | **fail/blocker** |
| Full conformer preflight | 69 / 71 candidates succeeded | technical pass |
| Rejected multicomponent candidates | 2 | correct fail-closed behavior |
| Conformer ensembles with >=2 retained members | 69 / 69 successes | technical pass |
| Variable invariant geometry features | 185 / 201 | technical pass only |
| Independent full conformer replay | byte-identical SHA-256 | pass |
| Locked external panel | not selected | **fail/blocker** |
| Gate-4A model fitting | not started | correct |

The Karaman primary supplement independently confirms the same censoring rule
for its declared 38-by-317 matrix. Its PDF matrix has not undergone a lossless
cell reconstruction, so no quantitative Karaman census is reported.

## Leakage and information-boundary diagnostics

- Gate-3 data are spent and absent from active selection inputs.
- Conformer production accepts only SMILES plus a frozen config. Tests cover API
  boundaries, deterministic repeatability, and rigid translation/rotation
  invariance.
- Integration tests verify that the conformer canonical heavy-atom order matches
  the hash-pinned gMolAI canonical order and `node_z` row count, including a
  stereochemical example.
- The preflight accessed no affinity values, protein/pocket inputs, crystal ligand
  coordinates, docked poses, or complex supervision.
- Candidate PubChem structures are explicitly marked forbidden for modeling.
- No desalt, tautomer, protomer, stereochemistry, or largest-fragment decision is
  made implicitly.

## Risks and required decisions

1. **Construct non-identifiability.** RefSeq accessions and assay labels do not
   reveal recombinant domain boundaries. Full-length predictions would not be
   exact-construct receptor views. The committee should either obtain construct
   records from the assay provider/authors or predeclare a weaker reference-domain
   receptor estimand. The latter materially changes the protocol.
2. **Compound identity.** Each candidate needs comparison to the publication's
   drawn structure and a second independent source. Salts, parent forms,
   stereochemistry, and the unnamed derivative require explicit dispositions.
3. **Confirmation governance.** A non-overlapping locked panel and custodian must
   be chosen before development outcomes are inspected. Karaman is too related to
   serve this role.
4. **Representational confounding.** Pharmacophore presence/counts are separated
   from geometry groups. Geometry-destruction and topology-only controls remain
   mandatory; free-3D gains alone will not establish pocket complementarity.

## Next bounded checkpoint

Proceed to dataset admission only if all of the following are achieved:

- high-confidence compound structures and stereochemistry for a preregistered
  inclusion set, with unresolved compounds quarantined;
- a frozen salt/protomer/tautomer policy and successful conformer generation;
- documented assay construct or defensible reference-domain mappings, including
  mutation and phosphorylation state;
- query-independent apo/predicted receptor and pocket provenance;
- a signed development/test separation and locked external panel;
- scaffold/target-component effective sample-size and power audits supporting the
  three preregistered contrasts.

If exact or defensible receptor-site mapping cannot be established, terminate the
`delta_3d_x_pocket` experiment on Davis rather than substituting undocumented
full-length receptors or escalating architecture complexity. A ligand-only
`delta_3d_ligand` experiment could still be evaluated separately, but it would not
test iScore3.0's proposed complementarity hypothesis.
