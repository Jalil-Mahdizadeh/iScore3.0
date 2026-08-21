# Gate-4A data and test governance

## Frozen dataset roles

- **Davis 2011:** development source for censor-aware dense kinase-panel work.
- **Karaman 2008:** related technical evidence only; not an independent final test.
- **OKL 2026:** conditional affinity-like confirmation candidate. Use raw dose data
  with interval/censor-aware endpoints, never all reported posterior estimates as
  exact Kd. Davis-overlapping compounds and leakage-connected components are
  ineligible.
- **KIRHub 2026:** conditional orthogonal functional ranking/classification
  candidate, not an absolute-pKd panel.
- **Gate-3:** permanently quarantined for architecture or threshold selection.

Confirmation labels remain unavailable to model development. Eligible-pair
ledgers must be chemically and structurally adjudicated, hash-frozen, and held by
an independent custodian before either external source can be released once.

## Source and label contract

For every matrix cell preserve publication, table, source row/column, assay,
construct annotation, relation, value, unit, censor limit, and whether the pair
was assayed. Derivative ML matrices are comparison aids, never the source of
truth.

Davis Supplementary Table 4 has 31,824 tested pairs: 9,424 numeric Kd values and
22,400 blanks. A blank means Kd above 10 micromolar or no detection in the
10-micromolar screen and is encoded as `right_censored_kd`. No numeric source
cell equals 10,000 nM. After chemical and receptor admission, the provisional
matrix is 69 ligands by 338 receptors, with 6,581 exact and 16,741 censored
cells.

## Chemical identity outcome

All 72 Davis compounds have explicit dispositions in
`data/processed/gate4a/davis-compound-adjudication-v1.tsv`. Sixty-nine parent
structures are admitted. `BIBF-1120 (derivative)` is unidentifiable; Ki-20227
has unspecified carbon stereochemistry; and SB-203580 has unspecified sulfoxide
stereochemistry. Those three are quarantined. Two hydrochloride registry records
are converted to the publication parent. The INCB018424/INCB18424 ruxolitinib
alias inconsistency is explicitly resolved, never silently normalized.

The review was label-blind and checked publication drawings, PubChem records,
ChEMBL exact-synonym evidence, and an independent NCATS record where needed.
The automated integrity pass is complete, but self-review is not independent.
Committee or other named external secondary review of the frozen 69-row packet is
required before any label-bearing fit.

## Standardized receptor estimand

The exact KINOMEscan recombinant constructs are not reported and cannot be
reconstructed. The primary estimand is therefore deliberately different:

> reviewed human canonical UniProt wild-type core kinase domain, represented at
> fixed KLIFS alignment positions 1--85 on the canonical AlphaFold structure.

Rows with mutations, phosphorylation-state qualifiers, partners, atypical or
incomplete sites, non-human proteins, or non-unique mappings are excluded.
This admits 338/442 Davis assay labels. It must never be described as an exact
assay-construct analysis.

## Ligand-independent structure boundary

The pocket is the ordered set of KLIFS alignment columns 1--85. A column may be an
explicit reference-alignment gap and then has no physical residue or coordinate.
Selection may not use a
query ligand, affinity label, ligand contact, holo template, pharmacophore, or
complex-derived exclusion volume. The canonical predicted view is primary; a
strict apo X-ray view and a separately labelled binding-site-unoccupied X-ray tier
are replication views. Holo structures are sensitivity-only and cannot define a
pocket.

Metadata coverage alone does not admit coordinates. Each predicted or apo view
must have a unique KLIFS-to-UniProt-to-coordinate residue mapping, complete required
coordinates, exact sequence agreement at every non-gap column, structure-file
hash, and chain/model provenance. The strict apo tier additionally requires zero
non-polymer entities. The binding-site-unoccupied tier may contain remote
nonpolymers but requires no non-water foreign heavy atom within 8 Å of the pocket,
exact WT pocket sequence, full C-alpha coverage, at least 90% side-chain-heavy-atom
completeness, and X-ray resolution at most 3 Å.

Leakage components are the union of ligand scaffold/similarity edges and receptor
family, pocket-sequence, and pocket-structure edges. Coordinate edges are frozen at
maximum length-normalized pocket US-align TM-score at least 0.75. The resulting
union has 10 receptor components and a 323-target dominant component; it may not be
repartitioned after label or model inspection.

## Noise, equivalence, and access control

The public repeat evidence does not identify a representative paired Kd error
distribution, so no numeric practical-equivalence margin is admitted. The future
margin is metric-specific: within exact ligand--construct independent repeats,
swap replicate A/B assignments, rerun the final multiway component bootstrap,
and set the half-width to the 95th percentile of the absolute metric contrast
caused by replicate choice. A universal pKd threshold is forbidden.

Development and confirmation manifests remain separate. The outcome-blind OKL and
KIRHub eligibility ledgers are frozen and currently release zero strict pairs.
Confirmation access is
one-shot after signed hashes for code, configuration, chemical/receptor ledgers,
features, splits, checkpoints, and analysis plan are recorded. The prior named-
reviewer requirement is superseded by the final owner-authorized 69-ligand ledger;
no additional external/manual chemical sign-off is required. The isolated
`Delta3D-ligand` analysis is complete with a no-go result. All receptor/additive and
`Delta3D x pocket` fitting remains blocked because the structural components are
non-informative for a Davis double-cold test, strict external pairs are unavailable,
and replicate-derived equivalence margins are unresolved.
