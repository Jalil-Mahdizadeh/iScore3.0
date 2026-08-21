# Gate-4A provenance-closure report

> Historical report. Superseded on 2026-08-21 for ligand QA and the isolated
> experiment: the project owner finalized all 69 admitted identities without
> additional external/manual sign-off, and `delta_3d_ligand` subsequently produced
> a no-go result. Receptor/additive and interaction blocks below remain in force.

**Decision: `delta_3d_ligand` is `BLOCKED_PENDING_INDEPENDENT_QA`; both
`delta_pocket_additive` and `delta_3d_x_pocket` are `BLOCKED`. No predictive model
was trained.**

The coordinate work produced a valid 288-target receptor subset and two defensible
non-holo tiers, but it also exposed a more consequential design failure: the frozen
structural-leakage rule connects 323/338 Davis estimands. Relaxing that rule after
seeing the component collapse would be threshold shopping. Davis therefore cannot
currently provide a credible double-cold interaction test under this protocol.

## Admission decisions

| Requirement | Decision | Audited result |
|---|---|---|
| Technical secondary chemistry checks | **PASS** | A second implementation round-trips all 69 provisional parents through RDKit/InChIKey checks and confirms explicit identity, stereochemistry, salt and protomer fields |
| Independent chemical sign-off/final ledger | **BLOCKED** | A 69-row packet is hash-frozen, but the original adjudicator cannot certify self-review as independent; named reviewer/date and four row-level decisions remain pending |
| Canonical predicted receptor view | **PASS for 288; BLOCKED for 50** | 211 direct exact mappings plus 77 reciprocal exact homology transfers; all 288 have exact residue identity and C-alpha coordinates at every non-gap KLIFS column |
| Strict global-zero-nonpolymer apo tier | **PASS as subset** | 34/288 targets |
| Binding-site-unoccupied apo tier | **PASS as subset** | 39/288 targets, including five entries with only remote nonpolymers; 42 searches hit the frozen 25-coordinate-attempt cap |
| Structural-similarity leakage closure | **PASS technically; FAIL identifiability** | 41,328 pair comparisons and 40,852 structure edges; union yields 10 components, largest 323/338 |
| OKL confirmation ledger | **FROZEN; RELEASE BLOCKED** | 117/192 non-Davis compounds await manual identity QA; no strict double-cold pairs survive current receptor provenance/components |
| KIRHub confirmation ledger | **FROZEN; RELEASE BLOCKED** | No machine-readable structures for 92 compounds; construct validation incomplete; no outcome table was read |
| Metric-specific practical equivalence | **BLOCKED** | Public cross-study/aggregate evidence is informative context but not a representative same-estimand repeat distribution |

## Chemical identity and experiment release

The automated second implementation found no technical inconsistency in the 69
admitted-parent rows. The immutable source ledger SHA-256 is
`e924f7f3fc95c25685857cbcc1bbdb8674cef40dce74ddf4ab82710c134c76f2`; the
row-complete reviewer packet SHA-256 is
`77bef1365eddf215f65d76e887ed48d2cdbaae6898494cba0cb94b5e64ae63b9`.
This is not independent scientific sign-off. Release requires all four secondary
decisions to be `PASS` for every row plus a named reviewer, date, and signed packet
hash in a new immutable manifest. Consequently, the optional isolated ligand-3D
experiment was not run.

## Receptor and pocket closure

The pocket is frozen as 85 ordered KLIFS alignment columns. Fifty-seven targets
contain at least one alignment-gap column; those gaps are explicit, not fabricated
residues. Among the 288 exact views, physical pocket size is 85 residues for 231,
84 for 56, and 82 for one. All coordinates come from canonical AlphaFold v6 mmCIF
files and are individually hashed. The lowest target-level mean pocket pLDDT is
78.72. Fifty estimands fail closed: 41 lack usable KLIFS API mapping evidence and
nine lack an exact transfer. Holo KLIFS coordinates were used only to transfer
alignment residue numbers, never as model input.

The strict apo tier retains the original whole-entry zero-nonpolymer rule. The new
binding-site-unoccupied tier permits remote additives/ions but requires exact WT
pocket sequence, full C-alpha coverage, at least 90% side-chain-heavy-atom
completeness, X-ray resolution no worse than 3 Å, and no non-water foreign heavy
atom within 8 Å. It adds exactly five views: FGFR4/4TYG, TGFBR2/5E8V, ERBB4/3BCE,
MAPK14/3OEF, and MKNK2/2AC3; their nearest foreign atoms are 8.11–19.10 Å away.
Candidate selection did not use ligand identity or affinity.

The source contracts follow the official [KLIFS](https://klifs.net/) mapping and
[AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) coordinate
records. Exact downloaded payloads and hashes are recorded under
`data/manifests/gate4a/`; rebuildable raw coordinates remain outside Git.

## Leakage and confirmation consequence

The predeclared receptor edge is the union of same KLIFS family, aligned KLIFS85
identity at least 0.70, and maximum length-normalized pocket US-align TM-score at
least 0.75. Sequential alignment is intentional because coordinates are serialized
in the same homologous KLIFS-column order. On the 288 admitted coordinate views,
40,852/41,328 pairs cross the structural threshold. Transitive closure with the
family/sequence edges leaves 10 components and a dominant 323-target component.

This is not a computational failure. It says that the selected conserved kinase
site and threshold encode nearly the whole kinome panel as one leakage group. The
current Davis matrix therefore has inadequate target-side component ESS for the
preregistered double-cold interaction contrast. A future design would need genuinely
independent target families/sites or a prospectively justified leakage estimand—not
post-hoc threshold relaxation.

OKL and KIRHub ledgers were built without reading cell-level outcomes and are frozen
without release. OKL is described in the [primary 2026 report](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015435/);
KIRHub in its [primary publication](https://doi.org/10.1038/s41587-026-03090-8).
Neither currently yields a strict confirmation panel under the frozen identity,
construct, and double-cold requirements.

## Noise evidence

The strongest newly acquired public evidence is the IDG-DREAM release comparing
Fabian and Davis values for 416 pairs. It gives contextual MAE 0.326 pKd, RMSE 0.517
pKd, Spearman 0.840, and mean Fabian-minus-Davis bias -0.149 pKd. However, 178 pairs
have at least one pKd=5 lower-bound/nonbinder value, constructs are not resolved to
the exact Gate-4A estimands, and these are cross-study averages rather than raw
independent repeats of the same assay reagent. The primary publication and exact
released analyses are preserved via [Nature Communications](https://doi.org/10.1038/s41467-021-23165-1)
and [Zenodo 4648011](https://doi.org/10.5281/zenodo.4648011).

No universal threshold is inferred from that RMSE, a vendor's less-than-four-fold
example, or published correlations. Numeric margins for censor-aware NLL, MAE,
RMSE, classification and ranking remain blocked until raw paired repeats support
the preregistered metric-specific replicate-swap procedure.

## Recommendation

Do not train receptor interaction models on this Davis design. The bounded next
actions are external signature of the ligand packet, resolution or explicit
exclusion of the 50 failed receptor mappings, and acquisition of a genuinely
independent crossed confirmation panel plus same-estimand repeat evidence. The
isolated `delta_3d_ligand` experiment becomes permissible only after chemical
sign-off; it cannot rescue receptor-interaction identifiability. If a future panel
does not restore adequate receptor-component ESS under prospectively frozen
leakage rules, terminate Gate-4A rather than increase model capacity.

All decisions and hashes are machine-readable under
[`evidence/`](evidence/); coordinate, apo and confirmation ledgers are under
[`data/processed/gate4a/`](../../data/processed/gate4a/); final receptor components
are under [`data/splits/gate4a/`](../../data/splits/gate4a/).
