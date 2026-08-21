# Gate-4A dataset-admission report

> Historical dataset-admission snapshot. The current decision and completed
> provenance audit are in
> [GATE4A_PROVENANCE_CLOSURE_REPORT.md](GATE4A_PROVENANCE_CLOSURE_REPORT.md).

**Decision: `Δ3D-ligand` is conditionally admissible for a later ligand-only
experiment; `Δ3D×pocket` remains `BLOCKED`. No predictive model was trained.**

This is the scientifically conservative outcome. Davis can support a standardized
wild-type reference-domain question, but it cannot establish the exact receptor
reagents used in the original assay. The ligand-independent site definition is
frozen, while coordinate-level pocket mapping, structural leakage edges, external
eligible-pair ledgers, and a numeric noise-derived equivalence region are not yet
qualified. Training an interaction model now would turn unresolved provenance into
apparent statistical signal.

## Admission decisions

| Item | Decision | Evidence |
|---|---|---|
| All 72 Davis compound identities | **PASS, committee QA required** | 69 parents admitted; Ki-20227 and SB-203580 quarantined for unresolved stereo; unnamed BIBF-1120 derivative quarantined |
| Davis assay censoring | **PASS** | 10,000 nM blanks retained as right-censored, never exact; no exact source value equals 10,000 nM |
| Exact KINOMEscan constructs | **BLOCKED** | Source does not report recombinant sequences/boundaries |
| Standardized WT reference-domain estimand | **PASS for 338/442 rows** | Reviewed canonical UniProt core domain, one unique complete 85-position KLIFS pocket, canonical AlphaFold availability; 104 state/mutant/atypical/ambiguous rows excluded |
| Ligand-independent pocket definition | **PASS** | Fixed KLIFS positions 1–85; no query ligand, affinity, contact, or holo-derived site selection |
| Predicted pocket coordinates | **BLOCKED** | Metadata coverage is 338/338, but unique residue-level KLIFS→UniProt→CIF mappings and coordinate hashes are not materialized |
| Strict apo candidates | **PASS as search; coordinates BLOCKED** | 160/338 accessions have 558 unique X-ray entries with zero non-polymer entities; WT pocket sequence/completeness audit remains |
| Davis admitted design | **PASS descriptively** | 69 ligands × 338 receptors = 23,322 tested pairs; 6,581 exact and 16,741 right-censored |
| Protein structural leakage edges | **BLOCKED** | KLIFS family/sequence components are frozen; coordinate TM-score/Foldseek edges await admitted pockets |
| External confirmation source | **PASS conditionally; release BLOCKED** | OKL 2026 for affinity-like confirmation plus KIRHub 2026 for orthogonal ranking/classification; both need eligible-pair adjudication and custodian lock |
| Practical-equivalence regions | **BLOCKED numerically** | Public repeat evidence is aggregate/correlational, not representative raw paired Kd repeats |
| Matched interaction budget | **PASS** | All modalities projected to 32 frozen PCA dimensions; MI2 rank 8 and MI23 rank 4+4 each have exactly 512 trainable bilinear-factor parameters |

## What the admitted data can and cannot identify

The admitted matrix contains a median 18.5 exact ligands per receptor; 316 targets
have at least eight exact ligands, 221 have at least fifteen, and 314 span at least
eight exact ligand leakage components. Exact pKd spans 5.004–10.620 (median 6.187).
This is adequate to ask whether deterministic free-conformer features improve a
ligand-only model, subject to final committee review of the chemical ledger.

It is not yet adequate to ask whether ligand 3D interacts with pocket 3D. The 23,322
cells are not independent. Outcome-independent grouping yields 66 ligand components
and 84 receptor components before coordinate structure edges; the smaller-axis ESS
upper bound is therefore 66 and can only decrease. A no-model power sensitivity
calculation gives a detectable paired mean contrast of about 0.069 when the
component-level contrast SD is 0.20, but this is neither a progression threshold nor
a substitute for the preregistered multiway bootstrap.

Free ligand geometry is non-degenerate: all 69 admitted parents generated
deterministically, all retained 2–31 conformers, 68/69 had mean pairwise heavy-atom
RMSD above 0.5 Å, and 185/201 geometry features varied. This passes only an input
informativeness test; it is not affinity evidence.

## Independent-panel audit

The 2026 Optimal Kinase Library is the best affinity-like candidate found: 192
compounds × 468 targets and 89,856 inferred pair estimates, with raw multi-dose
percent-control data in the public workbook. It is not a clean exact-Kd test.
Sixty-one compounds exactly overlap admitted Davis parents; only 74,070 pairs have
four raw doses and 15,786 have two; 61,565 Bayesian point estimates lie above the
highest tested dose and 2,117 below the lowest. The workbook's “KSKd (µM)” column is
numerically nM, an internal unit-label inconsistency. OKL must therefore use raw
dose/interval/censor-aware endpoints after scaffold-disjoint eligibility filtering,
not pretend that all inferred values are exact Kd. See the [OKL primary report](https://pmc.ncbi.nlm.nih.gov/articles/PMC13015435/).

KIRHub is the preferred orthogonal confirmation: 92 clinical inhibitors × 409
wild-type kinases in a duplicate HotSpot functional assay, with 36,593/37,628 public
residual-activity cells and construct metadata for 391 kinases. It also supplies a
14-compound × 369-kinase ten-dose subset. It is useful for ranking/classification,
not absolute pKd. Compound structures and raw replicate pairs are absent, 18 panel
kinases lack a construct-metadata row, and Reaction Biology reuse terms require
committee review. See the [KIRHub primary report](https://www.nature.com/articles/s41587-026-03090-8).

Karaman is not independent of Davis. IDG-DREAM is sparse and hit-preselected.
Anastassiadis/PKIS/PKIS2 are useful orthogonal screens but cannot serve as sole Kd
confirmation. The final policy is complementary one-shot ledgers: OKL for
affinity-like replication and KIRHub for orthogonal ranking/classification, with no
external label used for tuning.

## Noise and equivalence

The available vendor example reports less than four-fold Kd variation for one
compound–kinase pair over three lots; this is only a coarse 0.602-log10 span. KIRHub
reports duplicate-screen R²=0.99 but publishes averages, and OKL repeats are
QC-triggered rather than representative. None identifies the paired Davis Kd error
distribution. A numeric margin would therefore be invented.

The frozen remedy is metric-specific: obtain raw independent repeats, swap replicate
A/B labels within exact ligand–construct pairs, run the final multiway component
bootstrap, and set each practical-equivalence half-width to the 95th percentile of
the absolute metric contrast caused solely by replicate choice. A single shared pKd
threshold is forbidden.

## Recommendation

Do not train receptor interaction models. The next permissible bounded work is
provenance, not architecture: committee chemical QA; coordinate-qualified canonical
predicted pockets; strict apo chain audit; coordinate structural leakage closure;
OKL/KIRHub chemical and receptor eligibility ledgers held by a custodian; and raw
repeat acquisition or a prospective repeat panel. `Δ3D-ligand` may then proceed
independently. `Δ3D×pocket` must remain blocked until every listed condition passes.

Machine-readable evidence is under [`reports/gate4a/evidence/`](evidence/), the
admission ledgers under [`data/processed/gate4a/`](../../data/processed/gate4a/),
and leakage components under
[`data/splits/gate4a/davis-admission-components-v1.tsv`](../../data/splits/gate4a/davis-admission-components-v1.tsv).
