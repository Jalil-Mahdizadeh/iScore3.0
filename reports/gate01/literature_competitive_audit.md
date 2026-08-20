# Gate-0/1 Literature and Competitive Audit

**Cut-off:** 20 August 2026<br>
**Purpose:** update the original blueprint with missing direct competitors and classify effective ligand-coordinate use from papers and source code.

## Finding

The proposed modality is not novel. PSG-BAR, BANANA, CASTER-DTA, 3DProtDTA, AttentionMGT-DTA, PLMCA, HoloProt, and Graph_RG already establish protein-structure plus molecular-graph/SMILES affinity or activity prediction without a query docking pose. Cross-attention and structural pocket encoding are also precedented.

The viable research question is narrower: can an end-to-end ligand-3D-free, provenance-first system demonstrate reproducible incremental value from receptor-pocket geometry beyond ligand, sequence, metadata, and nearest-neighbour controls under similarity-component, conformation, temporal, and external tests?

## Gate-0/1 additions

| Work | Year | Effective modality | Class | Public-source conclusion | Effect on iScore3.0 |
|---|---:|---|---|---|---|
| PSG-BAR | 2022 | full protein 3D graph + 2D ligand graph | strict direct | pinned; k discrepancy, raw xyz/non-invariant GAT, no root licence | explicit direct prior art; clean-room comparator candidate |
| PLMCA | 2026 | PLMs + protein 3D/physicochemistry + ligand graph + assay context | provisional strict direct | no official source found | removes broad multimodal/cross-attention novelty; implementation boundary unverified |
| AttentionMGT-DTA | 2024 | pocket 3D/ESM + 2D ligand graph | strict direct | pinned; test-set-driven scheduling/selection and pocket-loop issue | direct comparator, but published protocol needs leakage-safe reproduction |
| BlendNet | 2025 | pocket/protein sequence + 2D ligand graph at inference | privileged training | pinned; PLIP complex labels train teacher | not strict training; useful upper comparator |
| AttentionSiteDTI | 2022 | structural binding-site graphs + 2D ligand | strict binary DTI | pinned; CC BY 4.0, old monolithic pipeline | modality precedent, not quantitative-affinity evidence |
| BindingSite-AugmentedDTA | 2023 | structural-site augmentation of DTA models | adjacent direct | pinned; repository mainly tables/README | site-information precedent; not turnkey |
| PGraphDTA | 2023 | PLM + contact map + ligand graph | structure proxy | cited repository unavailable | relevant contact-map baseline; weaker than local 3D chemistry |
| CSCo-DTA | 2024 | protein contact/network graphs + ligand graph | structure proxy | pinned; external contact maps and network inputs | transductive/proxy control, not strict pocket-geometry proof |
| MMPD-DTA | 2025 | sequence + ligand graph + pocket–drug complex graph | non-strict | pinned; graph construction absent | title is adjacent but model requires pose geometry |
| AlignNet | 2026 | ESM/GearNet + MolFormer/GraphMVP + complex graph | non-strict/privileged | pinned; ligand MOL2 coordinates and cross-complex edges confirmed | “structure-agnostic” is not ligand-3D-free |
| LigUnity | 2025 | shared 3D pocket–ligand Uni-Mol space | non-strict | pinned; ligand conformers and ligand-defined sites confirmed | high-priority docking alternative, outside strict claim |

## Source-code findings that affect scientific interpretation

- **PSG-BAR:** paper/source nearest-neighbour parameter mismatch (`5` versus `3`) and raw Cartesian coordinates in a standard attention graph mean reported gains do not establish rigid-motion-safe structural learning.
- **AttentionMGT-DTA:** `train_DTA.py` evaluates the test set every epoch, advances the scheduler with test MSE, and records the best test epoch. Its headline numbers should not be reused as if they came from a locked test.
- **BlendNet:** exact new-protein/new-compound/blind partitions are useful but do not enforce analogue or pocket-similarity component separation; the teacher’s PLIP labels are privileged complex information.
- **MMPD-DTA:** released code loads precomputed complex graphs and does not reproduce their geometric construction, preventing a complete feature-lineage audit.
- **AlignNet/LigUnity:** source inspection contradicts any strict interpretation based only on their “structure-agnostic” or docking-alternative framing.

All pinned revisions and exact inspected-file hashes are in [`third_party/reviewed_sources.tsv`](../../third_party/reviewed_sources.tsv). Source mirrors are intentionally ignored.

## Leakage perspective

Warm/random results cannot establish the intended deployment claim. PSG-BAR’s reported cold-both deterioration is an instructive example. Exact target/drug splits in several competitors still permit close ligand analogues, homologous proteins, similar pockets, shared publications, and pretrained-entity exposure. Conversely, similarity to training is not automatically illegitimate for an interpolation deployment. The required practice is:

1. define the intended novelty regime;
2. prevent exact and prohibited provenance overlap;
3. form union components over ligand, sequence, pocket, construct, and publication edges;
4. publish continuous nearest-neighbour distributions; and
5. keep random splits as shortcut diagnostics rather than headline evidence.

## Defensible claim boundary

Do not claim:

- invention of docking-free affinity prediction;
- invention of protein-3D plus ligand-2D fusion;
- invention of pocket–ligand cross-attention; or
- superior generalization from random or exact-ID cold splits.

A future claim becomes defensible only if the strict feature boundary is verified through all pretraining, preprocessing, training, and inference stages and receptor geometry passes a predeclared incremental-value test on independent, conformation-stressed and external data.

The complete 82-record review is in [`blueprint/evidence/publications.tsv`](../../blueprint/evidence/publications.tsv), with synthesis in [`blueprint/02_literature_and_software_review.md`](../../blueprint/02_literature_and_software_review.md).
