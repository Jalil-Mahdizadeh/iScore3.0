# Gate-0/1 Strict S0/S1 Pilot Data Card

## Identity and intended use

**Dataset:** `rcsb-bindingdb-kd-gate01-strict-v2`<br>
**Created:** 20 August 2026<br>
**Canonical manifest:** [`data/manifests/gate01-pilot-strict-v2.json`](../../data/manifests/gate01-pilot-strict-v2.json)

This is a small diagnostic dataset for validating provenance, ligand/protein feature boundaries, split construction, and shallow baselines. It is not a training corpus for the full iScore3.0 architecture, not a blind benchmark, and not evidence of clinical or prospective utility.

## Sources

- RCSB PDB Search/Data APIs supplied the candidate Kd annotations, construct sequences, entity/chain mappings, ligand components, experimental metadata, citations, dates, and mmCIF structures.
- BindingDB release 202608 supplied the record-level affinity relation/value, PDB/ligand/UniProt cross-references, measurement publication, and curation provenance.
- The ignored official BindingDB ZIP is 592,986,963 bytes with SHA-256 `8dc30541d668e5e403ba1fe9a682b76c42a50e1711aac44329c62ba24ccd57db`. Its 8,978,755,790-byte TSV contained 3,234,499 data rows.

Every tracked table and upstream acquisition/selection/structure/site manifest has an immutable hash in the canonical manifest. Raw coordinates and the BindingDB archive are intentionally not committed.

## Selection and correction history

The initial label-blind structure-quality selection yielded 80 supervised candidate complexes and six historical site references across six exact construct groups. gMolAI was audited on those 80 candidate ligands before affinity reconciliation.

BindingDB reconciliation required all of the following for a supervised row:

1. exact PDB identifier;
2. exact ligand component or matching InChIKey;
3. matching UniProt accession;
4. exact, uncensored Kd with compatible value/unit; and
5. an identified BindingDB measurement publication.

Twelve candidate values were actually `>`-qualified lower bounds even though the RCSB query representation appeared numeric. They were removed. Of the remaining 68 measurements, the WDR5 construct group retained only seven and failed the frozen minimum of eight, so the complete group and its reference were excluded. No target or ligand was added after examining baseline performance.

## Final composition

| Construct group | Protein | UniProt | Construct length | Supervised rows |
|---|---|---|---:|---:|
| `construct-001-c5158f5b655d` | nitric oxide synthase oxygenase | O34453 | 363 | 19 |
| `construct-003-fbf3f5c695a3` | bromodomain-containing protein 4 | O60885 | 127 | 15 |
| `construct-002-cb5cf2af5e4a` | Peregrin | P55201 | 116 | 10 |
| `construct-004-a02c0b37ee39` | cellular tumor antigen p53 | P04637 | 219 | 9 |
| `construct-005-6f4aa0a92e90` | mycocyclosin synthase | P9WPP7 | 396 | 8 |

There are 61 exact uncensored Kd labels, five blank-label site references, and 122 receptor-view feature rows. Supervised pKd spans 3.6253–8.3010 with mean 5.2361.

## Mapping tiers and structure views

All retained records are high-confidence S0 structure mappings at the measurement level: exact PDB, ligand, target, endpoint, and BindingDB record provenance. The feature table supplies two receptor views:

- **S0:** the query complex’s receptor coordinates, with query ligand atoms excluded from every feature. The site residue identities are defined through the frozen historical reference. S0 is still privileged because the receptor can carry query-induced conformation.
- **S1:** one earlier exact-construct historical receptor and one fixed residue set per group. Every query ligand in the group sees the same receptor features. This removes query-specific receptor conformation but remains a historical holo view, not apo/predicted deployment.

The historical ligand is used only once to identify residues within 6 Angstrom in the reference structure. Those residue identities are transferred by exact construct sequence. Reference affinity labels are blank and quarantined. Query ligand coordinates are never consumed.

## Feature boundary

Allowed ligand inputs are canonical SMILES and graph features derived from 2D connectivity. Query conformers, docked or co-crystal poses, ligand-centred grids, protein–ligand distances, atom maps, isotopes, and CXSMILES coordinate extensions are prohibited.

The Gate-0/1 protein representation is deliberately transparent rather than learned: 52 rigid-motion-invariant receptor descriptors covering residue/element composition and receptor-internal distance/shape summaries. This feature family is suitable for a geometry-value diagnostic, not a substitute for the proposed GVP/atom encoder.

## Leakage controls

The canonical split joins observations before splitting when any configured ligand, scaffold, construct, sequence, local-site, publication, or validated structure-similarity edge exists. Five disconnected union components result. No retained cross-construct pair exceeds the ECFP, full-sequence, or local-site thresholds.

The validated local pocket-structure edge is unavailable in this bounded pilot. Consequently, “zero component overlap” means zero overlap under the implemented graph, not proof that the data are globally leak-free.

## Known limitations

- five components are too few for confirmatory uncertainty or family-level generalization;
- all labels are structural-database-selected and can carry publication/chemotype selection bias;
- S0 and S1 retain different forms of holo privilege;
- source-paper experimental tables were not independently re-extracted;
- same-PDB mapping does not eliminate every ambiguity in constructs, protonation, oligomeric context, cofactors, or assay conditions;
- the pilot has no apo/predicted structure, temporal external test, or novel-pocket benchmark;
- gMolAI’s exact pretraining identity overlap is unknown; and
- corrected v2 labels must not be merged with the superseded v1 candidate table.

## Permitted and prohibited claims

Permitted: adapter validation, acquisition regression testing, split/leakage diagnostics, shallow-control comparison, and planning a higher-confidence data tranche.

Prohibited: state-of-the-art claims, external validity, prospective hit discovery, calibrated absolute affinity, general receptor-geometry benefit, or full-model architecture selection.
