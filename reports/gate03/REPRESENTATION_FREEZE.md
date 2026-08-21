# Gate-3 representation freeze

**Frozen:** 2026-08-21, before the first Gate-3 efficacy fit

**Contract:** `configs/gate03/representations-v1.yaml`

## Receptor views

- S1 historical fixed holo receptor: 34/34 series.
- S2 AlphaFold DB non-holo receptor: 32/34 series (94.12%), above the
  preregistered 80% requirement.
- S3 strict pocket-unoccupied X-ray receptor: 7/34 series, sensitivity only.
- Pocket descriptor v2: 116 rigid-invariant geometry, chemistry, radial-shell,
  topology, exposure/shape, missingness, and quality terms.

The two S2 exclusions are retained: one site has mean pLDDT 67.18 (<70), and one
viral AlphaFold record fails exact construct mapping. S3 selection uses exact
declared UniProt mapping, X-ray resolution <=2.8 A, zero nonpolymer entities,
and no foreign heavy atom within 8 A of the frozen site.

## Frozen pretrained encoders

ESM-2 `facebook/esm2_t30_150M_UR50D` is pinned to repository revision
`a695f6045e2e32885fa60af20c13cb35398ce30c` and checkpoint SHA-256
`c3f1da8aea53bddd32c246c86168c23b9fd72341fb9db9a94436f855f5053566`.
All 34 sequence vectors are finite, label-free, frozen, and bitwise repeat/order
deterministic.

ESM-IF1 `esm_if1_gvp4_t16_142M_UR50` uses the official checkpoint SHA-256
`be4ba36edec22a9bfaa4946ff6b2815f1f19d8a3d7e0eada8b796d5a0eae9fd4`.
Site-residue means and population standard deviations yield 1,024 dimensions.
The first uncentred encoding is retained as a failed audit: rotation was stable,
but a large translation caused 0.0142 maximum absolute drift. The selected v2
encoding deterministically centres finite backbone atoms and rounds to 0.0001 A
before inference. It is bitwise repeat/order deterministic and passes the same
rigid transform at `2.74e-5` maximum absolute difference.

## Ligand encoders

ECFP4/2,048 is primary and requires no pretrained corpus. The hash-pinned gMolAI
adapter accepted 657/661 observations (99.39%) and passed checkpoint/source,
atom-order mapping, dimensionality, repeat/order, equivalent-SMILES, and CPU/GPU
audits. Four disconnected ion-pair SMILES were rejected without salt stripping
or identity changes. Per amendment 03, gMolAI is a complete-case corroborating
sensitivity and cannot change the primary decision.

## Evidence

- `reports/gate03/evidence/receptor-view-audit-v1.json`
- `reports/gate03/evidence/apo-view-audit-v1.json`
- `reports/gate03/evidence/esm2-adapter-audit-v1.json`
- `reports/gate03/evidence/esm-if1-adapter-audit-v1.json` (failed, retained)
- `reports/gate03/evidence/esm-if1-adapter-audit-v2.json` (selected, pass)
- `reports/gate03/evidence/gmolai-adapter-audit-v1.json`

Raw model/receptor files and local packages remain ignored. Tracked manifests
contain their official URLs, versions, sizes, checksums, mapping ledgers, runtime
identity, and the hashes/shapes of every tracked feature array.
