#!/usr/bin/env python3
"""Write a deterministic hash inventory for the Gate-4A closure snapshot."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


TRACKED_GLOBS = (
    "configs/gate4a/protocol-v1.yaml",
    "data/manifests/gate4a/alphafold-pocket-coordinates-v1.json",
    "data/manifests/gate4a/apo-view-coordinates-v1.json",
    "data/manifests/gate4a/idg-dream-noise-evidence-v1.json",
    "data/manifests/gate4a/site-unoccupied-candidates-v1.json",
    "data/processed/gate4a/alphafold-pocket-admission-v1.tsv",
    "data/processed/gate4a/alphafold-pocket-structural-similarity-v1.tsv",
    "data/processed/gate4a/apo-view-admission-v1.tsv",
    "data/processed/gate4a/confirmation/*.tsv",
    "data/processed/gate4a/davis-ligand-secondary-qa-packet-v1.tsv",
    "data/splits/gate4a/davis-receptor-components-final-v1.tsv",
    "data/splits/gate4a/receptor-structural-leakage-edges-v1.tsv",
    "environments/gate4a-requirements.lock.txt",
    "reports/gate4a/GATE4A_PROVENANCE_CLOSURE_REPORT.md",
    "reports/gate4a/evidence/alphafold-pocket-admission-v1.json",
    "reports/gate4a/evidence/apo-view-admission-v1.json",
    "reports/gate4a/evidence/confirmation-ledger-freeze-v1.json",
    "reports/gate4a/evidence/davis-ligand-identity-freeze-v1.json",
    "reports/gate4a/evidence/pocket-provenance-contract-v2.json",
    "reports/gate4a/evidence/practical-equivalence-evidence-v2.json",
    "reports/gate4a/evidence/provenance-closure-replay-v1.json",
    "reports/gate4a/evidence/receptor-structural-leakage-v1.json",
    "scripts/gate4a/acquire_idg_dream_noise_evidence.py",
    "scripts/gate4a/acquire_site_unoccupied_candidates.py",
    "scripts/gate4a/audit_practical_equivalence_evidence.py",
    "scripts/gate4a/close_alphafold_receptors.py",
    "scripts/gate4a/finalize_structural_leakage.py",
    "scripts/gate4a/freeze_confirmation_ledgers.py",
    "scripts/gate4a/freeze_ligand_identity_qa.py",
    "scripts/gate4a/freeze_provenance_closure_manifest.py",
    "scripts/gate4a/qualify_apo_tiers.py",
    "src/iscore3/gate4a/receptor_closure.py",
    "tests/unit/gate4a/test_receptor_closure.py",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/gate4a/evidence/provenance-closure-reproducibility-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    paths = sorted(
        {
            path
            for pattern in TRACKED_GLOBS
            for path in root.glob(pattern)
            if path.is_file()
        }
    )
    missing = [pattern for pattern in TRACKED_GLOBS if not list(root.glob(pattern))]
    if missing:
        raise FileNotFoundError(f"missing closure artifacts: {missing}")
    payload = {
        "schema_version": 1,
        "phase": "gate4a_provenance_closure",
        "snapshot_date": "2026-08-21",
        "model_training_performed": False,
        "environment": "reports/gate4a/evidence/environment-v1.json",
        "artifacts": [
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in paths
        ],
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
