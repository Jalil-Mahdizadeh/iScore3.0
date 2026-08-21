#!/usr/bin/env python3
"""Audit the pinned gMolAI adapter and encode the bounded pilot."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from io import BytesIO
import csv
import hashlib
import json
from pathlib import Path
import platform
from typing import Any, Mapping, Sequence

import numpy as np

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
)
from iscore3.ligand.gmolai_adapter import (
    GmolaiAdapter,
    GmolaiEncoding,
    array_sha256,
    compare_encodings,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def npy_bytes(value: np.ndarray) -> bytes:
    buffer = BytesIO()
    np.save(buffer, value, allow_pickle=False)
    return buffer.getvalue()


def read_pilot(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    references = [row for row in rows if row["role"] == "site_reference_only"]
    if not references or any(
        row.get("pKd") or row.get("value_nm") for row in references
    ):
        raise RuntimeError("Site-reference label quarantine is absent or invalid")
    supervised = [row for row in rows if row["role"] == "supervised_s0"]
    if len({row["observation_id"] for row in supervised}) != len(supervised):
        raise RuntimeError("Supervised observation IDs are not unique")
    return supervised


def exact_comparison(value: Mapping[str, Any]) -> bool:
    return all(
        value[name]["exact"] for name in ("node_z", "graph_z", "released_molecule_z")
    )


def within_cross_device_tolerance(value: Mapping[str, Any]) -> bool:
    return all(
        value[name]["maximum_absolute_delta"] <= 2.0e-4
        and value[name]["relative_l2"] <= 2.0e-5
        and value[name]["cosine"] >= 0.999999
        for name in ("node_z", "graph_z", "released_molecule_z")
    )


def encode_pilot(
    adapter: GmolaiAdapter,
    rows: Sequence[Mapping[str, str]],
    feature_root: Path,
) -> tuple[list[GmolaiEncoding], list[dict[str, Any]], list[dict[str, Any]]]:
    encodings = adapter.encode_many(row["canonical_smiles"] for row in rows)
    node_offsets = [0]
    node_rows: list[np.ndarray] = []
    canonical_to_input: list[int] = []
    input_to_canonical: list[int] = []
    ledger: list[dict[str, Any]] = []
    for row, encoding in zip(rows, encodings, strict=True):
        node_rows.append(encoding.node_z.astype(np.float32, copy=False))
        canonical_to_input.extend(encoding.mapping.canonical_to_input)
        input_to_canonical.extend(encoding.mapping.input_to_canonical)
        node_offsets.append(node_offsets[-1] + encoding.node_z.shape[0])
        ledger.append(
            {
                "observation_id": row["observation_id"],
                "inchikey": row["inchikey"],
                "input_smiles": row["canonical_smiles"],
                "canonical_smiles": encoding.canonical_smiles,
                "molecule_hash": encoding.molecule_hash,
                "atom_count": int(encoding.node_z.shape[0]),
                "canonical_to_input_atom": list(encoding.mapping.canonical_to_input),
                "input_to_canonical_atom": list(encoding.mapping.input_to_canonical),
                "canonical_atom_symbols": list(encoding.mapping.canonical_atom_symbols),
                "node_features_sha256": array_sha256(encoding.node_features),
                "edge_index_sha256": array_sha256(encoding.edge_index),
                "edge_features_sha256": array_sha256(encoding.edge_features),
                "node_z_sha256": array_sha256(encoding.node_z),
                "graph_z_sha256": array_sha256(encoding.graph_z),
                "released_molecule_z_sha256": array_sha256(
                    encoding.released_molecule_z
                ),
            }
        )

    arrays = {
        "observation_ids.npy": np.asarray([row["observation_id"] for row in rows]),
        "canonical_smiles.npy": np.asarray(
            [encoding.canonical_smiles for encoding in encodings]
        ),
        "node_offsets.npy": np.asarray(node_offsets, dtype=np.int64),
        "canonical_to_input_atom.npy": np.asarray(canonical_to_input, dtype=np.int64),
        "input_to_canonical_atom.npy": np.asarray(input_to_canonical, dtype=np.int64),
        "node_z.npy": np.concatenate(node_rows, axis=0).astype(np.float32, copy=False),
        "graph_z.npy": np.stack([encoding.graph_z for encoding in encodings]).astype(
            np.float32, copy=False
        ),
        "released_molecule_z.npy": np.stack(
            [encoding.released_molecule_z for encoding in encodings]
        ).astype(np.float32, copy=False),
    }
    files: list[dict[str, Any]] = []
    for name, value in arrays.items():
        path = feature_root / name
        payload = npy_bytes(value)
        immutable_write(path, payload)
        files.append(
            {
                "name": name,
                "path": str(path.resolve()),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
                "dtype": str(value.dtype),
                "shape": list(value.shape),
            }
        )
    return encodings, ledger, files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--adapter-config", type=Path, required=True)
    parser.add_argument("--container", type=Path, required=True)
    parser.add_argument("--feature-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--audit-report", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--expected-supervised-count", type=int, required=True)
    args = parser.parse_args()

    rows = read_pilot(args.pilot)
    gpu = GmolaiAdapter(args.source_root, args.adapter_config, device=args.device)
    encodings, ledger, feature_files = encode_pilot(gpu, rows, args.feature_root)

    panel = [
        "OCC",
        "F[C@@H](O)C",
        rows[0]["canonical_smiles"],
        rows[len(rows) // 2]["canonical_smiles"],
    ]
    first = gpu.encode_many(panel)
    repeated = gpu.encode_many(panel)
    reverse = list(reversed(gpu.encode_many(reversed(panel))))
    exact_repeat = [
        compare_encodings(a, b) for a, b in zip(first, repeated, strict=True)
    ]
    order_invariance = [
        compare_encodings(a, b) for a, b in zip(first, reverse, strict=True)
    ]
    equivalent_first = gpu.encode("CCO")
    equivalent_second = gpu.encode("OCC")
    equivalent_smiles = compare_encodings(equivalent_first, equivalent_second)

    cpu = GmolaiAdapter(args.source_root, args.adapter_config, device="cpu")
    cpu_panel = cpu.encode_many(panel)
    cross_device = [
        compare_encodings(a, b) for a, b in zip(first, cpu_panel, strict=True)
    ]

    checks = {
        "supervised_count_matches_explicit_contract": (
            len(rows) == args.expected_supervised_count
        ),
        "all_dimensions": all(
            encoding.node_z.shape[1] == 128
            and encoding.graph_z.shape == (256,)
            and encoding.released_molecule_z.shape == (384,)
            and encoding.node_features.shape[1] == 48
            and encoding.edge_features.shape[1] == 15
            for encoding in encodings
        ),
        "canonical_atom_rows_equal_node_rows": all(
            len(encoding.mapping.canonical_to_input) == encoding.node_z.shape[0]
            for encoding in encodings
        ),
        "exact_repeat": all(exact_comparison(value) for value in exact_repeat),
        "input_order_invariance": all(
            exact_comparison(value) for value in order_invariance
        ),
        "equivalent_smiles_same_canonical_node_states": exact_comparison(
            equivalent_smiles
        ),
        "cpu_gpu_within_frozen_tolerance": all(
            within_cross_device_tolerance(value) for value in cross_device
        ),
        "no_reference_labels_encoded": len(rows) == len(encodings),
    }

    import rdkit
    import torch
    import torch_geometric

    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "pilot": {"path": str(args.pilot.resolve()), "sha256": sha256_file(args.pilot)},
        "adapter_config": {
            "path": str(args.adapter_config.resolve()),
            "sha256": sha256_file(args.adapter_config),
        },
        "container": {
            "path": str(args.container.resolve()),
            "sha256": sha256_file(args.container),
        },
        "identity": gpu.identity_summary(),
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "torch": torch.__version__,
            "torch_geometric": torch_geometric.__version__,
            "rdkit": rdkit.__version__,
            "cuda_available": bool(torch.cuda.is_available()),
            "requested_device": args.device,
            "resolved_device": str(gpu.device),
        },
        "counts": {
            "molecules": len(encodings),
            "expected_supervised_count": args.expected_supervised_count,
            "total_atoms": int(sum(encoding.node_z.shape[0] for encoding in encodings)),
            "unique_canonical_smiles": len(
                {encoding.canonical_smiles for encoding in encodings}
            ),
            "unique_inchikeys": len({row["inchikey"] for row in rows}),
        },
        "array_files": feature_files,
        "molecule_ledger": ledger,
    }
    preserve_manifest_timestamp(args.manifest, manifest, "created_utc")
    immutable_write(args.manifest, stable_json_bytes(manifest))

    audit = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "overall_status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "tolerances": {
            "cross_device_maximum_absolute_delta": 2.0e-4,
            "cross_device_relative_l2": 2.0e-5,
            "cross_device_minimum_cosine": 0.999999,
            "same_device_repeat": "bitwise_exact",
        },
        "panel_smiles": panel,
        "exact_repeat": exact_repeat,
        "input_order_invariance": order_invariance,
        "equivalent_smiles": equivalent_smiles,
        "cpu_gpu_comparison": cross_device,
        "manifest_sha256": sha256_file(args.manifest),
        "interpretation": {
            "node_z": "raw clean-view atom states ordered by canonical RDKit atom index",
            "mapping": "validated bidirectional canonical-index to original-RDKit-index permutation",
            "pretraining_entity_exposure": "unknown and plausibly present for public pilot ligands",
            "affinity_label_exposure": "not detected in reviewed source/config; full corpus ledger unavailable",
            "authority": (
                "project PI explicitly identified himself as the gMolAI developer and authorized "
                "this pinned local scientific use"
            ),
        },
    }
    preserve_manifest_timestamp(args.audit_report, audit, "created_utc")
    immutable_write(args.audit_report, stable_json_bytes(audit))
    print(
        json.dumps(
            {"overall_status": audit["overall_status"], "checks": checks}, indent=2
        )
    )
    if audit["overall_status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
