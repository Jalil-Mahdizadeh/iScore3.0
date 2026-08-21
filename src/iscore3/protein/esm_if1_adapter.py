"""Frozen residue-level ESM-IF1 encoder for reusable receptor views."""

from __future__ import annotations

from io import BytesIO
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np

from iscore3.artifacts import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.protein.pocket_features import read_mmcif_atoms


class EsmIf1Error(RuntimeError):
    """Raised when the structure-encoding contract is violated."""


def _install_native_scatter_compatibility() -> None:
    """Provide only the torch-scatter reductions imported by official ESM-IF1."""

    if "torch_scatter" in sys.modules:
        return
    import torch

    module = ModuleType("torch_scatter")

    def scatter_add(src, index, dim=-1, out=None, dim_size=None):
        axis = dim if dim >= 0 else src.dim() + dim
        if out is None:
            shape = list(src.shape)
            inferred = int(index.max().item() + 1) if index.numel() else 0
            shape[axis] = int(dim_size if dim_size is not None else inferred)
            out = src.new_zeros(shape)
        view = [1] * src.dim()
        view[axis] = -1
        expanded = index.view(view).expand_as(src)
        return out.scatter_add_(axis, expanded, src)

    def scatter(src, index, dim=-1, out=None, dim_size=None, reduce="sum"):
        if reduce not in {"sum", "add"}:
            raise EsmIf1Error(f"Unsupported scatter reduction: {reduce}")
        return scatter_add(src, index, dim, out, dim_size)

    module.scatter_add = scatter_add
    module.scatter = scatter
    module.__iscore3_native_compatibility__ = True
    sys.modules["torch_scatter"] = module


def _tree_sha256(root: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    count = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        payload = path.read_bytes()
        digest.update(str(path.relative_to(root)).encode("utf-8"))
        digest.update(b"\0")
        digest.update(hashlib.sha256(payload).digest())
        count += 1
    return digest.hexdigest(), count


def _npy_bytes(value: np.ndarray) -> bytes:
    buffer = BytesIO()
    np.save(buffer, value, allow_pickle=False)
    return buffer.getvalue()


def _array_sha256(value: np.ndarray) -> str:
    value = np.ascontiguousarray(value)
    return hashlib.sha256(
        str(value.dtype).encode("ascii")
        + json.dumps(list(value.shape), separators=(",", ":")).encode("ascii")
        + value.tobytes()
    ).hexdigest()


def _backbone(
    structure_path: Path,
    entity_id: str,
    asym_id: str,
    site_positions: Sequence[int],
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    atoms = [
        atom
        for atom in read_mmcif_atoms(structure_path)
        if atom.group == "ATOM"
        and atom.entity_id == str(entity_id)
        and atom.asym_id == str(asym_id)
        and atom.seq_id is not None
    ]
    if not atoms:
        raise EsmIf1Error(f"No mapped protein atoms in {structure_path}")
    present = sorted({int(atom.seq_id) for atom in atoms})
    positions = list(range(min(present), max(present) + 1))
    index_by_position = {position: index for index, position in enumerate(positions)}
    coords = np.full((len(positions), 3, 3), np.nan, dtype=np.float32)
    atom_index = {"N": 0, "CA": 1, "C": 2}
    for atom in atoms:
        name = atom.atom_name.strip()
        if name in atom_index:
            coords[index_by_position[int(atom.seq_id)], atom_index[name]] = atom.xyz
    site_indices = np.asarray(
        [index_by_position[position] for position in site_positions if position in index_by_position],
        dtype=np.int64,
    )
    if len(site_indices) == 0:
        raise EsmIf1Error(f"No site positions occur in mapped chain {structure_path}")
    complete = np.isfinite(coords).all(axis=(1, 2))
    complete_site = site_indices[complete[site_indices]]
    coverage = len(complete_site) / len(site_positions)
    return coords, complete_site, {
        "coordinate_positions_first": positions[0],
        "coordinate_positions_last": positions[-1],
        "coordinate_rows": len(positions),
        "complete_backbone_rows": int(np.sum(complete)),
        "requested_site_residues": len(site_positions),
        "complete_site_backbone_residues": len(complete_site),
        "complete_site_backbone_fraction": coverage,
    }


class EsmIf1Encoder:
    """Official frozen ESM-IF1 encoder with explicit site pooling."""

    def __init__(self, *, device: str, expected_checkpoint: Path, expected_sha256: str):
        _install_native_scatter_compatibility()
        import esm
        import torch

        if sha256_file(expected_checkpoint) != expected_sha256:
            raise EsmIf1Error("ESM-IF1 checkpoint hash mismatch")
        self.torch = torch
        self.device = torch.device(device)
        if self.device.type == "cuda" and not torch.cuda.is_available():
            raise EsmIf1Error("CUDA requested but unavailable")
        if self.device.type == "cuda" and os.environ.get("CUBLAS_WORKSPACE_CONFIG") not in {
            ":4096:8",
            ":16:8",
        }:
            raise EsmIf1Error("Deterministic CUDA requires CUBLAS_WORKSPACE_CONFIG")
        torch.manual_seed(20260821)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(20260821)
            torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.use_deterministic_algorithms(True)
        torch.set_float32_matmul_precision("highest")
        self.model, self.alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        self.model = self.model.to(self.device).eval().requires_grad_(False)
        if int(self.model.args.encoder_embed_dim) != 512:
            raise EsmIf1Error("Unexpected ESM-IF1 encoder dimension")
        from esm.inverse_folding.util import CoordBatchConverter

        self.batch_converter = CoordBatchConverter(self.alphabet)

    def encode(self, coords: np.ndarray, site_indices: np.ndarray) -> np.ndarray:
        torch = self.torch
        # ESM-IF1 is analytically translation invariant, but subtraction of
        # large float32 coordinate offsets can be numerically amplified by the
        # deep GVP stack.  Canonical centering in float64 and sub-coordinate-
        # precision quantization removes that irrelevant frame dependence.
        working = coords.astype(np.float64, copy=True)
        finite = np.isfinite(working).all(axis=2)
        points = working[finite]
        if not len(points):
            raise EsmIf1Error("Structure contains no finite backbone coordinates")
        centre = np.mean(points, axis=0)
        working -= centre
        working = np.round(working, decimals=4).astype(np.float32)
        batch = [(working, None, None)]
        values, confidence, _, _, padding = self.batch_converter(
            batch, device=self.device
        )
        with torch.inference_mode():
            output = self.model.encoder(
                values, padding, confidence, return_all_hiddens=False
            )["encoder_out"][0][1:-1, 0]
        if output.shape != (len(coords), 512):
            raise EsmIf1Error(f"Unexpected residue output shape: {tuple(output.shape)}")
        selected = output[torch.as_tensor(site_indices, device=self.device)]
        pooled = torch.cat((selected.mean(dim=0), selected.std(dim=0, unbiased=False)))
        value = pooled.float().cpu().numpy()
        if value.shape != (1024,) or not np.isfinite(value).all():
            raise EsmIf1Error("Invalid ESM-IF1 pooled site representation")
        return value


def _load_views(paths: Sequence[Path]) -> list[dict[str, Any]]:
    records = []
    for path in paths:
        manifest = json.loads(path.read_text(encoding="utf-8"))
        records.extend(
            row
            for row in manifest["views"]
            if row.get("status") == "eligible" and row.get("view") in {"S1", "S2", "S3"}
        )
    keys = [(row["series_id"], row["view"]) for row in records]
    if len(keys) != len(set(keys)):
        raise EsmIf1Error("Duplicate series/view record")
    return sorted(records, key=lambda row: (row["series_id"], row["view"]))


def encode_structure_views(
    *,
    view_manifests: Sequence[Path],
    config: Mapping[str, Any],
    checkpoint: Path,
    package_root: Path,
    feature_root: Path,
    manifest_path: Path,
    audit_path: Path,
    device: str,
) -> dict[str, Any]:
    """Encode all eligible S1/S2/S3 views and audit determinism/invariance."""

    spec = config["structure_encoder"]
    records = _load_views(view_manifests)
    encoder = EsmIf1Encoder(
        device=device,
        expected_checkpoint=checkpoint,
        expected_sha256=str(spec["checkpoint_sha256"]),
    )
    values = []
    ledgers = []
    prepared = []
    for row in records:
        path = Path(row["source_structure_path"])
        if sha256_file(path) != row["source_structure_sha256"]:
            raise EsmIf1Error(f"Structure hash mismatch: {row['series_id']}/{row['view']}")
        coords, site_indices, metrics = _backbone(
            path,
            row["protein_entity_id"],
            row["protein_asym_id"],
            [int(value) for value in row["source_site_positions"]],
        )
        if metrics["complete_site_backbone_fraction"] < float(
            spec["minimum_complete_site_backbone_fraction"]
        ):
            raise EsmIf1Error(
                f"Insufficient site backbone: {row['series_id']}/{row['view']}"
            )
        embedding = encoder.encode(coords, site_indices)
        values.append(embedding)
        prepared.append((coords, site_indices))
        ledgers.append(
            {
                "series_id": row["series_id"],
                "view": row["view"],
                "source_structure_path": str(path),
                "source_structure_sha256": sha256_file(path),
                "protein_entity_id": row["protein_entity_id"],
                "protein_asym_id": row["protein_asym_id"],
                **metrics,
                "embedding_sha256": _array_sha256(embedding),
            }
        )
    matrix = np.stack(values).astype(np.float32, copy=False)

    panel_indices = sorted({0, len(records) // 2, len(records) - 1})
    repeated = np.stack(
        [encoder.encode(*prepared[index]) for index in panel_indices]
    )
    expected = matrix[panel_indices]
    reverse = np.stack(
        [encoder.encode(*prepared[index]) for index in reversed(panel_indices)]
    )[::-1]
    repeat_delta = float(np.max(np.abs(expected - repeated)))
    order_delta = float(np.max(np.abs(expected - reverse)))

    rotation = np.asarray(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=np.float32
    )
    coords, site_indices = prepared[panel_indices[0]]
    transformed = coords @ rotation.T + np.asarray([31.0, -17.0, 8.0], dtype=np.float32)
    transformed_value = encoder.encode(transformed, site_indices)
    rigid_delta = float(np.max(np.abs(expected[0] - transformed_value)))

    arrays = {
        "series_ids.npy": np.asarray([row["series_id"] for row in records]),
        "views.npy": np.asarray([row["view"] for row in records]),
        "esm_if1_site_mean_std.npy": matrix,
    }
    array_files = []
    for name, value in arrays.items():
        payload = _npy_bytes(value)
        path = feature_root / name
        immutable_write(path, payload)
        array_files.append(
            {
                "name": name,
                "path": str(path),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
                "dtype": str(value.dtype),
                "shape": list(value.shape),
            }
        )
    esm_hash, esm_files = _tree_sha256(package_root / "esm")
    biotite_hash, biotite_files = _tree_sha256(package_root / "biotite")
    import biotite
    import esm
    import torch

    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "view_manifests": [
            {"path": str(path), "sha256": sha256_file(path)} for path in view_manifests
        ],
        "model": {
            "name": spec["model"],
            "official_repository": spec["official_repository"],
            "official_checkpoint_url": spec["official_checkpoint_url"],
            "checkpoint_path": str(checkpoint),
            "checkpoint_sha256": sha256_file(checkpoint),
            "encoder_hidden_dimension": 512,
            "pooled_dimension": 1024,
            "pooling": spec["pooling"],
            "coordinate_preprocessing": spec["coordinate_preprocessing"],
            "fine_tuned": False,
        },
        "packages": {
            "fair_esm": {
                "version": esm.__version__, "tree_sha256": esm_hash, "files": esm_files
            },
            "biotite": {
                "version": biotite.__version__,
                "tree_sha256": biotite_hash,
                "files": biotite_files,
            },
            "torch_scatter": spec["native_torch_scatter_compatibility"],
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "torch": torch.__version__,
            "device": str(encoder.device),
            "cuda_device": torch.cuda.get_device_name(encoder.device),
            "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
            "tf32_allowed": bool(torch.backends.cuda.matmul.allow_tf32),
            "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG", ""),
        },
        "counts": {
            "views": len(records),
            "S1": sum(row["view"] == "S1" for row in records),
            "S2": sum(row["view"] == "S2" for row in records),
            "S3": sum(row["view"] == "S3" for row in records),
        },
        "array_files": array_files,
        "view_ledger": ledgers,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    checks = {
        "checkpoint_hash_verified": True,
        "all_parameters_frozen": not any(p.requires_grad for p in encoder.model.parameters()),
        "model_eval_mode": not encoder.model.training,
        "expected_dimensions": matrix.shape == (len(records), 1024),
        "all_embeddings_finite": bool(np.isfinite(matrix).all()),
        "same_view_repeat_bitwise_exact": repeat_delta == 0.0,
        "input_order_invariance_bitwise_exact": order_delta == 0.0,
        "rigid_transform_max_abs_below_2e_minus_4": rigid_delta <= 2.0e-4,
        "labels_not_passed_to_encoder": True,
        "query_ligand_coordinates_not_passed": True,
    }
    audit = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "overall_status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "determinism": {
            "repeat_maximum_absolute_delta": repeat_delta,
            "input_order_maximum_absolute_delta": order_delta,
            "rigid_transform_maximum_absolute_delta": rigid_delta,
        },
        "manifest_sha256": sha256_file(manifest_path),
    }
    preserve_manifest_timestamp(audit_path, audit, "created_utc")
    immutable_write(audit_path, stable_json_bytes(audit))
    return audit
