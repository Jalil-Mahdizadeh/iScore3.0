"""Pinned ESM-2 construct encoder used only as a frozen Gate-2 baseline."""

from __future__ import annotations

from dataclasses import dataclass
from io import BytesIO
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
from typing import Any, Mapping, Sequence

import numpy as np

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)


class Esm2Error(RuntimeError):
    """Raised when the frozen ESM-2 encoding contract is violated."""


@dataclass(frozen=True, slots=True)
class ConstructSequence:
    group_id: str
    construct_sha256: str
    sequence: str


def _npy_bytes(value: np.ndarray) -> bytes:
    buffer = BytesIO()
    np.save(buffer, value, allow_pickle=False)
    return buffer.getvalue()


def _array_sha256(value: np.ndarray) -> str:
    contiguous = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(contiguous.dtype).encode("ascii"))
    digest.update(
        json.dumps(list(contiguous.shape), separators=(",", ":")).encode("ascii")
    )
    digest.update(contiguous.tobytes(order="C"))
    return digest.hexdigest()


def load_constructs(pilot: Path) -> tuple[list[ConstructSequence], int]:
    with pilot.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if rows and "role" in rows[0]:
        references = [row for row in rows if row["role"] == "site_reference_only"]
        if not references or any(
            row.get("pKd") or row.get("value_nm") for row in references
        ):
            raise Esm2Error("Historical reference label quarantine is absent or invalid")
        supervised = [row for row in rows if row["role"] == "supervised_s0"]
        fields = ("construct_group_id", "construct_sha256", "construct_sequence")
    elif rows and {"series_id", "target_sequence_sha256", "target_sequence"}.issubset(rows[0]):
        supervised = rows
        fields = ("series_id", "target_sequence_sha256", "target_sequence")
    else:
        raise Esm2Error("Unsupported construct table schema")
    by_group: dict[str, ConstructSequence] = {}
    for row in supervised:
        record = ConstructSequence(
            group_id=row[fields[0]],
            construct_sha256=row[fields[1]],
            sequence=row[fields[2]],
        )
        previous = by_group.setdefault(record.group_id, record)
        if previous != record:
            raise Esm2Error(f"Inconsistent construct definition: {record.group_id}")
    constructs = [by_group[group] for group in sorted(by_group)]
    if not constructs:
        raise Esm2Error("No supervised constructs")
    for record in constructs:
        if (
            hashlib.sha256(record.sequence.encode("ascii")).hexdigest()
            != record.construct_sha256
        ):
            raise Esm2Error(f"Construct sequence hash mismatch: {record.group_id}")
        if not record.sequence.isalpha() or record.sequence != record.sequence.upper():
            raise Esm2Error(f"Non-canonical sequence text: {record.group_id}")
    return constructs, len(supervised)


def acquire_snapshot(
    *, model_id: str, revision: str, cache_dir: Path, checkpoint_sha256: str
) -> tuple[Path, list[dict[str, Any]]]:
    from huggingface_hub import snapshot_download

    snapshot = Path(
        snapshot_download(
            repo_id=model_id,
            revision=revision,
            cache_dir=cache_dir,
            allow_patterns=[
                "config.json",
                "model.safetensors",
                "special_tokens_map.json",
                "tokenizer_config.json",
                "vocab.txt",
            ],
        )
    )
    checkpoint = snapshot / "model.safetensors"
    if not checkpoint.is_file() or sha256_file(checkpoint) != checkpoint_sha256:
        raise Esm2Error("Pinned ESM-2 checkpoint hash mismatch")
    files = []
    for path in sorted(value for value in snapshot.iterdir() if value.is_file()):
        files.append(
            {
                "name": path.name,
                "path": str(path.resolve()),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    return snapshot, files


class Esm2Encoder:
    """Frozen last-hidden-state residue mean with deterministic overlap handling."""

    def __init__(
        self,
        snapshot: Path,
        *,
        device: str,
        maximum_residues_per_window: int,
        overlap_residues: int,
        expected_hidden_dimension: int,
    ) -> None:
        import torch
        from transformers import AutoModel, AutoTokenizer

        if maximum_residues_per_window <= overlap_residues or overlap_residues < 0:
            raise Esm2Error("Invalid sequence-window contract")
        self.torch = torch
        self.maximum = maximum_residues_per_window
        self.overlap = overlap_residues
        self.device = torch.device(device)
        if self.device.type == "cuda" and not torch.cuda.is_available():
            raise Esm2Error("CUDA was requested but is unavailable")
        if self.device.type == "cuda" and os.environ.get(
            "CUBLAS_WORKSPACE_CONFIG"
        ) not in {
            ":4096:8",
            ":16:8",
        }:
            raise Esm2Error("CUDA deterministic mode requires CUBLAS_WORKSPACE_CONFIG")
        torch.manual_seed(20260820)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(20260820)
            torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        torch.use_deterministic_algorithms(True)
        torch.set_float32_matmul_precision("highest")
        self.tokenizer = AutoTokenizer.from_pretrained(
            snapshot, local_files_only=True, trust_remote_code=False
        )
        self.model = AutoModel.from_pretrained(
            snapshot,
            local_files_only=True,
            trust_remote_code=False,
            use_safetensors=True,
            add_pooling_layer=False,
        ).to(self.device)
        self.model.eval()
        self.model.requires_grad_(False)
        if int(self.model.config.hidden_size) != expected_hidden_dimension:
            raise Esm2Error(
                f"Unexpected hidden dimension: {self.model.config.hidden_size}"
            )
        self.hidden_dimension = expected_hidden_dimension

    def windows(self, length: int) -> list[tuple[int, int]]:
        if length <= self.maximum:
            return [(0, length)]
        stride = self.maximum - self.overlap
        starts = list(range(0, length - self.maximum + 1, stride))
        final = length - self.maximum
        if starts[-1] != final:
            starts.append(final)
        return [(start, min(length, start + self.maximum)) for start in starts]

    def encode(self, sequence: str) -> tuple[np.ndarray, dict[str, Any]]:
        torch = self.torch
        residue_sum = torch.zeros(
            (len(sequence), self.hidden_dimension), dtype=torch.float32, device="cpu"
        )
        residue_count = torch.zeros((len(sequence),), dtype=torch.int32, device="cpu")
        windows = self.windows(len(sequence))
        with torch.inference_mode():
            for start, stop in windows:
                segment = sequence[start:stop]
                tokens = self.tokenizer(
                    segment, return_tensors="pt", add_special_tokens=True, padding=False
                )
                input_ids = tokens["input_ids"].to(self.device)
                attention_mask = tokens["attention_mask"].to(self.device)
                if input_ids.shape[1] != len(segment) + 2:
                    raise Esm2Error("Tokenizer did not preserve one token per residue")
                hidden = self.model(
                    input_ids=input_ids, attention_mask=attention_mask
                ).last_hidden_state[0, 1 : len(segment) + 1]
                if hidden.shape != (len(segment), self.hidden_dimension):
                    raise Esm2Error(
                        f"Unexpected residue-state shape: {tuple(hidden.shape)}"
                    )
                residue_sum[start:stop] += hidden.float().cpu()
                residue_count[start:stop] += 1
        if bool(torch.any(residue_count == 0)):
            raise Esm2Error("Windowing left uncovered residues")
        averaged = residue_sum / residue_count.to(torch.float32).unsqueeze(1)
        embedding = averaged.mean(dim=0).numpy().astype(np.float32, copy=False)
        if (
            embedding.shape != (self.hidden_dimension,)
            or not np.isfinite(embedding).all()
        ):
            raise Esm2Error("Invalid pooled sequence embedding")
        return embedding, {
            "sequence_length": len(sequence),
            "window_count": len(windows),
            "windows": [[start, stop] for start, stop in windows],
            "minimum_residue_coverage": int(residue_count.min().item()),
            "maximum_residue_coverage": int(residue_count.max().item()),
        }

    def encode_many(
        self, constructs: Sequence[ConstructSequence]
    ) -> tuple[np.ndarray, list[dict[str, Any]]]:
        values = []
        records = []
        for record in constructs:
            embedding, windowing = self.encode(record.sequence)
            values.append(embedding)
            records.append(
                {
                    "construct_group_id": record.group_id,
                    "construct_sha256": record.construct_sha256,
                    **windowing,
                    "embedding_sha256": _array_sha256(embedding),
                }
            )
        return np.stack(values), records


def encode_sequences(
    *,
    pilot: Path,
    config: Mapping[str, Any],
    cache_dir: Path,
    feature_root: Path,
    manifest_path: Path,
    audit_path: Path,
    device: str,
) -> dict[str, Any]:
    import huggingface_hub
    import torch
    import transformers

    spec = config["sequence_encoder"]
    constructs, supervised_count = load_constructs(pilot)
    snapshot, snapshot_files = acquire_snapshot(
        model_id=str(spec["model"]),
        revision=str(spec["repository_revision"]),
        cache_dir=cache_dir,
        checkpoint_sha256=str(spec["checkpoint_sha256"]),
    )
    encoder = Esm2Encoder(
        snapshot,
        device=device,
        maximum_residues_per_window=int(spec["maximum_residues_per_window"]),
        overlap_residues=int(spec["window_overlap_residues"]),
        expected_hidden_dimension=int(spec["hidden_dimension"]),
    )
    embeddings, ledger = encoder.encode_many(constructs)
    panel_indices = sorted({0, len(constructs) // 2, len(constructs) - 1})
    panel = [constructs[index] for index in panel_indices]
    repeat, _ = encoder.encode_many(panel)
    reverse, _ = encoder.encode_many(list(reversed(panel)))
    reverse = reverse[::-1]
    expected = embeddings[panel_indices]
    repeat_max_abs = float(np.max(np.abs(expected - repeat)))
    order_max_abs = float(np.max(np.abs(expected - reverse)))

    arrays = {
        "construct_group_ids.npy": np.asarray(
            [record.group_id for record in constructs]
        ),
        "construct_sha256.npy": np.asarray(
            [record.construct_sha256 for record in constructs]
        ),
        "esm2_mean_last_hidden_state.npy": embeddings.astype(np.float32, copy=False),
    }
    array_files = []
    for name, value in arrays.items():
        payload = _npy_bytes(value)
        path = feature_root / name
        immutable_write(path, payload)
        array_files.append(
            {
                "name": name,
                "path": str(path.resolve()),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "dtype": str(value.dtype),
                "shape": list(value.shape),
            }
        )
    checks = {
        "all_construct_hashes_verified": True,
        "labels_not_passed_to_encoder": True,
        "reference_labels_quarantined": True,
        "model_eval_mode": not encoder.model.training,
        "all_parameters_frozen": not any(
            parameter.requires_grad for parameter in encoder.model.parameters()
        ),
        "hidden_dimension": embeddings.shape[1] == int(spec["hidden_dimension"]),
        "all_embeddings_finite": bool(np.isfinite(embeddings).all()),
        "same_sequence_repeat_bitwise_exact": repeat_max_abs == 0.0,
        "input_order_invariance_bitwise_exact": order_max_abs == 0.0,
        "checkpoint_hash_verified": True,
    }
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "pilot": {"path": str(pilot.resolve()), "sha256": sha256_file(pilot)},
        "model": {
            "repository": str(spec["model"]),
            "revision": str(spec["repository_revision"]),
            "snapshot_path": str(snapshot.resolve()),
            "checkpoint_sha256": str(spec["checkpoint_sha256"]),
            "files": snapshot_files,
            "hidden_dimension": int(spec["hidden_dimension"]),
            "layers": int(encoder.model.config.num_hidden_layers),
            "vocabulary_size": int(encoder.model.config.vocab_size),
            "pooling": str(spec["pooling"]),
            "fine_tuned": False,
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "torch": torch.__version__,
            "transformers": transformers.__version__,
            "huggingface_hub": huggingface_hub.__version__,
            "device": str(encoder.device),
            "cuda_device": (
                torch.cuda.get_device_name(encoder.device)
                if encoder.device.type == "cuda"
                else ""
            ),
            "deterministic_algorithms": torch.are_deterministic_algorithms_enabled(),
            "tf32_allowed": (
                bool(torch.backends.cuda.matmul.allow_tf32)
                if torch.cuda.is_available()
                else False
            ),
            "cublas_workspace_config": os.environ.get("CUBLAS_WORKSPACE_CONFIG", ""),
        },
        "counts": {
            "constructs": len(constructs),
            "supervised_observations": supervised_count,
            "minimum_sequence_length": min(
                len(record.sequence) for record in constructs
            ),
            "maximum_sequence_length": max(
                len(record.sequence) for record in constructs
            ),
            "multi_window_constructs": sum(row["window_count"] > 1 for row in ledger),
        },
        "window_contract": {
            "maximum_residues_per_window": int(spec["maximum_residues_per_window"]),
            "overlap_residues": int(spec["window_overlap_residues"]),
            "pooling": "overlap-average each residue, then unweighted residue mean",
        },
        "array_files": array_files,
        "construct_ledger": ledger,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    audit = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "overall_status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "repeat_maximum_absolute_delta": repeat_max_abs,
        "order_invariance_maximum_absolute_delta": order_max_abs,
        "panel_constructs": [record.group_id for record in panel],
        "manifest_sha256": sha256_file(manifest_path),
        "information_boundary": {
            "encoder_inputs": "construct amino-acid sequence only",
            "affinity_labels_read_by_encoder": False,
            "ligand_inputs_read_by_encoder": False,
            "fine_tuning": False,
        },
    }
    preserve_manifest_timestamp(audit_path, audit, "created_utc")
    immutable_write(audit_path, stable_json_bytes(audit))
    return audit
