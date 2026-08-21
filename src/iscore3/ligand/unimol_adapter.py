"""Hash-pinned Uni-Mol v1 adapter for explicit free-ligand conformers."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any, Sequence

import numpy as np

from iscore3.gate4a.ligand3d import FreeConformerEnsemble


class UniMolAdapterError(RuntimeError):
    """Raised when source, weights, input alignment, or output contract changes."""


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git_identity(path: Path) -> tuple[str, bool]:
    revision = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    dirty = bool(
        subprocess.run(
            ["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    )
    return revision, dirty


class FrozenUniMolV1Adapter:
    """Molecule-all-H CLS representations from caller-supplied conformers only."""

    output_width = 512

    def __init__(
        self,
        source_root: str | Path,
        manifest_path: str | Path,
        *,
        use_cuda: bool = True,
        batch_size: int = 64,
    ) -> None:
        self.source_root = Path(source_root).resolve()
        self.manifest_path = Path(manifest_path).resolve()
        self.repository_root = self.manifest_path.parents[3]
        self.manifest = json.loads(self.manifest_path.read_text(encoding="utf-8"))
        revision, dirty = _git_identity(self.source_root)
        expected_revision = self.manifest["software"]["commit"]
        if revision != expected_revision or dirty:
            raise UniMolAdapterError(
                f"Uni-Mol source identity mismatch: revision={revision}, dirty={dirty}"
            )
        artifacts = {Path(row["path"]).name: row for row in self.manifest["artifacts"]}
        checkpoint_row = artifacts["mol_pre_all_h_220816.pt"]
        dictionary_row = artifacts["mol.dict.txt"]
        self.checkpoint_path = self.repository_root / checkpoint_row["path"]
        self.dictionary_path = self.repository_root / dictionary_row["path"]
        for path, row in (
            (self.checkpoint_path, checkpoint_row),
            (self.dictionary_path, dictionary_row),
        ):
            if not path.is_file() or path.stat().st_size != int(row["bytes"]):
                raise UniMolAdapterError(f"missing or size-mismatched artifact: {path}")
            if _sha256(path) != row["sha256"]:
                raise UniMolAdapterError(f"artifact hash mismatch: {path}")

        from unimol_tools import UniMolRepr

        self.model = UniMolRepr(
            data_type="molecule",
            batch_size=batch_size,
            remove_hs=False,
            model_name="unimolv1",
            use_cuda=use_cuda,
            pretrained_model_path=str(self.checkpoint_path),
            pretrained_dict_path=str(self.dictionary_path),
            max_atoms=256,
        )
        self.device = str(self.model.device)

    @staticmethod
    def explicit_inputs(ensemble: FreeConformerEnsemble) -> dict[str, list[Any]]:
        atoms = [list(ensemble.atom_symbols) for _ in ensemble.coordinates_angstrom]
        coordinates = [
            np.asarray(value, dtype=np.float32)
            for value in ensemble.coordinates_angstrom
        ]
        if any(len(atom_row) != len(coord_row) for atom_row, coord_row in zip(atoms, coordinates, strict=True)):
            raise UniMolAdapterError("atom-coordinate alignment failed")
        return {"atoms": atoms, "coordinates": coordinates}

    def encode_conformers(self, ensemble: FreeConformerEnsemble) -> np.ndarray:
        values = self.model.get_repr(self.explicit_inputs(ensemble))
        array = np.asarray(values, dtype=np.float32)
        expected = (len(ensemble.coordinates_angstrom), self.output_width)
        if array.shape != expected:
            raise UniMolAdapterError(f"unexpected Uni-Mol output shape {array.shape}; expected {expected}")
        if not np.isfinite(array).all():
            raise UniMolAdapterError("Uni-Mol produced non-finite output")
        return array

    def encode_inputs(
        self,
        atoms: Sequence[Sequence[str]],
        coordinates: Sequence[np.ndarray],
    ) -> np.ndarray:
        values = self.model.get_repr(
            {
                "atoms": [list(value) for value in atoms],
                "coordinates": [np.asarray(value, dtype=np.float32) for value in coordinates],
            }
        )
        array = np.asarray(values, dtype=np.float32)
        if array.shape != (len(atoms), self.output_width) or not np.isfinite(array).all():
            raise UniMolAdapterError("invalid Uni-Mol explicit-coordinate output")
        return array

    @staticmethod
    def aggregate(
        conformer_representations: np.ndarray, weights: Sequence[float]
    ) -> np.ndarray:
        array = np.asarray(conformer_representations, dtype=np.float64)
        weight_array = np.asarray(weights, dtype=np.float64)
        if array.ndim != 2 or array.shape[0] != len(weight_array):
            raise UniMolAdapterError("representations and weights do not align")
        weight_array = weight_array / weight_array.sum()
        mean = np.sum(array * weight_array[:, None], axis=0)
        variance = np.sum(np.square(array - mean) * weight_array[:, None], axis=0)
        return np.concatenate([mean, np.sqrt(np.maximum(variance, 0.0))]).astype(np.float32)

    def provenance(self) -> dict[str, Any]:
        return {
            "source_commit": self.manifest["software"]["commit"],
            "checkpoint_sha256": _sha256(self.checkpoint_path),
            "dictionary_sha256": _sha256(self.dictionary_path),
            "variant": "unimolv1_molecule_all_h",
            "output_width_per_conformer": self.output_width,
            "ensemble_output_width": self.output_width * 2,
            "device": self.device,
            "supervised_finetuning": False,
            "input_kind": "caller_supplied_free_conformer_atoms_and_coordinates",
        }
