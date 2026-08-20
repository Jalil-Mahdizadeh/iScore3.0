"""Fail-closed adapter for the frozen gMolAI-v2.0 atom and molecule embeddings.

The adapter intentionally imports, rather than copies, the pinned upstream
canonicalizer, featurizer, model definition, and checkpoint loader. It adds the
atom-index ledger that the public molecule-vector CLI does not expose.

``node_z[i]`` always corresponds to atom ``i`` of ``canonical_smiles`` as parsed
by the pinned RDKit runtime. ``canonical_to_input_atom[i]`` maps that row back to
the RDKit atom index obtained by parsing the caller's original SMILES. This is a
graph mapping, not a source-string character offset.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import ModuleType
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml


class AdapterContractError(RuntimeError):
    """Raised when source, artifact, preprocessing, or output identity drifts."""


@dataclass(frozen=True, slots=True)
class AtomMapping:
    """Bidirectional RDKit atom-index mapping for one accepted SMILES."""

    canonical_to_input: tuple[int, ...]
    input_to_canonical: tuple[int, ...]
    input_atom_symbols: tuple[str, ...]
    canonical_atom_symbols: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class GmolaiEncoding:
    """One molecule encoded under the frozen Gate-0/1 contract."""

    input_smiles: str
    canonical_smiles: str
    molecule_hash: str
    mapping: AtomMapping
    node_features: np.ndarray
    edge_index: np.ndarray
    edge_features: np.ndarray
    node_z: np.ndarray
    graph_z: np.ndarray
    released_molecule_z: np.ndarray
    provenance: Mapping[str, Any]


def _sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AdapterContractError(f"Adapter configuration is not a mapping: {path}")
    if value.get("schema_version") != 1:
        raise AdapterContractError("Unsupported gMolAI adapter configuration schema")
    return value


def _git_revision(source_root: Path) -> tuple[str, bool]:
    try:
        revision = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        dirty = bool(
            subprocess.run(
                [
                    "git",
                    "-C",
                    str(source_root),
                    "status",
                    "--porcelain",
                    "--untracked-files=no",
                    "--",
                    "src/gmolai_retrain",
                    "inference/generate_embeddings.py",
                ],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise AdapterContractError(f"Cannot verify gMolAI Git identity: {error}") from error
    return revision, dirty


def _parse_output_order(molecule: Any) -> tuple[int, ...]:
    if not molecule.HasProp("_smilesAtomOutputOrder"):
        raise AdapterContractError("RDKit did not expose canonical SMILES atom output order")
    raw = molecule.GetProp("_smilesAtomOutputOrder")
    try:
        order = tuple(int(value) for value in ast.literal_eval(raw))
    except (SyntaxError, ValueError, TypeError) as error:
        raise AdapterContractError(f"Invalid RDKit atom output order {raw!r}") from error
    return order


def _atom_identity(atom: Any) -> tuple[int, int, int, bool]:
    return (
        int(atom.GetAtomicNum()),
        int(atom.GetFormalCharge()),
        int(atom.GetIsotope()),
        bool(atom.GetIsAromatic()),
    )


def canonical_atom_mapping(
    input_smiles: str,
    canonical_smiles: str,
    *,
    isomeric_smiles: bool = True,
) -> AtomMapping:
    """Map atoms of the upstream canonical reparse back to input parse indices.

    RDKit records the atom traversal used by ``MolToSmiles`` in the private but
    stable ``_smilesAtomOutputOrder`` molecule property. We validate that ledger
    against every atom and bond before exposing it. A failure is fatal; the
    adapter never guesses among graph automorphisms.
    """

    try:
        from rdkit import Chem
    except ImportError as error:  # pragma: no cover - environment specific
        raise AdapterContractError("RDKit is required for gMolAI atom mapping") from error

    input_molecule = Chem.MolFromSmiles(input_smiles)
    if input_molecule is None:
        raise AdapterContractError("RDKit cannot parse input SMILES for atom mapping")
    observed_canonical = Chem.MolToSmiles(
        input_molecule,
        canonical=True,
        isomericSmiles=bool(isomeric_smiles),
    )
    if observed_canonical != canonical_smiles:
        raise AdapterContractError(
            "Atom-mapping canonicalization differs from the accepted upstream identity: "
            f"{observed_canonical!r} != {canonical_smiles!r}"
        )
    canonical_to_input = _parse_output_order(input_molecule)
    canonical_molecule = Chem.MolFromSmiles(canonical_smiles)
    if canonical_molecule is None:
        raise AdapterContractError("RDKit cannot reparse accepted canonical SMILES")
    atom_count = int(canonical_molecule.GetNumAtoms())
    if len(canonical_to_input) != atom_count or set(canonical_to_input) != set(range(atom_count)):
        raise AdapterContractError("Canonical atom order is not a complete input-atom permutation")

    for canonical_index, input_index in enumerate(canonical_to_input):
        canonical_atom = canonical_molecule.GetAtomWithIdx(canonical_index)
        input_atom = input_molecule.GetAtomWithIdx(input_index)
        if _atom_identity(canonical_atom) != _atom_identity(input_atom):
            raise AdapterContractError(
                f"Atom identity mismatch at canonical atom {canonical_index} / input atom {input_index}"
            )

    for bond in canonical_molecule.GetBonds():
        first = canonical_to_input[int(bond.GetBeginAtomIdx())]
        second = canonical_to_input[int(bond.GetEndAtomIdx())]
        input_bond = input_molecule.GetBondBetweenAtoms(first, second)
        if input_bond is None or str(input_bond.GetBondType()) != str(bond.GetBondType()):
            raise AdapterContractError("Canonical atom ledger does not preserve molecular bonds")

    input_to_canonical_list = [-1] * atom_count
    for canonical_index, input_index in enumerate(canonical_to_input):
        input_to_canonical_list[input_index] = canonical_index
    return AtomMapping(
        canonical_to_input=canonical_to_input,
        input_to_canonical=tuple(input_to_canonical_list),
        input_atom_symbols=tuple(atom.GetSymbol() for atom in input_molecule.GetAtoms()),
        canonical_atom_symbols=tuple(atom.GetSymbol() for atom in canonical_molecule.GetAtoms()),
    )


def _import_upstream_api(source_root: Path) -> ModuleType:
    source_directory = source_root / "src"
    inference_directory = source_root / "inference"
    for directory in (source_directory, inference_directory):
        if str(directory) not in sys.path:
            sys.path.insert(0, str(directory))
    path = inference_directory / "generate_embeddings.py"
    module_name = "_iscore3_pinned_gmolai_generate_embeddings"
    existing = sys.modules.get(module_name)
    if existing is not None:
        existing_path = Path(str(getattr(existing, "__file__", ""))).resolve()
        if existing_path != path.resolve():
            raise AdapterContractError("A different gMolAI encoder API is already imported")
        return existing
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise AdapterContractError(f"Cannot construct an import spec for {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    try:
        spec.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        raise
    return module


class GmolaiAdapter:
    """Atom-level and released-vector access to one exact gMolAI release."""

    def __init__(
        self,
        source_root: str | Path,
        config_path: str | Path,
        *,
        device: str = "cpu",
    ) -> None:
        self.source_root = Path(source_root).resolve()
        self.config_path = Path(config_path).resolve()
        self.spec = _load_yaml(self.config_path)
        self._verify_static_identity()

        # Required before torch initializes cuBLAS. The public release CLI uses
        # the same value and deterministic settings.
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        self.api = _import_upstream_api(self.source_root)
        import torch

        torch.set_float32_matmul_precision("highest")
        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.benchmark = False
        self.torch = torch
        self.device = self.api.resolve_device(device)
        model_directory = self.source_root / self.spec["artifacts"]["directory"]
        self.bundle = self.api.load_model_bundle(model_directory, self.device)
        self._verify_loaded_contract()

    def _verify_static_identity(self) -> None:
        source = self.spec["source"]
        expected_revision = str(source["revision"])
        observed_revision, dirty = _git_revision(self.source_root)
        if observed_revision != expected_revision:
            raise AdapterContractError(
                f"gMolAI revision mismatch: expected {expected_revision}, observed {observed_revision}"
            )
        if dirty:
            raise AdapterContractError("Pinned gMolAI source contains tracked modifications")

        for relative, expected in source["required_file_sha256"].items():
            path = self.source_root / str(relative)
            if not path.is_file():
                raise AdapterContractError(f"Required gMolAI source file is missing: {path}")
            observed = _sha256(path)
            if observed != str(expected):
                raise AdapterContractError(
                    f"gMolAI source hash mismatch at {relative}: {observed} != {expected}"
                )
        model_directory = self.source_root / self.spec["artifacts"]["directory"]
        for name, expected in self.spec["artifacts"]["sha256"].items():
            path = model_directory / str(name)
            if not path.is_file():
                raise AdapterContractError(f"Required gMolAI artifact is missing: {path}")
            observed = _sha256(path)
            if observed != str(expected):
                raise AdapterContractError(
                    f"gMolAI artifact hash mismatch at {name}: {observed} != {expected}"
                )

    def _verify_loaded_contract(self) -> None:
        expected = self.spec["expected_contract"]
        observed = {
            "node_input_dimensions": int(self.bundle.model.feature_schema["node_input_dim"]),
            "edge_dimensions": int(self.bundle.model.feature_schema["edge_dim"]),
            "node_z_dimensions": int(self.bundle.model.node_latent_dim),
            "graph_z_dimensions": int(self.bundle.model.graph_latent_dim),
            "released_molecule_dimensions": int(self.bundle.embedding_dimensions),
            "checkpoint_global_step": int(self.bundle.checkpoint["global_step"]),
            "training_implementation_version": str(
                self.bundle.checkpoint["training_implementation_version"]
            ),
        }
        for key, value in observed.items():
            if value != expected[key]:
                raise AdapterContractError(
                    f"Loaded gMolAI contract mismatch at {key}: {value!r} != {expected[key]!r}"
                )
        if float(self.bundle.mean_node_weight) != float(expected["mean_node_weight"]):
            raise AdapterContractError("Loaded gMolAI mean-node weight changed")

    def _strict_precheck(self, smiles: str) -> None:
        from rdkit import Chem

        policy = self.spec["strict_input_policy"]
        if policy["reject_cxsmiles_extensions"] and "|" in smiles:
            raise AdapterContractError("CXSMILES extensions are outside the strict ligand boundary")
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return  # Upstream canonicalizer will return its stable reason.
        if policy["reject_input_conformers"] and molecule.GetNumConformers():
            raise AdapterContractError("Input SMILES unexpectedly materialized a conformer")
        if policy["reject_atom_maps"] and any(atom.GetAtomMapNum() for atom in molecule.GetAtoms()):
            raise AdapterContractError("Atom-map annotations are outside the strict ligand boundary")
        if policy["reject_isotopes"] and any(atom.GetIsotope() for atom in molecule.GetAtoms()):
            raise AdapterContractError("Isotope annotations are outside the frozen gMolAI feature schema")

    def encode(self, smiles: str) -> GmolaiEncoding:
        """Encode one SMILES and return a fully mapped, hash-bound record."""

        self._strict_precheck(smiles)
        canonical, reason = self.api.canonicalize_input(smiles, self.bundle.resolved_config)
        if reason is not None:
            raise AdapterContractError(f"gMolAI rejected the molecule: {reason}")
        assert canonical is not None
        policy = self.bundle.resolved_config["data"]["canonicalization"]
        mapping = canonical_atom_mapping(
            smiles,
            canonical.smiles,
            isomeric_smiles=bool(policy["isomeric_smiles"]),
        )

        from rdkit import Chem
        from gmolai_retrain.chem import featurize_molecule

        molecule = Chem.MolFromSmiles(canonical.smiles)
        if molecule is None or molecule.GetNumConformers() != 0:
            raise AdapterContractError("Canonical gMolAI molecule is invalid or contains coordinates")
        node_features, edge_index, edge_features = featurize_molecule(
            molecule,
            include_chirality=True,
            position_dim=0,
        )
        torch = self.torch
        x = torch.from_numpy(node_features).to(self.device)
        edges = torch.from_numpy(edge_index).to(self.device)
        edge_attr = torch.from_numpy(edge_features).to(self.device)
        batch = torch.zeros(len(node_features), dtype=torch.long, device=self.device)
        with torch.inference_mode():
            node_z, graph_z = self.bundle.model.encode(x, edges, edge_attr, batch)
            mean_node_z = node_z.mean(dim=0, keepdim=True)
            raw = torch.cat((graph_z, mean_node_z), dim=1)
            calibrated = self.bundle.model.apply_molecule_calibration(
                raw,
                self.bundle.coordinate_mean,
                self.bundle.coordinate_scale,
            )
            calibrated[:, int(self.bundle.graph_dimensions) :] *= float(
                self.bundle.mean_node_weight
            )

        expected = self.spec["expected_contract"]
        if tuple(node_z.shape) != (len(mapping.canonical_to_input), expected["node_z_dimensions"]):
            raise AdapterContractError(f"Unexpected node_z shape: {tuple(node_z.shape)}")
        if tuple(graph_z.shape) != (1, expected["graph_z_dimensions"]):
            raise AdapterContractError(f"Unexpected graph_z shape: {tuple(graph_z.shape)}")
        if tuple(calibrated.shape) != (1, expected["released_molecule_dimensions"]):
            raise AdapterContractError(
                f"Unexpected released molecule-vector shape: {tuple(calibrated.shape)}"
            )

        arrays = (
            node_features,
            edge_index,
            edge_features,
            node_z.detach().float().cpu().numpy(),
            graph_z.detach().float().cpu().numpy()[0],
            calibrated.detach().float().cpu().numpy()[0],
        )
        if not all(np.isfinite(value).all() for value in arrays):
            raise AdapterContractError("gMolAI adapter produced non-finite output")
        revision, dirty = _git_revision(self.source_root)
        provenance = {
            "source_url": self.spec["source"]["url"],
            "source_revision": revision,
            "source_dirty": dirty,
            "checkpoint_sha256": self.spec["artifacts"]["sha256"]["representation-best.pt"],
            "calibrator_sha256": self.spec["artifacts"]["sha256"][
                "representation-calibrator.pt"
            ],
            "selection_sha256": self.spec["artifacts"]["sha256"][
                "representation_selection.json"
            ],
            "feature_schema_hash": self.bundle.model.feature_schema["hash"],
            "checkpoint_global_step": int(self.bundle.checkpoint["global_step"]),
            "device": str(self.device),
            "torch": str(torch.__version__),
            "rdkit": str(Chem.rdBase.rdkitVersion),
            "node_z_definition": "raw deterministic clean-view atom states before molecule calibration",
            "released_molecule_definition": self.api.PUBLIC_EMBEDDING_DEFINITION,
            "pretraining_exact_entity_exposure": "unknown_likely_for_pubchem_compounds",
            "pretraining_affinity_label_exposure": "not_detected_but_corpus_lineage_incomplete",
        }
        return GmolaiEncoding(
            input_smiles=smiles,
            canonical_smiles=canonical.smiles,
            molecule_hash=canonical.molecule_hash,
            mapping=mapping,
            node_features=np.asarray(arrays[0]),
            edge_index=np.asarray(arrays[1]),
            edge_features=np.asarray(arrays[2]),
            node_z=np.asarray(arrays[3]),
            graph_z=np.asarray(arrays[4]),
            released_molecule_z=np.asarray(arrays[5]),
            provenance=provenance,
        )

    def encode_many(self, smiles: Iterable[str]) -> list[GmolaiEncoding]:
        """Audit-oriented bounded encoder; preserves caller order exactly."""

        return [self.encode(value) for value in smiles]

    def identity_summary(self) -> dict[str, Any]:
        revision, dirty = _git_revision(self.source_root)
        return {
            "source_url": self.spec["source"]["url"],
            "source_revision": revision,
            "source_dirty": dirty,
            "licence_status": self.spec["source"]["licence_status"],
            "artifact_sha256": dict(self.spec["artifacts"]["sha256"]),
            "expected_contract": dict(self.spec["expected_contract"]),
            "device": str(self.device),
            "pretraining_exposure": dict(self.spec["pretraining_exposure"]),
        }


def array_sha256(values: np.ndarray) -> str:
    """Stable hash over dtype, shape, and C-order bytes for audit ledgers."""

    contiguous = np.ascontiguousarray(values)
    header = json.dumps(
        {"dtype": str(contiguous.dtype), "shape": list(contiguous.shape)},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(header + b"\0" + contiguous.tobytes(order="C")).hexdigest()


def compare_encodings(first: GmolaiEncoding, second: GmolaiEncoding) -> dict[str, Any]:
    """Numerical comparison used by the Gate audit script."""

    if first.canonical_smiles != second.canonical_smiles:
        raise AdapterContractError("Cannot compare encodings with different canonical identities")

    def metrics(left: np.ndarray, right: np.ndarray) -> dict[str, Any]:
        if left.shape != right.shape:
            raise AdapterContractError(f"Array shapes differ: {left.shape} != {right.shape}")
        delta = left.astype(np.float64) - right.astype(np.float64)
        left_flat = left.astype(np.float64).ravel()
        right_flat = right.astype(np.float64).ravel()
        denominator = np.linalg.norm(left_flat) * np.linalg.norm(right_flat)
        cosine = 1.0 if denominator == 0 else float(np.dot(left_flat, right_flat) / denominator)
        return {
            "exact": bool(np.array_equal(left, right)),
            "maximum_absolute_delta": float(np.max(np.abs(delta), initial=0.0)),
            "relative_l2": float(np.linalg.norm(delta) / max(np.linalg.norm(left_flat), 1.0e-12)),
            "cosine": cosine,
            "first_sha256": array_sha256(left),
            "second_sha256": array_sha256(right),
        }

    return {
        "canonical_smiles": first.canonical_smiles,
        "node_z": metrics(first.node_z, second.node_z),
        "graph_z": metrics(first.graph_z, second.graph_z),
        "released_molecule_z": metrics(first.released_molecule_z, second.released_molecule_z),
    }
