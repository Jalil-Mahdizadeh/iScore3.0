"""Small, strict provenance primitives shared by Gate-4A workflows."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any


class ProvenanceError(RuntimeError):
    """Raised when a source file violates its immutable manifest contract."""


@dataclass(frozen=True)
class VerifiedFile:
    source_id: str
    path: str
    bytes: int
    sha256: str


def sha256_file(path: str | Path, chunk_bytes: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of *path* without loading it into memory."""

    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: str | Path) -> dict[str, Any]:
    value = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ProvenanceError(f"manifest must contain a JSON object: {path}")
    return value


def verify_source_manifest(
    manifest_path: str | Path,
    *,
    repository_root: str | Path,
) -> tuple[VerifiedFile, ...]:
    """Verify every manifest file is regular, in-root, byte-exact, and hash-exact."""

    root = Path(repository_root).resolve()
    manifest = load_json(manifest_path)
    entries = manifest.get("files")
    if not isinstance(entries, list) or not entries:
        raise ProvenanceError("manifest.files must be a non-empty list")

    verified: list[VerifiedFile] = []
    seen_ids: set[str] = set()
    seen_paths: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise ProvenanceError("each manifest file entry must be an object")
        source_id = str(entry.get("source_id", "")).strip()
        relative_path = str(entry.get("path", "")).strip()
        if not source_id or source_id in seen_ids:
            raise ProvenanceError(f"missing or duplicate source_id: {source_id!r}")
        if not relative_path or relative_path in seen_paths:
            raise ProvenanceError(f"missing or duplicate source path: {relative_path!r}")

        candidate = root / relative_path
        if candidate.is_symlink():
            raise ProvenanceError(f"manifest source cannot be a symbolic link: {relative_path}")
        try:
            resolved = candidate.resolve(strict=True)
            resolved.relative_to(root)
        except (FileNotFoundError, ValueError) as exc:
            raise ProvenanceError(f"source is missing or outside repository: {relative_path}") from exc
        if not resolved.is_file():
            raise ProvenanceError(f"source is not a regular file: {relative_path}")

        expected_bytes = int(entry.get("bytes", -1))
        observed_bytes = resolved.stat().st_size
        if observed_bytes != expected_bytes:
            raise ProvenanceError(
                f"byte-size mismatch for {source_id}: expected {expected_bytes}, got {observed_bytes}"
            )
        expected_sha = str(entry.get("sha256", "")).lower()
        observed_sha = sha256_file(resolved)
        if observed_sha != expected_sha:
            raise ProvenanceError(
                f"SHA-256 mismatch for {source_id}: expected {expected_sha}, got {observed_sha}"
            )

        seen_ids.add(source_id)
        seen_paths.add(relative_path)
        verified.append(
            VerifiedFile(
                source_id=source_id,
                path=relative_path,
                bytes=observed_bytes,
                sha256=observed_sha,
            )
        )
    return tuple(verified)
