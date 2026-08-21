"""Gate-neutral immutable artifact utilities."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

from iscore3.provenance import sha256_file


class ArtifactError(RuntimeError):
    """Raised when an immutable artifact contract is violated."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def stable_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def immutable_write(path: Path, payload: bytes) -> str:
    """Write once, or verify that an existing snapshot is byte-identical."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    expected = sha256_bytes(payload)
    if path.exists():
        observed = sha256_file(path)
        if observed != expected:
            raise ArtifactError(f"refusing to replace nonidentical immutable file: {path}")
        return observed
    temporary = path.with_name(f".{path.name}.partial")
    if temporary.exists():
        raise ArtifactError(f"stale partial output requires review: {temporary}")
    with temporary.open("xb") as handle:
        handle.write(payload)
        handle.flush()
    temporary.replace(path)
    return expected


def preserve_manifest_timestamp(path: Path, manifest: dict[str, Any], field: str) -> None:
    """Keep write-once manifests byte-stable when a command is rerun."""

    path = Path(path)
    if not path.exists():
        return
    previous = json.loads(path.read_text(encoding="utf-8"))
    if field not in previous:
        raise ArtifactError(f"existing manifest lacks required timestamp {field}: {path}")
    manifest[field] = previous[field]
