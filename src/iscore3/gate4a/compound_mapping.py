"""Candidate-only compound identity mapping for the Davis source panel."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Mapping, Sequence


_DERIVATIVE = re.compile(r"\bderivative\b", flags=re.IGNORECASE)


@dataclass(frozen=True)
class CompoundIdentity:
    source_row: int
    source_name: str
    alternative_name: str


def candidate_queries(identity: CompoundIdentity) -> tuple[str, ...]:
    """Return conservative PubChem queries, never erasing a derivative qualifier."""

    if _DERIVATIVE.search(identity.source_name):
        return ()
    candidates = [identity.alternative_name.strip(), identity.source_name.strip()]
    candidates.extend(part.strip() for part in identity.source_name.split("/"))
    output: list[str] = []
    for candidate in candidates:
        if candidate and candidate not in output:
            output.append(candidate)
    return tuple(output)


def extract_pubchem_properties(response: Mapping[str, Any]) -> tuple[dict[str, Any], ...]:
    table = response.get("PropertyTable")
    if not isinstance(table, Mapping):
        return ()
    properties = table.get("Properties")
    if not isinstance(properties, Sequence) or isinstance(properties, (str, bytes)):
        return ()
    output = []
    for property_record in properties:
        if not isinstance(property_record, Mapping):
            continue
        if "CID" not in property_record or "SMILES" not in property_record:
            continue
        output.append(dict(property_record))
    return tuple(output)


def candidate_mapping_state(
    identity: CompoundIdentity,
    properties: Sequence[Mapping[str, Any]],
) -> str:
    if _DERIVATIVE.search(identity.source_name):
        return "quarantined_ambiguous_source_identity"
    if not properties:
        return "unresolved"
    return "candidate_requires_manual_verification"
