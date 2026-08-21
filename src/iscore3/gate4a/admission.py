"""Frozen, outcome-blind Gate-4A dataset-admission rules."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Iterable, Mapping, Sequence


DAVIS_PDF_PAGE_BY_SOURCE_ROW = {
    **{row: 10 for row in range(2, 9)},
    **{row: 11 for row in range(9, 16)},
    **{row: 12 for row in range(16, 23)},
    **{row: 13 for row in range(23, 30)},
    **{row: 14 for row in range(30, 37)},
    **{row: 15 for row in range(37, 44)},
    **{row: 16 for row in range(44, 51)},
    **{row: 17 for row in range(51, 58)},
    **{row: 18 for row in range(58, 65)},
    **{row: 19 for row in range(65, 72)},
    **{row: 20 for row in range(72, 74)},
}

STANDARD_EPK_GROUPS = frozenset({"AGC", "CAMK", "CK1", "CMGC", "Other", "STE", "TK", "TKL"})
EXPLICIT_SOURCE_STEREO_ROWS = frozenset(
    {6, 13, 15, 19, 20, 25, 26, 30, 31, 35, 36, 38, 40, 43, 46, 53, 55, 65, 66, 67}
)
UNRESOLVED_STEREO_ROWS = frozenset({42, 61})
ALLOWED_PARENT_SALT_ROWS = frozenset({17, 41})
AMBIGUOUS_DERIVATIVE_ROWS = frozenset({14})
_MUTATION_TOKEN = re.compile(r"\([A-Z][0-9]+[A-Z](?:,[A-Z][0-9]+[A-Z])?\)")


@dataclass(frozen=True)
class ProjectionBudget:
    common_width: int = 32
    control_rank_2d: int = 8
    augmented_rank_2d: int = 4
    augmented_rank_3d: int = 4

    @property
    def control_interaction_parameters(self) -> int:
        return 2 * self.common_width * self.control_rank_2d

    @property
    def augmented_interaction_parameters(self) -> int:
        return 2 * self.common_width * (
            self.augmented_rank_2d + self.augmented_rank_3d
        )

    def validate(self) -> None:
        if self.common_width <= 0:
            raise ValueError("common projection width must be positive")
        if self.control_interaction_parameters != self.augmented_interaction_parameters:
            raise ValueError("interaction comparisons are not parameter matched")


def inchikey_connectivity(inchikey: str) -> str:
    if len(inchikey) < 14:
        raise ValueError("invalid InChIKey")
    return inchikey[:14]


def select_canonical_kinase_domain(
    assay_label: str, domain_features: Sequence[Mapping[str, Any]]
) -> tuple[str, Mapping[str, Any] | None]:
    """Select a UniProt core kinase domain without guessing assay constructs."""

    core = [
        feature
        for feature in domain_features
        if str(feature.get("description", ""))
        in {"Protein kinase", "Protein kinase 1", "Protein kinase 2"}
    ]
    if "JH1domain-catalytic" in assay_label:
        matches = [
            feature for feature in core if feature.get("description") == "Protein kinase 2"
        ]
        return (
            ("resolved_explicit_domain", matches[0])
            if len(matches) == 1
            else ("unresolved", None)
        )
    if "JH2domain-pseudokinase" in assay_label:
        matches = [
            feature for feature in core if feature.get("description") == "Protein kinase 1"
        ]
        return (
            ("resolved_explicit_domain", matches[0])
            if len(matches) == 1
            else ("unresolved", None)
        )
    if "Kin.Dom.1-N-terminal" in assay_label:
        matches = [
            feature for feature in core if feature.get("description") == "Protein kinase 1"
        ]
        return (
            ("resolved_explicit_domain", matches[0])
            if len(matches) == 1
            else ("unresolved", None)
        )
    if "Kin.Dom.2-C-terminal" in assay_label:
        matches = [
            feature for feature in core if feature.get("description") == "Protein kinase 2"
        ]
        return (
            ("resolved_explicit_domain", matches[0])
            if len(matches) == 1
            else ("unresolved", None)
        )
    if len(core) == 1:
        return "resolved_single_domain", core[0]
    return "unresolved", None


def receptor_exclusion_reasons(
    *,
    assay_label: str,
    mutant_flag: str,
    kinase_group: str,
    klifs_match_count: int,
    pocket_length: int,
    domain_status: str,
    alphafold_available: bool,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if mutant_flag != "NO" or _MUTATION_TOKEN.search(assay_label):
        reasons.append("mutant_or_mutation_label")
    if kinase_group not in STANDARD_EPK_GROUPS:
        reasons.append("outside_standard_epk_estimand")
    if "phosphorylated" in assay_label.lower():
        reasons.append("explicit_phosphorylation_state")
    if "cyclin" in assay_label.lower():
        reasons.append("explicit_binding_partner_state")
    if klifs_match_count != 1:
        reasons.append("klifs_mapping_not_unique")
    if pocket_length != 85:
        reasons.append("incomplete_fixed_klifs_pocket")
    if domain_status == "unresolved":
        reasons.append("canonical_domain_unresolved")
    if not alphafold_available:
        reasons.append("predicted_reference_view_unavailable")
    return tuple(reasons)


def connected_components(nodes: Iterable[str], edges: Iterable[tuple[str, str]]) -> list[list[str]]:
    parents = {node: node for node in nodes}

    def find(node: str) -> str:
        while parents[node] != node:
            parents[node] = parents[parents[node]]
            node = parents[node]
        return node

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    for left, right in edges:
        union(left, right)
    groups: dict[str, list[str]] = {}
    for node in parents:
        groups.setdefault(find(node), []).append(node)
    return sorted((sorted(group) for group in groups.values()), key=lambda group: group[0])


def ordered_identity_fraction(query: str, reference: str) -> float:
    """Return LCS(query, reference) / len(query) for pocket-to-domain auditing."""

    if not query:
        return 0.0
    previous = [0] * (len(reference) + 1)
    for query_residue in query:
        current = [0]
        for index, reference_residue in enumerate(reference, start=1):
            if query_residue == reference_residue:
                current.append(previous[index - 1] + 1)
            else:
                current.append(max(previous[index], current[-1]))
        previous = current
    return previous[-1] / len(query)
