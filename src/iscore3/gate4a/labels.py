"""Censor-aware affinity observation contract for Gate-4A."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import math
import re
from typing import Any


_NUMBER_WITH_RELATION = re.compile(
    r"^\s*(?P<relation>>=|>|<=|<)?\s*"
    r"(?P<value>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)\s*$"
)


class ObservationKind(str, Enum):
    EXACT = "exact"
    RIGHT_CENSORED_KD = "right_censored_kd"
    MISSING = "missing"
    INVALID = "invalid"


@dataclass(frozen=True)
class AffinityObservation:
    kind: ObservationKind
    source_text: str
    kd_nm: float | None = None
    kd_lower_bound_nm: float | None = None
    pkd: float | None = None
    pkd_upper_bound: float | None = None
    error: str | None = None

    def active_at_nm(self, threshold_nm: float) -> bool | None:
        """Return activity at Kd <= threshold, or None when censoring is insufficient."""

        if threshold_nm <= 0 or not math.isfinite(threshold_nm):
            raise ValueError("activity threshold must be finite and positive")
        if self.kind is ObservationKind.EXACT:
            assert self.kd_nm is not None
            return self.kd_nm <= threshold_nm
        if self.kind is ObservationKind.RIGHT_CENSORED_KD:
            assert self.kd_lower_bound_nm is not None
            return False if self.kd_lower_bound_nm >= threshold_nm else None
        return None

    def to_record(self) -> dict[str, Any]:
        record = asdict(self)
        record["kind"] = self.kind.value
        return record


def pkd_from_nm(kd_nm: float) -> float:
    if kd_nm <= 0 or not math.isfinite(kd_nm):
        raise ValueError("Kd must be finite and positive")
    return 9.0 - math.log10(kd_nm)


def parse_kd_cell(
    raw_value: Any,
    *,
    blank_is_censored: bool,
    censor_limit_nm: float = 10_000.0,
) -> AffinityObservation:
    """Parse a source Kd cell without converting censoring into an exact value."""

    if censor_limit_nm <= 0 or not math.isfinite(censor_limit_nm):
        raise ValueError("censor_limit_nm must be finite and positive")
    if raw_value is None:
        text = ""
    elif isinstance(raw_value, bool):
        return AffinityObservation(
            kind=ObservationKind.INVALID,
            source_text=str(raw_value),
            error="boolean is not an affinity",
        )
    else:
        text = str(raw_value).strip()

    if not text:
        if blank_is_censored:
            return AffinityObservation(
                kind=ObservationKind.RIGHT_CENSORED_KD,
                source_text="",
                kd_lower_bound_nm=censor_limit_nm,
                pkd_upper_bound=pkd_from_nm(censor_limit_nm),
            )
        return AffinityObservation(kind=ObservationKind.MISSING, source_text="")

    match = _NUMBER_WITH_RELATION.fullmatch(text)
    if match is None:
        return AffinityObservation(
            kind=ObservationKind.INVALID,
            source_text=text,
            error="unrecognized Kd cell",
        )
    value = float(match.group("value"))
    if value <= 0 or not math.isfinite(value):
        return AffinityObservation(
            kind=ObservationKind.INVALID,
            source_text=text,
            error="Kd must be finite and positive",
        )
    relation = match.group("relation")
    if relation in {">", ">="}:
        return AffinityObservation(
            kind=ObservationKind.RIGHT_CENSORED_KD,
            source_text=text,
            kd_lower_bound_nm=value,
            pkd_upper_bound=pkd_from_nm(value),
        )
    if relation in {"<", "<="}:
        return AffinityObservation(
            kind=ObservationKind.INVALID,
            source_text=text,
            error="left-censored Kd is not yet supported by the Gate-4A contract",
        )
    return AffinityObservation(
        kind=ObservationKind.EXACT,
        source_text=text,
        kd_nm=value,
        pkd=pkd_from_nm(value),
    )
