import math

import pytest

from iscore3.gate4a.labels import ObservationKind, parse_kd_cell, pkd_from_nm


def test_blank_davis_cell_is_censored_not_exact() -> None:
    observation = parse_kd_cell("", blank_is_censored=True, censor_limit_nm=10_000.0)
    assert observation.kind is ObservationKind.RIGHT_CENSORED_KD
    assert observation.kd_nm is None
    assert observation.kd_lower_bound_nm == 10_000.0
    assert observation.pkd_upper_bound == 5.0


def test_explicit_numeric_10000_is_not_silently_relabelled() -> None:
    observation = parse_kd_cell("10000", blank_is_censored=True)
    assert observation.kind is ObservationKind.EXACT
    assert observation.kd_nm == 10_000.0
    assert observation.pkd == 5.0


def test_censored_observation_supports_only_justified_activity_thresholds() -> None:
    observation = parse_kd_cell("", blank_is_censored=True)
    assert observation.active_at_nm(100.0) is False
    assert observation.active_at_nm(1_000.0) is False
    assert observation.active_at_nm(20_000.0) is None


def test_invalid_or_unsupported_cells_remain_visible() -> None:
    assert parse_kd_cell("not measured", blank_is_censored=True).kind is ObservationKind.INVALID
    assert parse_kd_cell("<3", blank_is_censored=True).kind is ObservationKind.INVALID


def test_pkd_conversion_uses_nanomolar_units() -> None:
    assert pkd_from_nm(1.0) == 9.0
    assert pkd_from_nm(100.0) == 7.0
    assert math.isclose(pkd_from_nm(3.0), 9.0 - math.log10(3.0))
    with pytest.raises(ValueError):
        pkd_from_nm(0.0)
