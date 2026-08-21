import numpy as np

from iscore3.gate03.evaluation import (
    Dataset,
    Fold,
    Plan,
    _fixed_projection,
    _knn_observation,
    _predict_plan,
    pairwise_concordance,
)


def test_fixed_projection_is_deterministic_and_label_independent():
    value = np.arange(30, dtype=float).reshape(6, 5)
    assert np.array_equal(_fixed_projection(value, 3, 7), _fixed_projection(value, 3, 7))
    assert not np.array_equal(_fixed_projection(value, 3, 7), _fixed_projection(value, 3, 8))


def test_pairwise_concordance_handles_reverse_and_ties():
    y = np.asarray([1.0, 2.0, 3.0])
    assert pairwise_concordance(y, y) == 1.0
    assert pairwise_concordance(y, -y) == 0.0
    assert pairwise_concordance(y, np.zeros(3)) == 0.5


def test_observation_knn_maps_fold_local_responses_to_global_rows():
    similarity = np.eye(5)
    similarity[4, 3] = similarity[3, 4] = 0.9
    similarity[4, 1] = similarity[1, 4] = 0.8
    fitting = np.asarray([1, 3])
    prediction = _knn_observation(
        similarity,
        np.asarray([10.0, 30.0]),
        fitting,
        np.asarray([4]),
        1,
    )
    assert prediction.tolist() == [30.0]


def test_complete_case_scaffold_exclusion_returns_empty_evaluation_contract():
    dataset = Dataset(
        rows=tuple({} for _ in range(8)),
        ids=np.asarray([str(index) for index in range(8)]),
        series=np.asarray(["series"] * 8),
        components=np.asarray(["component"] * 8),
        scaffolds=np.asarray(["train"] * 7 + ["test"]),
        scaffold_eligible=np.ones(8, dtype=bool),
        y=np.arange(8, dtype=float),
        features={},
        available={"gmolai": np.ones(8, dtype=bool)},
        similarities={},
        series_order=np.asarray(["series"]),
        series_similarities={},
        derangements={},
    )
    fold = Fold(
        "scaffold",
        "one-row-test",
        "S1",
        np.arange(7),
        np.asarray([7]),
        "series",
        "component",
        "test",
    )
    prediction, _, _, evaluation, _ = _predict_plan(
        dataset,
        Plan("complete-case", "ridge", ("gmolai",), complete_case="gmolai"),
        fold,
        {},
    )
    assert prediction.size == 0
    assert evaluation.size == 0
