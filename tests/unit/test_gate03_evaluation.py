import numpy as np

from iscore3.gate03.evaluation import (
    _fixed_projection,
    _knn_observation,
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
