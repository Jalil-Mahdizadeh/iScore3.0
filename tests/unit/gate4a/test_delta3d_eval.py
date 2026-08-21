import numpy as np

from iscore3.gate4a.delta3d_eval import (
    BranchProjector,
    aggregate_label_matrix,
    deterministic_group_folds,
    fit_linear_tobit,
    per_ligand_nll,
)


def test_group_folds_are_deterministic_balanced_and_keep_groups_whole() -> None:
    groups = ["A", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    first = deterministic_group_folds(groups, 3)
    second = deterministic_group_folds(groups, 3)
    np.testing.assert_array_equal(first, second)
    assert first[0] == first[1]
    loads = np.bincount(first, minlength=3)
    assert int(loads.max() - loads.min()) <= 1


def test_projector_is_training_determined_and_has_frozen_width() -> None:
    rng = np.random.default_rng(91)
    training = rng.normal(size=(45, 60))
    projector = BranchProjector.fit(training, width=32)
    assert projector.transform(training).shape == (45, 32)
    assert projector.transform(rng.normal(size=(7, 60))).shape == (7, 32)
    pivots = np.argmax(np.abs(projector.components), axis=1)
    assert np.all(projector.components[np.arange(32), pivots] >= 0.0)


def test_tobit_uses_censoring_without_substituting_exact_limits() -> None:
    x = np.asarray([[-2.0], [-1.0], [0.0], [1.0], [2.0]])
    exact = np.asarray(
        [[0.0, 0.0], [0.0, 0.0], [5.5, 0.0], [6.5, 7.0], [8.0, 8.5]]
    )
    mask = exact > 0.0
    labels = aggregate_label_matrix(exact, mask)
    model = fit_linear_tobit(x, labels, alpha=0.001)
    prediction = model.predict(x)
    assert model.sigma > 0.0
    assert prediction[-1] > prediction[0]
    loss, counts = per_ligand_nll(prediction, model.sigma, labels)
    assert np.isfinite(loss).all()
    assert int(counts.sum()) == 10
