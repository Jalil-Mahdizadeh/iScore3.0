import sys

import torch

from iscore3.protein.esm_if1_adapter import _install_native_scatter_compatibility


def test_native_scatter_add_matches_explicit_sum():
    sys.modules.pop("torch_scatter", None)
    _install_native_scatter_compatibility()
    from torch_scatter import scatter_add

    source = torch.tensor([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    index = torch.tensor([1, 0, 1])
    observed = scatter_add(source, index, dim=0, dim_size=3)
    expected = torch.tensor([[3.0, 4.0], [6.0, 8.0], [0.0, 0.0]])
    assert torch.equal(observed, expected)
