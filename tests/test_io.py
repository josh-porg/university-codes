import numpy as np
from scipy.io import savemat

from unicodes.io import load_mat


def test_load_mat_roundtrip(tmp_path):
    path = tmp_path / "data.mat"
    savemat(path, {"x": np.arange(5.0), "s": {"alt": 3.0, "name": "flight"}})

    data = load_mat(path)

    np.testing.assert_array_equal(data["x"], np.arange(5.0))
    assert data["s"]["alt"] == 3.0
    assert data["s"]["name"] == "flight"
