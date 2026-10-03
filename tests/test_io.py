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


def _unit_square_mesh():
    # 3x2 nodes -> 2 quad cells
    xyz = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], dtype=float)
    conn = np.array([[0, 1, 4, 3], [1, 2, 5, 4]])
    return xyz, conn


def test_tecplot_roundtrip_cell_centred(tmp_path):
    from unicodes.io import read_block_data, read_grid, write_block_data

    xyz, conn = _unit_square_mesh()
    values = np.array([[1.5, 10.0], [2.5, 20.0]])
    path = tmp_path / "flow.dat"
    write_block_data(path, xyz, values, ["P", "T"], conn)

    text = path.read_text()
    assert "VARLOCATION=([3-4]=CELLCENTERED)" in text
    header = text[: text.index("DT=(")].count("\n") + 1  # lines up to and including DT=

    # Stored blocks: x, y (nodal), then P, T (cell-centred)
    out = read_block_data(path, header, variables=[2, 3], n_nodes=6, n_elements=2,
                          cell_centered=[2, 3])
    np.testing.assert_allclose(out, values)


def test_tecplot_read_grid(tmp_path):
    from unicodes.io import read_grid

    xyz, conn = _unit_square_mesh()
    path = tmp_path / "grid.dat"
    lines = ["header"] + [" ".join(map(str, xyz[:, j])) for j in range(2)]
    lines += [" ".join(str(n + 1) for n in row) for row in conn]
    path.write_text("\n".join(lines) + "\n")

    grid_xyz, grid_conn = read_grid(path, 1, ndim=2, n_nodes=6, n_elements=2)
    np.testing.assert_allclose(grid_xyz, xyz)
    np.testing.assert_array_equal(grid_conn, conn)


def test_node_to_cell_average():
    from unicodes.io import node_to_cell

    _, conn = _unit_square_mesh()
    np.testing.assert_allclose(node_to_cell(np.arange(6.0), conn), [2.0, 3.0])
