"""DMD pipeline for the combustor LES snapshots (``load_snapshots*``, ``load_driver``, ``DMD``,
``DMD_driver``, ``DMD_function_driver_*``, ``DMD_combustor_driver``, ``export_snapshots``,
``export_driver``, ``saveModesForExport``, ``convert_mat_two_bin_*``,
``create_bin_from_mat_data``, ``undersample_data_for_combustor``).

Steps (each a sub-command; the Tecplot files are not in the repository)::

    python combustor_dmd.py load --grid grid.dat --pattern "test_file_ncons_{}.dat" --start 150000 --end 220000 --variable 0 -o P.npy
    python combustor_dmd.py dmd P.npy --dt 1e-7 -o P_dmd.npz [--svht] [--undersample 2]
    python combustor_dmd.py dmd combustorData150000_169999.mat -o comb_dmd.npz   (DMD_combustor)
    python combustor_dmd.py export P_dmd.npz --grid grid.dat --what modes --prefix dmd700000_P_modes
    python combustor_dmd.py bin P.npy -o P.bin

Variables in the files: 0 P, 1 U, 2 V, 3 T, 4 Y_CH4, 5 Y_O2, 6 Y_H2O, 7 Y_CO2.
Snapshot arrays are stored ``(n_cells, n_time)``. The MATLAB grid had
39065 nodes and 38523 quadrilateral cells (pass ``--nodes/--cells``).
"""

import argparse

import numpy as np

from unicodes.decomposition import DMD
from unicodes.io import read_block_data, read_grid, write_block_data

VARIABLES = ["P", "U", "V", "T", "Y_CH4", "Y_O2", "Y_H2O", "Y_CO2"]

p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
sub = p.add_subparsers(dest="cmd", required=True)
for name in ("load", "export"):
    s = sub.add_parser(name)
    s.add_argument("--grid", required=True)
    s.add_argument("--nodes", type=int, default=39065)
    s.add_argument("--cells", type=int, default=38523)
    s.add_argument("--grid-header", type=int, default=9)
s = sub.choices["load"]
s.add_argument("--pattern", required=True, help="file name with {} for the time index")
s.add_argument("--start", type=int, required=True)
s.add_argument("--end", type=int, required=True)
s.add_argument("--step", type=int, default=1)
s.add_argument("--variable", type=int, default=0)
s.add_argument("--header", type=int, default=15)
s.add_argument("-o", "--output", required=True)
s = sub.choices["export"]
s.add_argument("results")
s.add_argument("--what", choices=["modes", "reconstruction"], default="modes")
s.add_argument("--prefix", required=True)
s.add_argument("--count", type=int, default=None, help="number of columns to write")
s = sub.add_parser("dmd")
s.add_argument("snapshots")
s.add_argument("--dt", type=float, default=1e-7)
s.add_argument("--svht", action="store_true", help="Gavish-Donoho truncation instead of 0.5 % of the singular-value sum")
s.add_argument("--undersample", type=int, default=1)
s.add_argument("-o", "--output", required=True)
s = sub.add_parser("bin")
s.add_argument("snapshots")
s.add_argument("-o", "--output", required=True)
a = p.parse_args()

if a.cmd == "load":
    xyz, conn = read_grid(a.grid, a.grid_header, 2, a.nodes, a.cells)
    cols = []
    for i in range(a.start, a.end + 1, a.step):
        data = read_block_data(a.pattern.format(i), a.header, [a.variable], a.nodes, a.cells,
                               cell_centered=range(8), connectivity=conn)
        cols.append(np.ravel(data))
    np.save(a.output, np.column_stack(cols))
    print(f"saved {len(cols)} snapshots of {VARIABLES[a.variable]} to {a.output}")
elif a.cmd == "dmd":
    if a.snapshots.endswith(".mat"):  # MATLAB files store p as (time, cells)
        from unicodes.io import load_mat

        X = np.asarray(load_mat(a.snapshots)["p"], float).T[:, :: a.undersample]
    else:
        X = np.load(a.snapshots)[:, :: a.undersample]
    dmd = DMD(X, dt=a.dt * a.undersample, truncation="svht" if a.svht else 0.005)
    np.savez(a.output, modes=dmd.modes, eigenvalues=dmd.eigenvalues, omega=dmd.omega, amplitudes=dmd.amplitudes,
             time_dynamics=dmd.time_dynamics(), reconstruction=dmd.reconstruct().real,
             power=dmd.growth_adjusted_power())
    print(f"rank {dmd.rank}; strongest frequencies (Hz):",
          np.round(np.sort(np.abs(dmd.frequencies[np.argsort(dmd.growth_adjusted_power())[-5:]])), 1))
elif a.cmd == "export":
    xyz, conn = read_grid(a.grid, a.grid_header, 2, a.nodes, a.cells)
    data = np.load(a.results)[a.what]
    data = data.real if np.iscomplexobj(data) else data
    for k in range(data.shape[1] if a.count is None else a.count):
        write_block_data(f"{a.prefix}_flowfield_{k + 1}.dat", xyz, data[:, k:k + 1], [a.what], conn)
    print("written")
else:
    np.load(a.snapshots).T.astype("<f8").tofile(a.output)  # MATLAB fwrite of p (time x space), column-major
    print(f"wrote {a.output}")
