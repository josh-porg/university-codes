# unicodes

Python ports of my university MATLAB codes, packaged as a library so any
future project can `import` them instead of copy-pasting.

## Using it in another project

Install straight from GitHub (run this inside the other project's virtual
environment):

```bash
pip install "git+https://github.com/josh-porg/university-codes.git"
```

Then in that project's Python code:

```python
from unicodes.io import load_mat

flight = load_mat("Maiden Voyage.mat")
```

To pin a specific version, so later changes here can't break the other
project, add a tag to the URL: `...university-codes.git@v0.1.0`. To list
the library as a dependency, put the same line in that project's
`requirements.txt`:

```
unicodes @ git+https://github.com/josh-porg/university-codes.git@v0.1.0
```

## Working on the library itself

```bash
git clone https://github.com/josh-porg/university-codes.git
cd university-codes
python -m venv .venv && source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -e ".[dev,plot,mat73]"   # -e = "editable": your edits take effect immediately
pytest
```

## Layout

```
src/unicodes/              the library: one module or subpackage per topic
    aero/                  wing, airfoil, NACA, thin airfoil, boundary layer, pressure taps,
                           performance/constraints, V-n loads, soaring, stability, propeller,
                           sizing, PAH camber surfaces, XFOIL runner
    atmosphere.py          1976 standard atmosphere, viscosity, speed of sound
    cfd/                   Euler fluxes, 2-D structured solver, quasi-1D nozzle, advection
                           schemes, grids, flux-reconstruction operators
    composites/            classical lamination theory (Ply, Laminate, micromechanics)
    controls.py            linear flight-dynamics models, modes, doublets, linear MPC
    decomposition/         DMD, DMDc, POD, SINDy, randomized SVD, SVHT, DMD/FFT spectra
    flight_dynamics/       6-DOF model with passive aero-compliant flaps, Dryden turbulence,
                           ISO 2631 ride quality, kinematics/quaternions
    fluids.py              pipe sections, friction factors
    gasdynamics.py         isentropic flow, normal/oblique shocks, expansions, inlets
    io/                    MATLAB .mat, Tecplot ASCII, LaTeX tables
    materials.py           isotropic material records
    neural.py              small sigmoid network trained by backpropagation
    numerics.py            Jacobians, TV-regularised derivatives
    ode.py                 fixed-step Euler, Heun, SSP-RK3, RK4
    optimize/              GA, PSO, cross-entropy, constrained steepest descent, simulated
                           annealing, line searches, benchmark functions
    orbital.py             two-body propagation, Kepler, orbital elements, transfers
    propulsion.py          rocket-assisted projectile record
    propulsion_cycles.py   gas-turbine cycle analysis: turbojet/afterburner, turbofans (separate and
                           mixed exhaust), turboprop, nozzles, compressor stage
    remote_sensing.py      Planck/Wien, diffraction-limited optics, SAR radar equation
    structures/            wing-box stresses and buckling, beam finite elements
    thermo/                combustion and flame temperature, equations of state, relation solver
    units.py               unit conversion constants
examples/<course>/         the homework, labs, exams and projects as runnable scripts
docs/CONVERSION.md         every MATLAB file and where its port lives (or why it was not ported)
tests/                     pytest tests
matlab_source/             the original .m/.mlx files the ports are checked against
pyproject.toml             package name, version and dependencies
```

Run an example from its folder, e.g. `python examples/ae746_cfd/project2_nozzle.py --case 2`.
Scripts that need data files the repository does not hold (spreadsheets,
`.mat` logs, simulation output) take the path as a command-line argument;
each script's docstring says which file it expects. A few examples use
optional packages: `sympy` (symbolic homework), `pandas` + `openpyxl`
(`.xlsx` input), `h5py` (MATLAB v7.3 / HDF5 files), `imageio` (video).

## How a MATLAB file becomes library code

| MATLAB | Python here |
|---|---|
| `function y = f(x)` in `f.m` | `def f(x):` in a topic module, e.g. `src/unicodes/atmosphere.py` |
| script `.m` with hard-coded inputs | a function taking those inputs as arguments, plus an example under `examples/` |
| `struct` | `dict` or a `@dataclass` |
| 1-based indexing `x(1)`, `x(end)` | 0-based `x[0]`, `x[-1]` |
| `A*B` (matrix), `A.*B` (elementwise) | `A @ B`, `A * B` |
| `ode45`, `fzero`, `fsolve` | `scipy.integrate.solve_ivp`, `scipy.optimize.brentq`, `scipy.optimize.fsolve` |
| `plot`, `figure` | `matplotlib.pyplot` |

To make a new function importable, add it to its subpackage's
`__init__.py` (see `src/unicodes/io/__init__.py`).
