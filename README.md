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
src/unicodes/          the library: one subpackage per topic
    io/                loading .mat files and other data
tests/                 pytest tests, one file per subpackage
matlab_source/         the original .m files the ports are checked against
pyproject.toml         package name, version and dependencies
```

## How a MATLAB file becomes library code

| MATLAB | Python here |
|---|---|
| `function y = f(x)` in `f.m` | `def f(x):` in a topic module, e.g. `src/unicodes/aero/atmosphere.py` |
| script `.m` with hard-coded inputs | a function taking those inputs as arguments, plus an example under `examples/` |
| `struct` | `dict` or a `@dataclass` |
| 1-based indexing `x(1)`, `x(end)` | 0-based `x[0]`, `x[-1]` |
| `A*B` (matrix), `A.*B` (elementwise) | `A @ B`, `A * B` |
| `ode45`, `fzero`, `fsolve` | `scipy.integrate.solve_ivp`, `scipy.optimize.brentq`, `scipy.optimize.fsolve` |
| `plot`, `figure` | `matplotlib.pyplot` |

To make a new function importable, add it to its subpackage's
`__init__.py` (see `src/unicodes/io/__init__.py`).
