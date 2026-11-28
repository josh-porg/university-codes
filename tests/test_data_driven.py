import numpy as np
from scipy.integrate import solve_ivp

from unicodes.decomposition import (
    DMD,
    SINDy,
    dmd_spectrum,
    dmdc,
    exact_dmd,
    fft_amplitude,
    hankel,
    library,
    library_names,
    lorenz,
    stls,
    stridge,
)
from unicodes.numerics import tv_derivative


def test_library_shapes_and_names():
    X = np.random.default_rng(0).random((5, 3))
    terms = ("const", "linear", "poly2", "poly3", "sin")
    Theta = library(X, terms, n_harmonics=2)
    names = library_names(["x", "y", "z"], terms, n_harmonics=2)
    assert Theta.shape == (5, 1 + 3 + 6 + 10 + 12) == (5, len(names))
    assert names[:5] == ["1", "x", "y", "z", "xx"]
    assert np.allclose(Theta[:, names.index("xyz")], X.prod(axis=1))


def test_sindy_recovers_lorenz():
    t = np.arange(0, 10, 0.002)
    X = solve_ivp(lorenz, (0, 10), [0, 1, 20], t_eval=t, rtol=1e-10, atol=1e-10).y.T
    dX = np.array([lorenz(0, x) for x in X])
    m = SINDy(("const", "linear", "poly2"), threshold=0.5).fit(X, dX)
    names = library_names(["x", "y", "z"], m.terms)
    expected = np.zeros_like(m.coefficients)
    expected[names.index("x"), 0], expected[names.index("y"), 0] = -10, 10
    expected[names.index("x"), 1], expected[names.index("y"), 1], expected[names.index("xz"), 1] = 28, -1, -1
    expected[names.index("z"), 2], expected[names.index("xy"), 2] = -8 / 3, 1
    assert np.allclose(m.coefficients, expected, atol=1e-6)


def test_stls_equals_stridge_without_ridge():
    rng = np.random.default_rng(1)
    Theta = rng.standard_normal((200, 6))
    Xi = np.array([[1.0, 0], [0, 0], [0, -2], [0.05, 0], [0, 0], [3, 1]])
    dX = Theta @ Xi + 1e-3 * rng.standard_normal((200, 2))
    assert np.allclose(stls(Theta, dX, 0.2), stridge(Theta, dX, 0.2), atol=1e-10)


def test_dmdc_recovers_linear_system():
    rng = np.random.default_rng(2)
    A = np.array([[0.9, 0.1, 0], [-0.1, 0.9, 0], [0, 0, 0.5]])
    B = np.array([[1.0], [0], [0.5]])
    X = np.zeros((3, 60))
    U = rng.standard_normal((1, 60))
    X[:, 0] = [1, 0, 0]
    for k in range(59):
        X[:, k + 1] = A @ X[:, k] + B @ U[:, k]
    r = dmdc(X, U, tol=1e-10)
    assert np.allclose(r.A, A, atol=1e-8) and np.allclose(r.B, B, atol=1e-8)


def test_spectra_find_tones():
    dt = 1e-3
    t = np.arange(0, 1, dt)
    x = np.sin(2 * np.pi * 50 * t) + 0.5 * np.sin(2 * np.pi * 120 * t)
    f, P = fft_amplitude(x, dt)
    assert set(np.round(f[np.argsort(P)[-2:]])) == {50, 120}
    fd, Pd = dmd_spectrum(x, dt, delays=50, rank=4)
    assert np.allclose(np.sort(np.abs(fd)), [50, 50, 120, 120], atol=1e-6)
    assert hankel(np.arange(5), 2).tolist() == [[0, 1, 2, 3], [1, 2, 3, 4]]


def test_tv_derivative_small_and_large():
    x = np.linspace(0, 2 * np.pi, 300)
    dx = x[1] - x[0]
    y = np.sin(x)
    u_small = tv_derivative(y, 10, 1e-3, dx)
    edges = np.r_[x - dx / 2, x[-1] + dx / 2]
    assert np.sqrt(np.mean((u_small[3:-3] - np.cos(edges[3:-3])) ** 2)) < 0.05
    u_large = tv_derivative(y, 10, 1e-2, dx, scale="large")
    assert np.sqrt(np.mean((u_large[3:-3] - np.cos(x[3:-3])) ** 2)) < 0.05


def test_exact_dmd_and_growth_power():
    A = np.array([[0.95, 0.2], [-0.2, 0.95]])
    X = np.zeros((2, 40))
    X[:, 0] = [1, 0]
    for k in range(39):
        X[:, k + 1] = A @ X[:, k]
    lam, modes, _ = exact_dmd(X[:, :-1], X[:, 1:])
    assert np.allclose(np.sort_complex(lam), np.sort_complex(np.linalg.eigvals(A)))
    d = DMD(X, truncation=None)
    assert d.growth_adjusted_power().shape == (2,)
    assert np.allclose(np.abs(d.frequencies), np.angle(np.linalg.eigvals(A)[0]) / (2 * np.pi))
