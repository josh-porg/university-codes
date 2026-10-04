import numpy as np
import pytest
from scipy.fft import idct

from unicodes.decomposition import (cosamp, dmd_eigenvalues, era, fast_ica, kernel_dmd, mrdmd, optimal_amplitudes,
                                    time_delay_stack)

A_TRUE = np.array([[1, 1], [-1, 2]]) / np.sqrt(3)


def lds(m=100):
    X = np.zeros((2, m))
    X[:, 0] = [0.5, 1]
    for k in range(1, m):
        X[:, k] = A_TRUE @ X[:, k - 1]
    return X


@pytest.mark.parametrize("flavor", ["dmd", "fbdmd", "tlsdmd"])
def test_noise_free_eigenvalues(flavor):
    lam = np.sort_complex(dmd_eigenvalues(lds(), flavor))
    np.testing.assert_allclose(lam, np.sort_complex(np.linalg.eigvals(A_TRUE)), atol=1e-8)


def test_fbdmd_reduces_noise_bias():
    rng = np.random.default_rng(1)
    X = lds()
    true = np.linalg.eigvals(A_TRUE)
    err = {f: [] for f in ("dmd", "fbdmd")}
    for _ in range(200):
        Y = X + 0.2 * rng.standard_normal(X.shape)
        for f in err:
            lam = dmd_eigenvalues(Y, f)
            err[f].append(np.min(np.abs(lam[:, None] - true[None, :]), axis=0).mean())
    assert np.mean(err["fbdmd"]) < np.mean(err["dmd"])


def test_era_recovers_two_decay_rates():
    dt = 0.01
    t = np.arange(0, 100 + dt / 2, dt)
    x = np.exp(-0.01 * t) + np.exp(-10 * t)
    A, B, C, D, hsv = era(x[None, None, :], 100, 100, 2)
    rates = np.sort(np.log(np.linalg.eigvals(A)).real / dt)
    np.testing.assert_allclose(rates, [-10, -0.01], rtol=1e-6)
    y = [float((C @ np.linalg.matrix_power(A, k) @ B).item()) for k in (0, 50)]
    np.testing.assert_allclose(y, x[[1, 51]], rtol=1e-8)


def test_optimal_amplitudes_fit_noise_free_data():
    X = lds(30)
    U, S, Vh = np.linalg.svd(X[:, :-1], full_matrices=False)
    lam, W = np.linalg.eig(U.T @ X[:, 1:] @ Vh.T / S)
    Phi = X[:, 1:] @ Vh.T / S @ W
    b = optimal_amplitudes(lam, W, S, Vh.T)
    recon = U @ W @ (b[:, None] * lam[:, None] ** np.arange(29))  # projected modes U W (Jovanovic et al.)
    np.testing.assert_allclose(recon, X[:, :-1], atol=1e-10)
    np.testing.assert_allclose(Phi @ (b[:, None] * lam[:, None] ** np.arange(29)), X[:, 1:], atol=1e-10)  # exact modes: one step on


def test_time_delay_stack_shape():
    H = time_delay_stack(np.arange(10.0), 3)
    assert H.shape == (3, 8) and H[2, 0] == 2 and H[0, -1] == 7


def test_cosamp_recovers_sparse_dct_signal():
    rng = np.random.default_rng(0)
    n, p = 1024, 128
    s_true = np.zeros(n)
    s_true[[40, 300]] = [1.0, -0.7]
    Psi = idct(np.eye(n), norm="ortho", axis=0)
    perm = rng.choice(n, p, replace=False)
    s = cosamp(Psi[perm], Psi[perm] @ s_true, 4, 1e-10, 50)
    np.testing.assert_allclose(s, s_true, atol=1e-8)


def test_fast_ica_separates_two_sources():
    t = np.linspace(0, 10, 4000)
    S = np.vstack([np.sin(2 * t), np.sign(np.sin(3 * t))])
    X = np.array([[1.0, 0.5], [0.4, 1.0]]) @ S
    sources, mixing = fast_ica(X, 2, rng=0)
    corr = np.abs(np.corrcoef(np.vstack([sources, S]))[:2, 2:])
    assert np.all(corr.max(axis=1) > 0.99)


def test_kernel_dmd_linear_kernel_matches_dmd():
    X = lds(20)
    lam_k = kernel_dmd(X[:, :-1], X[:, 1:], lambda a, b: a @ b)
    np.testing.assert_allclose(np.sort_complex(lam_k), np.sort_complex(np.linalg.eigvals(A_TRUE)), atol=1e-8)


def test_mrdmd_finds_slow_mode_at_top_level():
    dt = 0.01
    t = np.arange(0, 10, dt)
    x = np.outer([1.0, 0.5, -0.3], np.cos(2 * np.pi * 0.15 * t)) + np.outer([0.2, 1.0, 0.1], np.cos(2 * np.pi * 5.0 * t))
    tree = mrdmd(x, dt, rank=6, max_cycles=2, levels=3)
    top = tree.level(0)[0]
    assert top.rho == pytest.approx(0.2)
    assert np.any(np.isclose(np.abs(top.omega.imag), 0.15, atol=0.01))
    amp, low_f = tree.amplitude_map()
    assert amp.shape == (3, 4) and low_f[1] == pytest.approx(0.2)
