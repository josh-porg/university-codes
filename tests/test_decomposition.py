import numpy as np
import pytest

from unicodes.decomposition import (
    DMD,
    POD,
    data_to_snapshots,
    optimal_svht_coef,
    remove_mean,
    rsvd,
    snapshots_to_data,
)


def travelling_waves(nt=80, dt=0.05):
    x = np.linspace(-5, 5, 64)
    y = np.linspace(-3, 3, 32)
    t = np.arange(nt) * dt
    X, Y, T = np.meshgrid(x, y, t, indexing="ij")
    # Each oscillation needs two spatial shapes (cos and sin parts) for real data
    g = np.exp(-(X**2) - Y**2)
    f1 = g * np.cos(2.0 * T) + X * g * np.sin(2.0 * T)  # omega = ±2i
    h = np.tanh(X) / np.cosh(Y) * np.exp(-0.3 * T)
    f2 = h * (np.sin(5.0 * T) + Y * np.cos(5.0 * T))  # omega = -0.3 ± 5i
    return f1 + f2, dt


def test_snapshot_roundtrip():
    data = np.random.default_rng(0).normal(size=(4, 5, 6, 7))
    snaps = data_to_snapshots(data)
    assert snaps.shape == (120, 7)
    np.testing.assert_array_equal(snapshots_to_data(snaps, (4, 5, 6)), data)
    fluct, mean = remove_mean(snaps)
    np.testing.assert_allclose(fluct.mean(axis=1), 0, atol=1e-12)
    np.testing.assert_allclose(fluct + mean, snaps)


def test_dmd_recovers_frequencies_and_reconstructs():
    data, dt = travelling_waves()
    dmd = DMD(data, rank=4, dt=dt, truncation=None)
    expected = np.array([2j, -2j, -0.3 + 5j, -0.3 - 5j])
    for w in expected:
        assert np.min(np.abs(dmd.omega - w)) < 1e-6
    np.testing.assert_allclose(dmd.reconstruction, data, atol=1e-6)
    assert dmd.get_modes().shape == (64, 32, 4)


def test_dmd_background_separation():
    data, dt = travelling_waves()
    still = np.ones_like(data) * np.linspace(0, 1, data.shape[0])[:, None, None]
    dmd = DMD(data + still, rank=5, dt=dt, truncation=None, background_threshold=1e-3)
    np.testing.assert_allclose(dmd.background, still, atol=1e-6)
    np.testing.assert_allclose(dmd.foreground, data, atol=1e-6)


def test_pod_is_orthonormal_and_exact():
    data, _ = travelling_waves()
    pod = POD(data)
    k = 4
    Phi = pod.spatial_modes[:, :k]
    np.testing.assert_allclose(Phi.T @ Phi, np.eye(k), atol=1e-10)
    assert pod.energy_fraction[:4].sum() == pytest.approx(1.0)
    np.testing.assert_allclose(pod.reconstruct(k), data, atol=1e-8)


def test_rsvd_matches_svd_for_low_rank_matrix():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(500, 5)) @ rng.normal(size=(5, 60))
    U, S, Vh = rsvd(X, 5, rng=2)
    np.testing.assert_allclose(S, np.linalg.svd(X, compute_uv=False)[:5], rtol=1e-8)
    np.testing.assert_allclose(U * S @ Vh, X, atol=1e-8)


def test_optimal_svht_coef_known_values():
    # Square matrix: 4/sqrt(3) with known noise; ~2.858 with unknown noise (Gavish & Donoho)
    assert optimal_svht_coef(1.0, sigma_known=True) == pytest.approx(4 / np.sqrt(3))
    assert optimal_svht_coef(1.0) == pytest.approx(2.858, abs=1e-3)
