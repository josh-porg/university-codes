import numpy as np
import pytest

from unicodes.cfd import StructuredMesh2D, euler, residual, solve_steady

GAMMA = 1.4


def freestream(M, ny=1):
    return euler.conserved(1.0, [M * np.sqrt(GAMMA), 0.0], 1.0, GAMMA)  # rho = p = 1


def test_primitive_roundtrip():
    Q = euler.conserved(1.2, [100.0, -20.0], 9e4)
    s = euler.primitives(Q)
    assert s.rho == pytest.approx(1.2)
    np.testing.assert_allclose(s.velocity, [100.0, -20.0])
    assert s.p == pytest.approx(9e4)


def test_rusanov_is_consistent():
    Q = euler.conserved(1.0, [2.0, 0.5], 1.0)
    n = np.array([0.6, 0.8])
    np.testing.assert_allclose(euler.rusanov_flux(Q, Q, n), euler.normal_flux(Q, n))


def ramp_mesh(ni=60, nj=40, angle_deg=10.0):
    x = np.linspace(-0.5, 1.5, ni + 1)
    wall = np.where(x > 0, x * np.tan(np.radians(angle_deg)), 0.0)
    eta = np.linspace(0, 1, nj + 1)
    X = np.repeat(x[:, None], nj + 1, axis=1)
    Y = wall[:, None] + eta[None, :] * (1.5 - wall[:, None])
    return StructuredMesh2D(np.stack([X, Y], axis=-1))


def test_mesh_metrics():
    mesh = ramp_mesh()
    assert mesh.volume.sum() == pytest.approx(2.0 * 1.5 - 0.5 * 1.5**2 * np.tan(np.radians(10)), rel=1e-12)
    np.testing.assert_allclose(np.linalg.norm(mesh.i_faces[..., 1:], axis=-1), 1)


def test_freestream_is_preserved_on_skewed_mesh():
    mesh = ramp_mesh()
    Q = np.broadcast_to(freestream(2.0), mesh.shape + (4,)).copy()
    bcs = {"i_min": ("inlet", freestream(2.0)), "i_max": "exit", "j_min": ("inlet", freestream(2.0)), "j_max": "exit"}
    np.testing.assert_allclose(residual(mesh, Q, bcs, order=2), 0, atol=1e-10)


@pytest.mark.parametrize("order, rel", [(1, 0.03), (2, 0.005)])
def test_oblique_shock_pressure_ratio(order, rel):
    # Mach 2 over a 10 degree ramp: theory gives p2/p1 = 1.7066
    mesh = ramp_mesh()
    Q_in = freestream(2.0)
    Q0 = np.broadcast_to(Q_in, mesh.shape + (4,)).copy()
    bcs = {"i_min": ("inlet", Q_in), "i_max": "exit", "j_min": "wall", "j_max": "exit"}
    result = solve_steady(mesh, Q0, bcs, order=order, cfl=0.5, tol=1e-5, max_iterations=6000)
    assert result.converged
    p = euler.primitives(result.Q).p
    wall_p = p[(mesh.centroid[:, 0, 0] > 0.6) & (mesh.centroid[:, 0, 0] < 1.0), 0]
    assert wall_p.mean() == pytest.approx(1.7066, rel=rel)


def box_mesh(nx, ny, L=1.0, H=1.0):
    x, y = np.linspace(0, L, nx + 1), np.linspace(0, H, ny + 1)
    return StructuredMesh2D(np.stack(np.meshgrid(x, y, indexing="ij"), axis=-1))


def test_viscous_terms_vanish_for_uniform_flow():
    from unicodes.cfd import Viscosity

    mesh = ramp_mesh(30, 20)
    Q_in = freestream(2.0)
    Q = np.broadcast_to(Q_in, mesh.shape + (4,)).copy()
    bcs = {"i_min": ("inlet", Q_in), "i_max": "exit", "j_min": ("inlet", Q_in), "j_max": "exit"}
    visc = Viscosity(mu_ref=0.01)
    np.testing.assert_allclose(residual(mesh, Q, bcs, order=2, viscosity=visc), 0, atol=1e-10)


def test_shear_wave_decays_at_the_viscous_rate():
    # u = A cos(pi y / H) between slip walls decays as exp(-nu (pi / H)^2 t)
    from unicodes.cfd import Viscosity, solve_unsteady

    mesh = box_mesh(4, 40)
    y = mesh.centroid[..., 1]
    A, nu, t_end = 0.05, 0.05, 1.0
    Q0 = euler.conserved(np.ones_like(y), np.stack([A * np.cos(np.pi * y), 0 * y], -1), np.ones_like(y) / GAMMA)
    bcs = {"i_min": "exit", "i_max": "exit", "j_min": "symmetry", "j_max": "symmetry"}
    visc = Viscosity(mu_ref=nu, T_ref=1 / GAMMA, omega=0.0)  # rho = 1, constant viscosity
    out = solve_unsteady(mesh, Q0, bcs, t_end, order=2, cfl=0.4, viscosity=visc)
    u = euler.primitives(out.Q).velocity[..., 0]
    amplitude = np.sum(u * np.cos(np.pi * y)) / np.sum(np.cos(np.pi * y) ** 2)
    assert out.time == pytest.approx(t_end)
    assert amplitude / A == pytest.approx(np.exp(-nu * np.pi**2 * t_end), rel=0.02)


def test_noslip_wall_ghost_reverses_velocity():
    from unicodes.cfd.solver2d import _ghost

    Q = euler.conserved(1.2, [3.0, -1.0], 2.0)[None, None]
    faces = np.array([[[1.0, 0.0, 1.0]]])
    g = euler.primitives(_ghost("noslip", Q, faces, GAMMA), GAMMA)
    np.testing.assert_allclose(g.velocity, [[[-3.0, 1.0]]])
    assert g.p.item() == pytest.approx(2.0) and g.rho.item() == pytest.approx(1.2)
    gT = euler.primitives(_ghost(("noslip", 1.0), Q, faces, GAMMA), GAMMA)
    T_cell = 2.0 / 1.2
    assert (gT.p / gT.rho).item() == pytest.approx(2 * 1.0 - T_cell)  # wall temperature is the face average
