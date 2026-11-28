"""Linear flight-dynamics models, mode analysis and a small linear MPC (AE 551).

Replaces the Control System and MPC toolbox calls in the AE 551 homework and
final project (``ss``, ``ss2tf``, ``lsim``, ``mpc``) with scipy.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal
from scipy.linalg import expm
from scipy.optimize import lsq_linear

G_FT = 32.2


def wrap_to_pi(angle):
    """Wrap angles to [-pi, pi) (``MyPiWrap``/``wrapToPi``)."""
    return (np.asarray(angle) + np.pi) % (2 * np.pi) - np.pi


@dataclass
class Mode:
    eigenvalue: complex
    natural_frequency: float  # rad/s
    damping_ratio: float
    damped_frequency: float  # rad/s

    @property
    def time_constant(self) -> float:
        return -1 / self.eigenvalue.real if self.eigenvalue.real else np.inf


def modes(A):
    """Eigenvalues of ``A`` with natural frequency, damping ratio and damped frequency (AE 551 HW 6)."""
    out = []
    for lam in np.linalg.eigvals(A):
        wn = abs(lam)
        out.append(Mode(lam, wn, -lam.real / wn if wn else 0.0, abs(lam.imag)))
    return sorted(out, key=lambda m: -m.natural_frequency)


def longitudinal_state_space(U1, theta1, X_u, X_alpha, Z_u, Z_alpha, M_u, M_alpha, M_q, Z_de, M_de,
                             X_q=0.0, Z_q=0.0, Z_alpha_dot=0.0, M_alpha_dot=0.0, X_de=0.0, g=G_FT):
    """Roskam dimensional longitudinal model, states ``[u, alpha, q, theta]``, input elevator (HW 6).

    The alpha-dot derivatives enter through ``C' x_dot = A' x + B' u`` and are
    eliminated as in the MATLAB (``A = C' \\ A'``).
    """
    C = np.array([[1, 0, 0, 0], [0, U1 - Z_alpha_dot, 0, 0], [0, -M_alpha_dot, 1, 0], [0, 0, 0, 1]], dtype=float)
    A = np.array([
        [X_u, X_alpha, X_q, -g * np.cos(theta1)],
        [Z_u, Z_alpha, U1 + Z_q, -g * np.sin(theta1)],
        [M_u, M_alpha, M_q, 0],
        [0, 0, 1, 0],
    ], dtype=float)  # fmt: skip
    B = np.array([[X_de], [Z_de], [M_de], [0]], dtype=float)
    return np.linalg.solve(C, A), np.linalg.solve(C, B)


def lateral_state_space(U1, theta1, L_p, L_beta, L_r, Y_p, Y_beta, Y_r, N_p, N_beta, N_r,
                        L_da, L_dr, N_da, N_dr, Y_da=0.0, Y_dr=0.0, g=G_FT):
    """Lateral-directional model, states ``[p, phi, beta, r]``, inputs ``[aileron, rudder]`` (final project)."""
    A = np.array([
        [L_p, 0, L_beta, L_r],
        [1, 0, 0, 0],
        [Y_p / U1, g / U1 * np.cos(theta1), Y_beta / U1, Y_r / U1 - 1],
        [N_p, 0, N_beta, N_r],
    ], dtype=float)  # fmt: skip
    B = np.array([[L_da, L_dr], [0, 0], [Y_da / U1, Y_dr / U1], [N_da, N_dr]], dtype=float)
    return A, B


def transfer_functions(A, B, input_index=0):
    """``(num, den)`` from input ``input_index`` to every state (``ss2tf``)."""
    C = np.eye(A.shape[0])
    D = np.zeros((A.shape[0], B.shape[1]))
    return signal.ss2tf(A, B, C, D, input=input_index)


def step_input(t, start, amplitude=1.0):
    t = np.asarray(t)
    return np.where(t >= start, amplitude, 0.0)


def doublet(t, start, width, amplitude=1.0):
    """+amplitude for ``width`` seconds then -amplitude for ``width`` seconds."""
    t = np.asarray(t)
    return amplitude * (((t >= start) & (t < start + width)).astype(float) - ((t >= start + width) & (t < start + 2 * width)))


def simulate(A, B, u, t, x0=None):
    """State response to input history ``u`` (``lsim`` with C = I)."""
    n = A.shape[0]
    sys = signal.StateSpace(A, B, np.eye(n), np.zeros((n, B.shape[1])))
    _, y, _ = signal.lsim(sys, np.asarray(u, dtype=float), t, X0=x0)
    return y


def discretize(A, B, dt):
    """Zero-order-hold discretisation."""
    n, m = B.shape
    M = expm(np.block([[A, B], [np.zeros((m, n + m))]]) * dt)
    return M[:n, :n], M[:n, n:]


@dataclass
class LinearMPC:
    """Condensed linear MPC with box-constrained inputs (stand-in for the MPC Toolbox controllers of ``mpc1_script``/``mpcscript_final``).

    Minimises ``sum ||W_y (y - r)||^2 + ||W_du du||^2 + ||W_u (u - u_nom)||^2``
    over ``control_horizon`` moves (held to the end of the prediction
    horizon) subject to ``u_min <= u <= u_max``. Output constraints are not
    enforced (the toolbox softened them).
    """

    A: np.ndarray  # continuous
    B: np.ndarray
    C: np.ndarray
    dt: float
    prediction_horizon: int = 10
    control_horizon: int = 2
    w_y: np.ndarray | float = 1.0
    w_du: np.ndarray | float = 0.1
    w_u: np.ndarray | float = 0.0
    u_min: np.ndarray | float = -np.inf
    u_max: np.ndarray | float = np.inf
    u_nominal: np.ndarray | float = 0.0

    def __post_init__(self):
        self.Ad, self.Bd = discretize(np.asarray(self.A, float), np.asarray(self.B, float), self.dt)
        self.C = np.atleast_2d(self.C)
        self.ny, self.nu = self.C.shape[0], self.Bd.shape[1]
        P, Mc = self.prediction_horizon, self.control_horizon
        # y_k = C A^k x0 + sum C A^(k-1-j) B u_j, with u_j = u_{Mc-1} for j >= Mc
        self._F = np.vstack([self.C @ np.linalg.matrix_power(self.Ad, k + 1) for k in range(P)])
        G = np.zeros((P * self.ny, Mc * self.nu))
        for k in range(P):
            for j in range(k + 1):
                block = self.C @ np.linalg.matrix_power(self.Ad, k - j) @ self.Bd
                col = min(j, Mc - 1)
                G[k * self.ny:(k + 1) * self.ny, col * self.nu:(col + 1) * self.nu] += block
        self._G = G
        D = np.eye(Mc * self.nu) - np.eye(Mc * self.nu, k=-self.nu)
        self._D = D

    def _weights(self, w, n, reps):
        w = np.broadcast_to(np.asarray(w, dtype=float), (n,))
        return np.tile(np.sqrt(w), reps)

    def control(self, x, reference, u_prev):
        """Optimal first move given state ``x``, output reference and the previous input."""
        P, Mc = self.prediction_horizon, self.control_horizon
        r = np.tile(np.broadcast_to(reference, (self.ny,)), P)
        wy = self._weights(self.w_y, self.ny, P)
        wdu = self._weights(self.w_du, self.nu, Mc)
        wu = self._weights(self.w_u, self.nu, Mc)
        u_prev = np.broadcast_to(np.asarray(u_prev, float), (self.nu,))
        u_nom = np.tile(np.broadcast_to(np.asarray(self.u_nominal, float), (self.nu,)), Mc)
        e0 = np.zeros(Mc * self.nu)
        e0[: self.nu] = u_prev
        A_ls = np.vstack([wy[:, None] * self._G, wdu[:, None] * self._D, wu[:, None] * np.eye(Mc * self.nu)])
        b_ls = np.concatenate([wy * (r - self._F @ x), wdu * e0, wu * u_nom])
        lo = np.tile(np.broadcast_to(np.asarray(self.u_min, float), (self.nu,)), Mc)
        hi = np.tile(np.broadcast_to(np.asarray(self.u_max, float), (self.nu,)), Mc)
        res = lsq_linear(A_ls, b_ls, bounds=(lo, hi))
        return res.x[: self.nu]

    def simulate(self, x0, reference, steps, u0=0.0):
        """Closed loop on the discretised plant; ``reference(k)`` or a constant. Returns ``(x, u, y)`` histories."""
        x = np.asarray(x0, float)
        u = np.broadcast_to(np.asarray(u0, float), (self.nu,)).copy()
        X, U, Y = [x], [], [self.C @ x]
        for k in range(steps):
            r = reference(k) if callable(reference) else reference
            u = self.control(x, r, u)
            x = self.Ad @ x + self.Bd @ u
            X.append(x)
            U.append(u)
            Y.append(self.C @ x)
        return np.array(X), np.array(U), np.array(Y)
