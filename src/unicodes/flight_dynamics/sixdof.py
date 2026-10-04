"""Nonlinear 6-DOF aircraft simulation with linear stability derivatives (AE 722 PAH studies).

Port of Ben Mays' ``BenMays6DOF`` and its descendants (``BenMays6DOF_PAH*``,
``Aether_6DOF``, ``AETHER_Ben_6DOF_with_optimizer*``, ``Roskam_6DOF_C182``,
``Cessna_182_PAH``). One right-hand side covers all of them:

* states ``[V, alpha, beta, phi, theta, psi, P, Q, R, x_N, y_E, z_D,
  dt, de, da, dr]`` plus ``[delta, delta_dot]`` for every passive
  aero-compliant hinged (PAH) flap;
* controls ``[throttle, elevator, aileron, rudder]`` drive first-order
  actuators; flaps take spring stiffness ``k``, damping ``c_r``/``c_c`` and
  neutral angle ``delta_0`` as their "controls";
* optional gusts ``[u_g, v_g, w_g(, p_g, q_g, r_g)]`` in body axes are added
  to the air-relative velocity and rates, as in ``Aether_6DOF``.

Units are whatever the model is given (the course used ft, slug, s). The
MATLAB mixed metre flap chords with ft/s airspeeds; keep the flap in the
same system as the aircraft.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

STATE_NAMES = [
    "V", "alpha", "beta", "phi", "theta", "psi", "P", "Q", "R",
    "x_N", "y_E", "z_D", "throttle", "elevator", "aileron", "rudder",
]  # fmt: skip


@dataclass
class StabilityDerivatives:
    """Non-dimensional stability and control derivatives (per rad), Roskam convention."""

    C_L1: float = 0.0
    C_m1: float = 0.0
    C_mT1: float = 0.0
    C_D0: float = 0.03
    e: float = 0.8
    C_L_alpha: float = 4.5
    C_L_q: float = 0.0
    C_L_u: float = 0.0
    C_L_de: float = 0.0
    C_m_alpha: float = -0.6
    C_m_q: float = -12.0
    C_m_u: float = 0.0
    C_m_de: float = -1.0
    C_mT_u: float = 0.0
    C_mT_alpha: float = 0.0
    C_y_beta: float = 0.0
    C_y_p: float = 0.0
    C_y_r: float = 0.0
    C_y_da: float = 0.0
    C_y_dr: float = 0.0
    C_l_beta: float = 0.0
    C_l_p: float = 0.0
    C_l_r: float = 0.0
    C_l_da: float = 0.0
    C_l_dr: float = 0.0
    C_n_beta: float = 0.0
    C_n_p: float = 0.0
    C_n_r: float = 0.0
    C_n_da: float = 0.0
    C_n_dr: float = 0.0


@dataclass
class Aircraft:
    S: float
    c_bar: float
    b: float
    mass: float
    I_xx: float
    I_yy: float
    I_zz: float
    I_xz: float
    rho: float
    g: float
    thrust: float  # constant thrust
    derivatives: StabilityDerivatives
    V_trim: float
    alpha_trim: float = 0.0
    actuator_rates: tuple = (7.0, 30.0, 30.0, 30.0)  # first-order lags for throttle, elevator, aileron, rudder

    @property
    def aspect_ratio(self) -> float:
        return self.b**2 / self.S


@dataclass
class PAHFlap:
    """A passive aero-compliant hinged flap (``Adaptive_Aerocompliant_Airfoil_Dynamics``).

    Position vectors are from the aircraft c.g. in body axes: hinge-line
    inboard/outboard ends, inboard trailing edge and the hinge-line station
    of the flap c.g. ``tau`` is the flap effectiveness; the flap adds
    ``tau C_L_alpha / 2 * delta`` to the lift and the matching rolling moment.
    """

    chord: float
    r_inboard: np.ndarray
    r_outboard: np.ndarray
    r_inboard_te: np.ndarray
    r_cg: np.ndarray
    I_h: float
    k: float
    c_r: float
    delta_0: float
    c_c: float = 0.0
    C_h0: float = 0.0
    C_h_alpha: float = -0.0584
    C_h_delta: float = -0.363
    tau: float = 0.45
    delta_min: float = -np.deg2rad(20)
    delta_max: float = np.deg2rad(20)
    k_stop: float = 1000.0
    disabled: bool = False

    def __post_init__(self):
        for name in ("r_inboard", "r_outboard", "r_inboard_te", "r_cg"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        x = -(self.r_inboard_te - self.r_inboard)
        x /= np.linalg.norm(x)
        y = -(self.r_inboard - self.r_outboard)
        y /= np.linalg.norm(y)
        self._H_BF = np.vstack([x, y, np.cross(x, y)])  # body -> flap axes

    def hinge_acceleration(self, delta, delta_dot, rho, V_body, omega_body):
        """Angular acceleration of the flap about its hinge."""
        if self.disabled:
            return 0.0
        P, Q, R = omega_body
        omega_tilde = np.array([[0, -R, Q], [R, 0, -P], [-Q, P, 0]])
        V_f = self._H_BF @ (V_body + omega_tilde @ self.r_cg)
        alpha_f = np.arctan2(V_f[2], V_f[0])
        q_f = 0.5 * rho * V_f @ V_f
        M_h = q_f * self.chord**2 * (self.C_h0 + self.C_h_alpha * alpha_f + self.C_h_delta * delta)
        M_stop = -self.k_stop * (min(delta - self.delta_min, 0.0) + max(delta - self.delta_max, 0.0))
        M = -self.k * (delta - self.delta_0) - self.c_r * delta_dot + M_h + M_stop
        if abs(M) <= self.c_c:  # Coulomb friction holds the flap
            return 0.0
        return (M - self.c_c * np.sign(delta_dot)) / self.I_h


def critical_damping(I_h, k, q_bar, chord, C_h_delta, zeta=0.7):
    """Rate damping giving damping ratio ``zeta`` for the flap (the ``c_d`` formula in ``BenMays6DOF_PAH*``)."""
    mu = 4 * I_h * (q_bar * chord**2 * C_h_delta - k)
    return 2 * I_h * zeta * np.sqrt(-mu * (2 * I_h * zeta - 1) * (2 * I_h * zeta + 1)) / (4 * I_h**2 * zeta**2 - 1)


def derivatives(t, x, controls, ac: Aircraft, flaps=(), gust=None):
    """State derivative of the 6-DOF model (the ``aircraft`` nested function)."""
    d = ac.derivatives
    V, alpha, beta, phi, theta, psi, P, Q, R = x[:9]
    dt_x, de_x, da_x, dr_x = x[12:16]
    u = V * np.cos(alpha) * np.cos(beta)
    v = V * np.sin(beta)
    w = V * np.sin(alpha) * np.cos(beta)
    Pa, Qa, Ra = P, Q, R
    if gust is not None:
        u, v, w = u + gust[0], v + gust[1], w + gust[2]
        if len(gust) == 6:
            Pa, Qa, Ra = P + gust[3], Q + gust[4], R + gust[5]

    u_trim = ac.V_trim * np.cos(ac.alpha_trim)
    da_alpha = alpha - ac.alpha_trim
    q_bar = 0.5 * ac.rho * V**2

    flap_lift = sum(0.5 * d.C_L_alpha * f.tau * x[16 + 2 * i] for i, f in enumerate(flaps))
    flap_roll = -sum(0.5 * d.C_L_alpha * f.tau * x[16 + 2 * i] * f.r_cg[1] for i, f in enumerate(flaps))

    C_L = (d.C_L1 + d.C_L_alpha * da_alpha + d.C_L_q * Qa * ac.c_bar / 2 / u_trim
           + d.C_L_u * (u - u_trim) / u_trim + d.C_L_de * de_x + flap_lift)
    C_D = d.C_D0 + C_L**2 / (np.pi * ac.aspect_ratio * d.e)
    hb = ac.b / 2 / u_trim
    C_Y = d.C_y_beta * beta + d.C_y_p * Pa * hb + d.C_y_r * Ra * hb + d.C_y_da * da_x + d.C_y_dr * dr_x
    C_ls = d.C_l_beta * beta + d.C_l_p * Pa * hb + d.C_l_r * Ra * hb + d.C_l_da * da_x + d.C_l_dr * dr_x + flap_roll
    du = (u - u_trim) / u_trim
    C_ms = (d.C_m1 + d.C_m_alpha * da_alpha + d.C_m_q * Qa * ac.c_bar / 2 / u_trim + d.C_m_u * du
            + d.C_m_de * de_x + 2 * d.C_m1 * du + (d.C_mT_u + 2 * d.C_mT1) * du + d.C_mT_alpha * da_alpha)
    C_ns = d.C_n_beta * beta + d.C_n_p * Pa * hb + d.C_n_r * Ra * hb + d.C_n_da * da_x + d.C_n_dr * dr_x

    sa, ca = np.sin(alpha), np.cos(alpha)
    C_x = C_L * sa - C_D * ca
    C_z = -C_L * ca - C_D * sa
    C_l = C_ls * ca - C_ns * sa
    C_n = C_ls * sa + C_ns * ca

    sp, cp, st, ct, ss, cs = np.sin(phi), np.cos(phi), np.sin(theta), np.cos(theta), np.sin(psi), np.cos(psi)
    g = ac.g
    ax = (q_bar * ac.S * C_x + ac.thrust) / ac.mass
    ay = q_bar * ac.S * C_Y / ac.mass
    az = q_bar * ac.S * C_z / ac.mass
    L = C_l * q_bar * ac.S * ac.b
    M = C_ms * q_bar * ac.S * ac.c_bar
    N = C_n * q_bar * ac.S * ac.b

    udot = ax - g * st + R * v - Q * w
    vdot = ay + g * sp * ct - R * u + P * w
    wdot = az + g * cp * ct + Q * u - P * v
    Vdot = (u * udot + v * vdot + w * wdot) / V
    alphadot = (u * wdot - w * udot) / (u**2 + w**2)
    betadot = (V * vdot - v * Vdot) / (u**2 + w**2) * np.cos(beta)

    Ixx, Iyy, Izz, Ixz = ac.I_xx, ac.I_yy, ac.I_zz, ac.I_xz
    den = Ixx * Izz - Ixz**2
    Pdot = (Izz * L + Ixz * N - (Ixz * (Iyy - Ixx - Izz) * P + (Ixz**2 + Izz * (Izz - Iyy)) * R) * Q) / den
    Qdot = (M - (Ixx - Izz) * P * R - Ixz * (P**2 - R**2)) / Iyy
    Rdot = (Ixz * L + Ixx * N + (Ixz * (Iyy - Ixx - Izz) * R + (Ixz**2 + Ixx * (Ixx - Iyy)) * P) * Q) / den

    phidot = P + (R * cp + Q * sp) * np.tan(theta)
    thetadot = Q * cp - R * sp
    psidot = (Q * sp + R * cp) / ct
    xdot = ct * cs * u + (-cp * ss + sp * st * cs) * v + (sp * ss + cp * st * cs) * w
    ydot = ct * ss * u + (sp * st * ss + cp * cs) * v + (cp * st * ss - sp * cs) * w
    zdot = -st * u + sp * ct * v + cp * ct * w

    rates = np.asarray(ac.actuator_rates)
    act = rates * (np.asarray(controls[:4]) - x[12:16])

    out = [Vdot, alphadot, betadot, phidot, thetadot, psidot, Pdot, Qdot, Rdot, xdot, ydot, zdot, *act]
    for i, f in enumerate(flaps):
        delta, delta_dot = x[16 + 2 * i], x[17 + 2 * i]
        if f.disabled:
            out += [0.0, 0.0]
        else:
            out += [delta_dot, f.hinge_acceleration(delta, delta_dot, ac.rho, np.array([u, v, w]), (Pa, Qa, Ra))]
    return np.array(out, dtype=float)


@dataclass
class SimulationResult:
    t: np.ndarray
    x: np.ndarray  # (n_steps, n_states)
    names: list = field(default_factory=list)

    def __getitem__(self, name):
        return self.x[:, self.names.index(name)]

    def inertial_acceleration(self):
        """Second difference of the N/E/D position, the quantity the PAH optimiser minimised."""
        dt = np.diff(self.t).mean()
        return np.diff(self.x[:, 9:12], 2, axis=0) / dt**2


def simulate(ac: Aircraft, t, x0, controls, flaps=(), gusts=None, rtol=1e-6, atol=1e-8):
    """Integrate with zero-order-hold inputs between the samples ``t`` (the MATLAB ``ode45`` loop).

    ``controls`` is ``(len(t), 4)``; ``gusts`` is ``None`` or ``(len(t), 3|6)``.
    """
    t = np.asarray(t, dtype=float)
    controls = np.asarray(controls, dtype=float)
    X = np.zeros((len(t), len(x0)))
    X[0] = x0
    for i in range(len(t) - 1):
        g = None if gusts is None else gusts[i]
        sol = solve_ivp(derivatives, (t[i], t[i + 1]), X[i], args=(controls[i], ac, flaps, g), rtol=rtol, atol=atol)
        X[i + 1] = sol.y[:, -1]
    names = STATE_NAMES + [n for i in range(len(flaps)) for n in (f"delta_f{i + 1}", f"delta_f{i + 1}_dot")]
    return SimulationResult(t, X, names)


def initial_state(ac: Aircraft, theta=None, flaps=(), elevator=0.0):
    """Trimmed-looking initial state: V_trim, alpha_trim, level wings, flaps at their neutral angle."""
    x = np.zeros(16 + 2 * len(flaps))
    x[0], x[1] = ac.V_trim, ac.alpha_trim
    x[4] = ac.alpha_trim if theta is None else theta
    x[13] = elevator
    for i, f in enumerate(flaps):
        x[16 + 2 * i] = f.delta_0
    return x


def with_pah(flaps, k=None, c_r=None, c_c=None, delta_0=None, I_h=None, disabled=None):
    """Copies of ``flaps`` with new spring/damper settings (the optimiser's design variables)."""
    changes = {key: val for key, val in dict(k=k, c_r=c_r, c_c=c_c, delta_0=delta_0, I_h=I_h, disabled=disabled).items() if val is not None}
    return tuple(replace(f, **changes) for f in flaps)


# Reference aircraft used in the course scripts


def cessna_182() -> Aircraft:
    """Cessna 182 at 220 ft/s (``BenMays6DOF``; Roskam Airplane Flight Dynamics Part I appendix)."""
    g = 32.14741
    d = StabilityDerivatives(
        C_L1=0.307, C_D0=0.032, e=0.75, C_L_alpha=4.41, C_L_q=3.9, C_L_de=0.43,
        C_m_alpha=-0.613, C_m_q=-12.4, C_m_de=-1.122,
        C_y_beta=-0.393, C_y_p=-0.075, C_y_r=0.214, C_y_dr=0.187,
        C_l_beta=-0.0923, C_l_p=-0.484, C_l_r=0.0798, C_l_da=0.229, C_l_dr=0.0147,
        C_n_beta=0.0587, C_n_p=-0.0278, C_n_r=-0.0937, C_n_da=-0.0216, C_n_dr=-0.0645,
    )  # fmt: skip
    return Aircraft(S=174, c_bar=4.9, b=36, mass=2650 / g, I_xx=948, I_yy=1346, I_zz=1967, I_xz=0.0,
                    rho=0.002048, g=g, thrust=322.5669, derivatives=d, V_trim=220.1)


def aether() -> Aircraft:
    """AETHER business jet at 46,000 ft, 487.54 kt (``AETHER_Ben_6DOF_with_optimizer_AAA_*``, derivatives from AAA)."""
    g = 31.7153794531126
    d = StabilityDerivatives(
        C_L1=0.532856354160566, C_m1=0.0148998516426128, C_mT1=-0.014899851642622,
        C_D0=0.0168477084000695, e=0.830582171991287,
        C_L_alpha=4.23019118620399, C_L_q=4.52284048284355, C_L_u=0.388809603754571, C_L_de=0.0558642102755981,
        C_m_alpha=-0.243334502306824, C_m_q=-14.7581643366953, C_m_u=-0.195742486769034,
        C_m_de=-0.202966499104865, C_mT_u=0.0452685785311446,
        C_y_beta=-0.419347403822892, C_y_p=-0.0676930444441548, C_y_r=0.204240184135172, C_y_dr=0.0480026076052465,
        C_l_beta=-0.140578157040709, C_l_p=-0.371485760250744, C_l_r=0.217266887233459,
        C_l_da=0.0210299331545958, C_l_dr=0.217266887233459,
        C_n_beta=-0.0190716004532104, C_n_p=-0.0822407440867264, C_n_r=-0.115902396101941,
        C_n_da=-0.00105750416876261, C_n_dr=-0.115902396101941,
    )  # fmt: skip
    kt = 6076.12 / 3600
    return Aircraft(S=932.57, c_bar=14.92, b=76.80, mass=74753 / g, I_xx=510417.5747, I_yy=342708.943,
                    I_zz=856772.3576, I_xz=23697.95883, rho=0.00044, g=g, thrust=6382.09494583848,
                    derivatives=d, V_trim=487.54 * kt, alpha_trim=np.deg2rad(5.4851042558604))


def wing_flaps(semi_span, chord, tau, **kwargs):
    """A right and a left full-semi-span PAH flap (the two-flap layout of the course scripts)."""
    right = PAHFlap(chord, [0, 0, 0], [0, semi_span, 0], [-chord, 0, 0], [0, semi_span / 2, 0], tau=tau, **kwargs)
    left = PAHFlap(chord, [0, -semi_span, 0], [0, 0, 0], [-chord, -semi_span, 0], [0, -semi_span / 2, 0], tau=tau, **kwargs)
    return right, left


def linearize(ac: Aircraft, x0, controls, flaps=(), step=0.01):
    """State Jacobian of the 6-DOF model at ``x0`` (``BenMays_6DOF_jacobian``)."""
    from ..numerics import jacobian

    return jacobian(lambda x: derivatives(0.0, x, controls, ac, flaps), x0, step, relative=False)


def flap_eigenvalues(I_h, k, c_r, q_bar, chord, C_h_delta):
    """Eigenvalues of the linearised flap ``[[0, 1], [-(k - q c^2 C_h_delta)/I_h, -c_r/I_h]]``
    (``adaptive_airfoil_dynamics_linearizedEvaluation``)."""
    A = np.array([[0.0, 1.0], [-(k - q_bar * chord**2 * C_h_delta) / I_h, -c_r / I_h]])
    return np.linalg.eigvals(A)
