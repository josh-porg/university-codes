import numpy as np
from scipy.integrate import quad

from unicodes import remote_sensing as rs
from unicodes.neural import HIGHAM_X, HIGHAM_Y, SigmoidNetwork
from unicodes.optimize import test_functions as tf


def test_known_minima():
    assert tf.rosenbrock(np.ones(3)) == 0
    assert np.isclose(tf.goldstein_price(np.array([0.0, -1.0])), 3)
    assert np.isclose(tf.branin(np.array([np.pi, 2.275])), 0.397887, atol=1e-6)
    assert np.isclose(tf.camel6(np.array([0.0898, -0.7126])), -1.0316, atol=1e-4)
    assert np.isclose(tf.ackley(np.zeros(2)), 0, atol=1e-12)
    assert np.isclose(tf.styblinski_tang(np.full(2, -2.903534)), -39.16599 * 2, atol=1e-4)
    X = np.random.default_rng(0).random((2, 3, 4))
    for fs in tf.CATEGORIES.values():
        for f in fs:
            assert f(X).shape == (3, 4)


def test_planck_and_wien():
    T = 300.0
    lam = np.linspace(2e-6, 40e-6, 2000)
    assert abs(lam[np.argmax(rs.planck_spectral_radiance(T, lam))] - rs.wien_peak_wavelength(T)) < 0.05e-6
    total, _ = quad(lambda l: rs.planck_spectral_radiance(T, l), 1e-7, 1e-3, limit=200)
    assert np.isclose(np.pi * total, 5.670374e-8 * T**4, rtol=1e-3)  # Stefan-Boltzmann


def test_radar_equation_inverse():
    args = (100.0, 1e3, 1e3, 0.15, 5e3, 1e-8, 150.0)
    P = rs.radar_received_power(args[0], args[1], args[2], args[3], 10.0, args[4], args[5], args[6])
    assert np.isclose(rs.rcs_for_received_power(P, *args[:4], *args[4:]), 10.0)


def test_backprop_gradient_matches_finite_difference():
    net = SigmoidNetwork.random(rng=1)
    dW, db = net.gradients(HIGHAM_X, HIGHAM_Y)
    g = np.concatenate([w.ravel(order="F") for w in dW] + [b.ravel() for b in db])
    p0 = net.pack()
    h = 1e-6
    fd = np.array([(net.unpack(p0 + h * e).cost(HIGHAM_X, HIGHAM_Y) - net.unpack(p0 - h * e).cost(HIGHAM_X, HIGHAM_Y))
                   / (2 * h) for e in np.eye(p0.size)])
    assert np.allclose(2 * g, fd, atol=1e-6)  # cost is the full sum of squares
    assert np.allclose(net.unpack(p0).pack(), p0)
