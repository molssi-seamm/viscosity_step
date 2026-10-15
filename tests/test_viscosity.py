#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of the viscosity analysis, with synthetic stresses of known viscosity."""

import numpy as np
import pytest
from scipy.signal import lfilter

from seamm_md_analysis import H5MDWriter
from seamm_util import Q_

from viscosity_step.viscosity import Viscosity

# A water-like system: 30 000 Å^3 at 300 K, the stress decorrelating in 0.2 ps
T = 300.0
V = 30000.0
TAU = 0.2e-12  # s
ETA = 0.9e-3  # Pa s
DT = Q_(5.0, "fs")


def ou(rng, n, dt, tau, sigma):
    """An Ornstein-Uhlenbeck process: <x(0)x(t)> = sigma^2 exp(-t/tau)."""
    phi = np.exp(-dt / tau)
    noise = rng.normal(scale=sigma * np.sqrt(1 - phi**2), size=n)
    noise[0] = rng.normal(scale=sigma)
    return lfilter([1.0], [1.0, -phi], noise)


def pressure_tensor(rng, n):
    """A pressure tensor (Pa) whose five shear components each give ETA.

    Each component's Green-Kubo integral is sigma^2 tau = ETA k_B T / V.
    """
    k_B = 1.380649e-23
    sigma = np.sqrt(ETA * k_B * T / (V * 1e-30) / TAU)
    X = [ou(rng, n, DT.m_as("s"), TAU, sigma) for _ in range(5)]
    P = np.zeros((n, 3, 3))
    P[:, 0, 1] = P[:, 1, 0] = X[0]
    P[:, 0, 2] = P[:, 2, 0] = X[1]
    P[:, 1, 2] = P[:, 2, 1] = X[2]
    # (Pxx - Pyy)/2 = X3 and (Pyy - Pzz)/2 = X4, about a mean pressure of 1 atm
    P[:, 0, 0] = 2 * X[3]
    P[:, 2, 2] = -2 * X[4]
    P[:, [0, 1, 2], [0, 1, 2]] += 101325.0
    return P


def node():
    step = Viscosity.__new__(Viscosity)
    step._reset()
    return step


def test_known_viscosity_from_blocks():
    rng = np.random.default_rng(2026)
    step = node()
    data = (pressure_tensor(rng, 200_000), DT, T, V)
    blocks = step._split_into_blocks(data, 5)
    assert len(blocks) == 5
    for block in blocks:
        step.process_run(block)
    results = step._fit()
    gk = results["Green-Kubo"][-1]
    helfand = results["Helfand"][-1]
    assert gk["value"] == pytest.approx(0.9, rel=0.1)  # mPa*s
    assert helfand["value"] == pytest.approx(0.9, rel=0.15)
    # The error from the spread of the blocks is sensible
    assert 0.0 < gk["stderr_runs"] < 0.1
    # Each of the five components on its own
    for fit in results["Green-Kubo"][:5]:
        assert fit["value"] == pytest.approx(0.9, rel=0.25)


def test_isotropic_pressure_does_not_matter():
    # The shear components do not see the hydrostatic pressure
    rng = np.random.default_rng(7)
    P = pressure_tensor(rng, 20_000)
    a = node()
    a.process_run((P, DT, T, V))
    b = node()
    b.process_run((P + 5.0e7 * np.eye(3), DT, T, V))
    np.testing.assert_allclose(a.integrals[0], b.integrals[0], rtol=1e-6, atol=1e-12)
    np.testing.assert_allclose(a.Ms[0], b.Ms[0], rtol=1e-6, atol=1e-12)


def test_too_few_samples_for_blocks():
    step = node()
    data = (np.zeros((5000, 3, 3)), DT, T, V)
    with pytest.raises(RuntimeError, match="too few for 5 blocks"):
        step._split_into_blocks(data, 5)
    assert len(step._split_into_blocks(data, 2)) == 2


def write(path, periodic=True, tensor=True, volume=False):
    steps = np.arange(0, 500, 5)
    times = steps * 1.0
    cell = np.diag([30.0, 20.0, 50.0])
    with H5MDWriter(path) as w:
        w.particles(ids=[1, 2], cell=cell if periodic else None, periodic=periodic)
        if tensor:
            P = np.zeros((len(steps), 3, 3))
            P[:, 0, 1] = P[:, 1, 0] = steps
            w.observable("pressure_tensor", steps, times, P, "atm")
        w.observable("temperature", steps, times, np.full(len(steps), 310.0), "K")
        if volume:
            w.observable(
                "volume", steps, times, np.full(len(steps), 28000.0), "Angstrom3"
            )


def test_reads_tensor_in_pascal(tmp_path):
    write(tmp_path / "md.h5md")
    tensor, dt, T_, V_ = node()._read_md_file(tmp_path / "md.h5md")
    assert tensor.shape == (100, 3, 3)
    assert tensor[2, 0, 1] == pytest.approx(10 * 101325.0)
    assert dt.m_as("fs") == pytest.approx(5.0)
    assert T_ == pytest.approx(310.0)
    assert V_ == pytest.approx(30000.0)


def test_npt_volume_is_used(tmp_path):
    write(tmp_path / "md.h5md", volume=True)
    assert node()._read_md_file(tmp_path / "md.h5md")[3] == pytest.approx(28000.0)


def test_non_periodic_is_refused(tmp_path):
    write(tmp_path / "md.h5md", periodic=False)
    with pytest.raises(RuntimeError, match="not periodic"):
        node()._read_md_file(tmp_path / "md.h5md")


def test_missing_tensor(tmp_path):
    write(tmp_path / "md.h5md", tensor=False)
    with pytest.raises(RuntimeError, match="no pressure tensor"):
        node()._read_md_file(tmp_path / "md.h5md")


def test_requirements():
    req = node().md_requirements()
    assert req.periodic
    (need,) = req.needs
    assert need.names == ["pressure_tensor"]
    assert need.max_time.m_as("fs") == 5
