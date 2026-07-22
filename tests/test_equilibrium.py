"""Tests for the equilibrium module.

The first test is the permanent form of the original smoke test: it pins the
adiabatic flame temperature of stoichiometric hydrogen-air against the accepted
literature value, and therefore fails loudly if the mechanism, the solver or
the reactant setup ever changes behaviour.
"""

import math

import cantera as ct
import pytest

from h2thermo import (
    adiabatic_flame_temperature,
    create_gas,
    equilibrium_properties,
)

AMBIENT_TEMPERATURE = 298.15  # K
AMBIENT_PRESSURE = ct.one_atm  # Pa

#: Accepted literature value for stoichiometric hydrogen-air combustion
#: initiated at ambient conditions.
REFERENCE_FLAME_TEMPERATURE = 2400.0  # K
FLAME_TEMPERATURE_TOLERANCE = 100.0  # K


@pytest.fixture(scope="module")
def gas() -> ct.Solution:
    """Provide a single solution object shared by all tests in this module."""
    return create_gas()


def test_stoichiometric_flame_temperature_matches_literature(gas):
    """Stoichiometric hydrogen-air must ignite to roughly 2400 K."""
    flame_temperature = adiabatic_flame_temperature(
        AMBIENT_TEMPERATURE, AMBIENT_PRESSURE, equivalence_ratio=1.0, gas=gas
    )
    assert flame_temperature == pytest.approx(
        REFERENCE_FLAME_TEMPERATURE, abs=FLAME_TEMPERATURE_TOLERANCE
    )


def test_flame_temperature_peaks_near_stoichiometric(gas):
    """Lean and rich mixtures must burn cooler than the stoichiometric one."""
    temperatures = {
        phi: adiabatic_flame_temperature(
            AMBIENT_TEMPERATURE, AMBIENT_PRESSURE, equivalence_ratio=phi, gas=gas
        )
        for phi in (0.5, 1.0, 1.5)
    }
    assert temperatures[0.5] < temperatures[1.0]
    assert temperatures[1.5] < temperatures[1.0]


def test_lean_combustion_leaves_excess_oxygen(gas):
    """A fuel lean mixture must retain unburnt oxygen at equilibrium."""
    state = equilibrium_properties(
        1500.0, AMBIENT_PRESSURE, equivalence_ratio=0.5, gas=gas
    )
    assert state.mole_fractions["O2"] > 1.0e-2
    # The fuel is essentially fully consumed; only a trace remains.
    assert state.mole_fractions["H2"] < 1.0e-4


def test_mole_fractions_sum_to_unity(gas):
    """The composition returned by the solver must be normalised."""
    state = equilibrium_properties(
        2000.0, 10.0e5, equivalence_ratio=0.8, gas=gas
    )
    assert math.isclose(sum(state.mole_fractions.values()), 1.0, rel_tol=1.0e-9)


def test_properties_are_physically_plausible(gas):
    """Returned properties must lie within physically meaningful bounds."""
    state = equilibrium_properties(
        1800.0, 20.0e5, equivalence_ratio=0.6, gas=gas
    )
    assert state.cp > 0.0
    assert state.cv > 0.0
    assert 1.0 < state.gamma < 1.7
    assert 0.0 < state.mean_molecular_weight < 50.0
    assert state.density > 0.0


def test_dissociation_increases_with_temperature(gas):
    """Radical concentrations must grow as the product temperature rises."""
    cool = equilibrium_properties(
        1500.0, AMBIENT_PRESSURE, equivalence_ratio=1.0, gas=gas
    )
    hot = equilibrium_properties(
        2800.0, AMBIENT_PRESSURE, equivalence_ratio=1.0, gas=gas
    )
    assert hot.mole_fractions["OH"] > cool.mole_fractions["OH"]
    assert hot.mole_fractions["H"] > cool.mole_fractions["H"]


@pytest.mark.parametrize(
    "kwargs",
    [
        {"temperature": -1.0, "pressure": 1.0e5, "equivalence_ratio": 1.0},
        {"temperature": 1000.0, "pressure": 0.0, "equivalence_ratio": 1.0},
        {"temperature": 1000.0, "pressure": 1.0e5, "equivalence_ratio": -0.5},
    ],
)
def test_invalid_inputs_are_rejected(kwargs, gas):
    """Non-physical inputs must raise instead of returning silent garbage."""
    with pytest.raises(ValueError):
        equilibrium_properties(gas=gas, **kwargs)
