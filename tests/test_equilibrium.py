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
from h2thermo.equilibrium import equilibrium_specific_heats

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


class TestPureOxidizer:
    """Equivalence ratio of zero: pure oxidizer, no fuel present.

    This is the state pyCycle's tabular thermo format calls its ``FAR = 0``
    row, needed for engine sections such as an inlet or compressor that
    never see fuel.
    """

    #: Accepted handbook value for the specific heat of dry air near room
    #: temperature. The simplified two-component oxidizer used here (O2/N2
    #: only, no argon) agrees to within half a percent.
    REFERENCE_AIR_CP = 1005.0  # J/(kg K)

    def test_matches_the_oxidizer_composition_exactly(self, gas):
        """With no fuel, equilibrium composition is just the oxidizer."""
        state = equilibrium_properties(
            300.0, AMBIENT_PRESSURE, equivalence_ratio=0.0, gas=gas
        )
        assert state.mole_fractions["O2"] == pytest.approx(1.0 / 4.76, rel=1.0e-6)
        assert state.mole_fractions["N2"] == pytest.approx(3.76 / 4.76, rel=1.0e-6)
        for species in ("H2", "H", "OH", "H2O", "HO2", "H2O2"):
            assert state.mole_fractions[species] < 1.0e-12

    def test_specific_heat_matches_air_near_room_temperature(self, gas):
        """Below the onset of dissociation this is just the cp of air."""
        state = equilibrium_properties(
            300.0, AMBIENT_PRESSURE, equivalence_ratio=0.0, gas=gas
        )
        assert state.cp == pytest.approx(self.REFERENCE_AIR_CP, rel=1.0e-2)

    def test_oxygen_dissociates_at_high_temperature(self, gas):
        """Even fuel-free air dissociates once it is hot enough."""
        cool = equilibrium_properties(
            300.0, AMBIENT_PRESSURE, equivalence_ratio=0.0, gas=gas
        )
        hot = equilibrium_properties(
            2900.0, AMBIENT_PRESSURE, equivalence_ratio=0.0, gas=gas
        )
        assert hot.mole_fractions["O"] > cool.mole_fractions["O"]
        assert hot.mole_fractions["O2"] < cool.mole_fractions["O2"]

    def test_equilibrium_specific_heats_accepts_zero(self, gas):
        """The shifting-composition path must accept a fuel-free mixture too."""
        shifting = equilibrium_specific_heats(
            2900.0, AMBIENT_PRESSURE, equivalence_ratio=0.0, gas=gas
        )
        assert shifting.cp > 0.0
        assert shifting.cv > 0.0
        assert 1.0 < shifting.isentropic_exponent < 1.7

    def test_adiabatic_flame_temperature_equals_inlet(self, gas):
        """With no fuel, nothing burns, so the temperature cannot change."""
        inlet_temperature = 500.0
        result = adiabatic_flame_temperature(
            inlet_temperature, AMBIENT_PRESSURE, equivalence_ratio=0.0, gas=gas
        )
        assert result == pytest.approx(inlet_temperature, abs=1.0e-6)

    def test_negative_equivalence_ratio_is_still_rejected(self, gas):
        """Zero is a valid boundary; negative fuel content is not."""
        with pytest.raises(ValueError):
            equilibrium_properties(
                1000.0, AMBIENT_PRESSURE, equivalence_ratio=-0.1, gas=gas
            )


class TestEquilibriumSpecificHeats:
    """Specific heats that account for a shifting chemical equilibrium."""

    def test_equals_the_frozen_value_when_nothing_dissociates(self, gas):
        """Below the onset of dissociation the two definitions must coincide."""
        state = equilibrium_properties(
            800.0, AMBIENT_PRESSURE, equivalence_ratio=0.5, gas=gas
        )
        shifting = equilibrium_specific_heats(
            800.0, AMBIENT_PRESSURE, equivalence_ratio=0.5, gas=gas
        )
        assert shifting.cp == pytest.approx(state.cp, rel=1.0e-3)
        assert shifting.cv == pytest.approx(state.cv, rel=1.0e-3)
        assert shifting.isentropic_exponent == pytest.approx(
            state.gamma, rel=1.0e-3
        )

    def test_exceeds_the_frozen_value_once_dissociation_matters(self, gas):
        """Breaking bonds absorbs energy, so the equilibrium value is larger."""
        state = equilibrium_properties(
            2800.0, AMBIENT_PRESSURE, equivalence_ratio=1.0, gas=gas
        )
        shifting = equilibrium_specific_heats(
            2800.0, AMBIENT_PRESSURE, equivalence_ratio=1.0, gas=gas
        )
        assert shifting.cp > 2.0 * state.cp

    def test_contribution_falls_with_pressure(self, gas):
        """Pressure suppresses dissociation, and with it the extra capacity."""
        ratios = []
        for pressure in (1.0e5, 10.0e5, 60.0e5):
            state = equilibrium_properties(
                2800.0, pressure, equivalence_ratio=1.0, gas=gas
            )
            shifting = equilibrium_specific_heats(
                2800.0, pressure, equivalence_ratio=1.0, gas=gas
            )
            ratios.append(shifting.cp / state.cp)
        assert ratios == sorted(ratios, reverse=True)

    @pytest.mark.parametrize("step", [0.2, 1.0, 5.0])
    def test_result_is_insensitive_to_the_temperature_step(self, gas, step):
        """A derivative that moved with the step size would be unreliable."""
        reference = equilibrium_specific_heats(
            2400.0, AMBIENT_PRESSURE, equivalence_ratio=1.0, gas=gas
        )
        varied = equilibrium_specific_heats(
            2400.0,
            AMBIENT_PRESSURE,
            equivalence_ratio=1.0,
            gas=gas,
            temperature_step=step,
        )
        assert varied.cp == pytest.approx(reference.cp, rel=1.0e-3)

    def test_specific_heats_stay_ordered(self, gas):
        shifting = equilibrium_specific_heats(
            2200.0, 20.0e5, equivalence_ratio=0.7, gas=gas
        )
        assert shifting.cp > shifting.cv > 0.0
        assert 1.0 < shifting.isentropic_exponent < 1.7

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"temperature": 0.5, "pressure": 1.0e5, "equivalence_ratio": 1.0},
            {"temperature": 1000.0, "pressure": -1.0, "equivalence_ratio": 1.0},
            {
                "temperature": 1000.0,
                "pressure": 1.0e5,
                "equivalence_ratio": 1.0,
                "relative_pressure_step": 2.0,
            },
        ],
    )
    def test_invalid_inputs_are_rejected(self, gas, kwargs):
        with pytest.raises(ValueError):
            equilibrium_specific_heats(gas=gas, **kwargs)
