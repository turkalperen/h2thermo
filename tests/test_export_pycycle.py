"""Tests for the pyCycle tabular thermo export adapter.

pyCycle is not installed in this environment; these tests check that the
written file matches the schema measured in
``scripts/probe_pycycle_definitions.py`` (dictionary keys, axis order, unit
conventions), not that pyCycle itself can load it. See docs/validation.md
section 5 for how that schema and the property mapping were determined.
"""

import pickle

import numpy as np
import pytest

from h2thermo.export.pycycle import (
    PYCYCLE_PROPERTY_SOURCES,
    stoichiometric_fuel_air_ratio,
    write_pycycle_table,
)
from h2thermo.table import GridSpecification, ThermoTable

TEST_SHAPE = (4, 3, 3)


@pytest.fixture(scope="module")
def table() -> ThermoTable:
    """A small table spanning the intended operating envelope."""
    grid = GridSpecification.linear(
        temperature_range=(500.0, 2800.0),
        pressure_range=(1.0e5, 40.0e5),
        equivalence_ratio_range=(0.2, 1.0),
        shape=TEST_SHAPE,
    )
    return ThermoTable.generate(grid)


@pytest.fixture(scope="module")
def table_with_pure_air() -> ThermoTable:
    """A table whose equivalence_ratio axis includes the pure-air row."""
    grid = GridSpecification.linear(
        temperature_range=(300.0, 2800.0),
        pressure_range=(1.0e5, 40.0e5),
        equivalence_ratio_range=(0.0, 1.0),
        shape=(4, 3, 3),
    )
    return ThermoTable.generate(grid)


class TestStoichiometricFuelAirRatio:
    """Conversion from equivalence ratio to pyCycle's fuel-air ratio."""

    def test_hydrogen_dry_air_matches_the_measured_value(self):
        # Measured directly from the mechanism: gas.set_equivalence_ratio(1.0,
        # "H2", DRY_AIR) gives a hydrogen mass fraction of 0.028522..., so
        # FAR_stoichiometric = Y_H2 / (1 - Y_H2).
        assert stoichiometric_fuel_air_ratio() == pytest.approx(0.029360, rel=1.0e-4)

    def test_scales_linearly_with_a_second_call(self):
        # FAR_stoichiometric depends only on the mechanism and reactants, not
        # on any grid, so repeated calls must agree exactly.
        first = stoichiometric_fuel_air_ratio()
        second = stoichiometric_fuel_air_ratio()
        assert first == second


class TestWritePycycleTable:
    """Schema and content of the written pickle file."""

    def test_appends_the_expected_suffix(self, table, tmp_path):
        written = write_pycycle_table(table, tmp_path / "table_without_suffix")
        assert written.suffix == ".pkl"
        assert written.exists()

    def test_payload_has_the_keys_pycycle_expects(self, table, tmp_path):
        path = write_pycycle_table(table, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        expected_keys = {"T", "P", "FAR", "R"} | set(PYCYCLE_PROPERTY_SOURCES)
        assert set(payload) == expected_keys

    def test_axes_are_one_dimensional_and_properties_are_far_p_t(
        self, table, tmp_path
    ):
        path = write_pycycle_table(table, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        n_t, n_p, n_far = table.grid.shape
        assert payload["T"].shape == (n_t,)
        assert payload["P"].shape == (n_p,)
        assert payload["FAR"].shape == (n_far,)
        for name in set(PYCYCLE_PROPERTY_SOURCES) | {"R"}:
            assert payload[name].shape == (n_far, n_p, n_t)

    def test_fuel_air_ratio_is_the_equivalence_ratio_scaled_linearly(
        self, table, tmp_path
    ):
        path = write_pycycle_table(table, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        far_stoichiometric = stoichiometric_fuel_air_ratio(
            fuel=table.metadata["fuel"],
            oxidizer=table.metadata["oxidizer"],
            mechanism=table.metadata["mechanism"],
        )
        expected = table.grid.equivalence_ratio * far_stoichiometric
        assert np.allclose(payload["FAR"], expected)

    def test_properties_are_reordered_and_mapped_correctly(self, table, tmp_path):
        path = write_pycycle_table(table, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        for pycycle_name, h2thermo_name in PYCYCLE_PROPERTY_SOURCES.items():
            expected = np.transpose(
                table.properties[h2thermo_name], axes=(2, 1, 0)
            )
            assert np.allclose(payload[pycycle_name], expected), pycycle_name

    def test_specific_gas_constant_matches_the_equation_of_state(
        self, table, tmp_path
    ):
        path = write_pycycle_table(table, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        mean_molecular_weight = np.transpose(
            table.properties["mean_molecular_weight"], axes=(2, 1, 0)
        )
        expected = 8314.46261815324 / mean_molecular_weight
        assert np.allclose(payload["R"], expected, rtol=1.0e-8)

    def test_non_convergent_table_is_rejected(self, table, tmp_path):
        import copy

        broken = copy.deepcopy(table)
        broken.properties["cp"][0, 0, 0] = np.nan

        with pytest.raises(ValueError, match="non-convergent"):
            write_pycycle_table(broken, tmp_path / "table.pkl")


class TestPureAirRowExport:
    """The FAR=0 row pyCycle's own air/Jet-A table has, and a full engine
    model needs for unburned sections such as an inlet or compressor."""

    def test_far_axis_starts_at_zero(self, table_with_pure_air, tmp_path):
        path = write_pycycle_table(table_with_pure_air, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        assert payload["FAR"][0] == 0.0

    def test_far_zero_slice_has_no_negative_specific_heats(
        self, table_with_pure_air, tmp_path
    ):
        # A minimal sanity check that the pure-air row survived the export
        # mapping intact, without re-deriving the full physics already
        # covered by TestPureOxidizer in tests/test_equilibrium.py.
        path = write_pycycle_table(table_with_pure_air, tmp_path / "table.pkl")
        with path.open("rb") as handle:
            payload = pickle.load(handle)

        assert np.all(payload["Cp"][0, :, :] > 0.0)
        assert np.all(payload["Cv"][0, :, :] > 0.0)
        assert np.all(payload["R"][0, :, :] > 0.0)
