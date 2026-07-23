"""Tests for the thermodynamic table module."""

import numpy as np
import pytest

from h2thermo.equilibrium import (
    equilibrium_properties,
    equilibrium_specific_heats,
)
from h2thermo.table import (
    EQUILIBRIUM_PROPERTY_SOURCES,
    PROPERTY_NAMES,
    GridSpecification,
    ThermoTable,
)

TEST_SHAPE = (4, 3, 3)


@pytest.fixture(scope="module")
def grid() -> GridSpecification:
    """Provide a small grid that spans the intended operating envelope."""
    return GridSpecification.linear(
        temperature_range=(500.0, 2800.0),
        pressure_range=(1.0e5, 40.0e5),
        equivalence_ratio_range=(0.2, 1.0),
        shape=TEST_SHAPE,
    )


@pytest.fixture(scope="module")
def table(grid: GridSpecification) -> ThermoTable:
    """Generate the table once and share it across tests in this module."""
    return ThermoTable.generate(grid)


class TestGridSpecification:
    """Validation and construction of the sampling grid."""

    def test_linear_produces_requested_shape(self, grid):
        assert grid.shape == TEST_SHAPE
        assert grid.size == int(np.prod(TEST_SHAPE))

    def test_axes_span_the_requested_bounds(self, grid):
        assert grid.temperature[0] == pytest.approx(500.0)
        assert grid.temperature[-1] == pytest.approx(2800.0)
        assert grid.equivalence_ratio[-1] == pytest.approx(1.0)

    @pytest.mark.parametrize(
        "temperature",
        [
            [1000.0, 900.0, 1100.0],  # not increasing
            [-100.0, 500.0],  # not positive
            [],  # empty
        ],
    )
    def test_invalid_axes_are_rejected(self, temperature):
        with pytest.raises(ValueError):
            GridSpecification(
                temperature=np.asarray(temperature, dtype=float),
                pressure=np.array([1.0e5]),
                equivalence_ratio=np.array([1.0]),
            )

    def test_equivalence_ratio_may_start_at_zero(self):
        """Zero (pure oxidizer) is a valid lower bound, unlike temperature
        and pressure, which must stay strictly positive."""
        grid = GridSpecification(
            temperature=np.array([300.0]),
            pressure=np.array([1.0e5]),
            equivalence_ratio=np.array([0.0, 0.5, 1.0]),
        )
        assert grid.equivalence_ratio[0] == 0.0

    def test_negative_equivalence_ratio_is_rejected(self):
        with pytest.raises(ValueError):
            GridSpecification(
                temperature=np.array([300.0]),
                pressure=np.array([1.0e5]),
                equivalence_ratio=np.array([-0.1, 0.5]),
            )


class TestGeneration:
    """Population of the table from the equilibrium solver."""

    def test_every_node_converges(self, table):
        assert table.failed_node_count == 0

    def test_property_arrays_match_grid_shape(self, table):
        for name in PROPERTY_NAMES:
            assert table.properties[name].shape == TEST_SHAPE

    def test_mole_fraction_array_has_species_axis(self, table):
        assert table.mole_fractions.shape == TEST_SHAPE + (
            len(table.species_names),
        )

    def test_mole_fractions_sum_to_unity_at_every_node(self, table):
        totals = table.mole_fractions.sum(axis=-1)
        assert np.allclose(totals, 1.0, rtol=1.0e-9)

    def test_tabulated_values_match_a_direct_solver_call(self, table, grid):
        """A spot check that the table stores what the solvers returned."""
        i, j, k = 2, 1, 2
        state = (
            float(grid.temperature[i]),
            float(grid.pressure[j]),
            float(grid.equivalence_ratio[k]),
        )
        reference = equilibrium_properties(*state)
        shifting = equilibrium_specific_heats(*state)

        for name in PROPERTY_NAMES:
            source = EQUILIBRIUM_PROPERTY_SOURCES.get(name)
            expected = (
                getattr(shifting, source)
                if source is not None
                else getattr(reference, name)
            )
            assert table.properties[name][i, j, k] == pytest.approx(expected), name

    def test_equilibrium_specific_heat_exceeds_the_frozen_one(self, table):
        """Shifting equilibrium can only add heat capacity, never remove it."""
        frozen = table.properties["cp"]
        equilibrium = table.properties["cp_equilibrium"]
        assert np.all(equilibrium >= frozen - 1.0e-6)

    def test_isentropic_exponent_stays_physical(self, table):
        exponent = table.properties["isentropic_exponent"]
        assert np.all(exponent > 1.0)
        assert np.all(exponent < 1.7)

    def test_metadata_records_provenance(self, table):
        assert table.metadata["fuel"] == "H2"
        assert table.metadata["mechanism"] == "h2o2.yaml"
        assert "generated_at" in table.metadata


class TestPhysicalBehaviour:
    """Trends the tabulated data must reproduce."""

    def test_gamma_stays_within_physical_bounds(self, table):
        gamma = table.properties["gamma"]
        assert np.all(gamma > 1.0)
        assert np.all(gamma < 1.7)

    def test_water_content_rises_towards_stoichiometric(self, table):
        """Richer mixtures up to phi = 1 must produce more water vapour."""
        water = table.mole_fraction_of("H2O")
        # Compare the leanest and the stoichiometric slice at every T and P.
        assert np.all(water[..., -1] > water[..., 0])

    def test_dissociation_rises_with_temperature(self, table):
        """Radical content must increase along the temperature axis."""
        hydroxyl = table.mole_fraction_of("OH")
        assert np.all(hydroxyl[-1, ...] > hydroxyl[0, ...])


class TestPureAirRow:
    """Generation over a grid whose equivalence_ratio axis includes zero."""

    @staticmethod
    @pytest.fixture(scope="class")
    def grid_with_pure_air() -> GridSpecification:
        return GridSpecification.linear(
            temperature_range=(300.0, 2900.0),
            pressure_range=(1.0e5, 20.0e5),
            equivalence_ratio_range=(0.0, 1.0),
            shape=(3, 2, 3),
        )

    @staticmethod
    @pytest.fixture(scope="class")
    def table_with_pure_air(grid_with_pure_air: GridSpecification) -> ThermoTable:
        return ThermoTable.generate(grid_with_pure_air)

    def test_every_node_still_converges(self, table_with_pure_air):
        assert table_with_pure_air.failed_node_count == 0

    def test_pure_air_slice_has_no_fuel_derived_species(self, table_with_pure_air):
        """The k=0 slice is equivalence_ratio=0.0: no hydrogen anywhere."""
        for species in ("H2", "H", "OH", "H2O"):
            fraction = table_with_pure_air.mole_fraction_of(species)
            assert np.all(fraction[:, :, 0] < 1.0e-12)

    def test_pure_air_slice_is_mostly_oxygen_and_nitrogen(self, table_with_pure_air):
        oxygen = table_with_pure_air.mole_fraction_of("O2")[:, :, 0]
        nitrogen = table_with_pure_air.mole_fraction_of("N2")[:, :, 0]
        # At 2900 K a few per cent of O2 dissociates into atomic O, so the
        # bound is loosened at the hot end rather than tightened artificially.
        assert np.all(oxygen + nitrogen > 0.95)


class TestPersistence:
    """Round trip through the on-disk format."""

    def test_save_appends_the_expected_suffix(self, table, tmp_path):
        written = table.save(tmp_path / "table_without_suffix")
        assert written.suffix == ".npz"
        assert written.exists()

    def test_round_trip_preserves_all_data(self, table, tmp_path):
        restored = ThermoTable.load(table.save(tmp_path / "table.npz"))

        assert restored.species_names == table.species_names
        assert restored.metadata == table.metadata
        assert np.allclose(restored.grid.temperature, table.grid.temperature)
        assert np.allclose(restored.grid.pressure, table.grid.pressure)
        assert np.allclose(
            restored.grid.equivalence_ratio, table.grid.equivalence_ratio
        )
        assert np.allclose(restored.mole_fractions, table.mole_fractions)
        for name in PROPERTY_NAMES:
            assert np.allclose(restored.properties[name], table.properties[name])

    def test_incompatible_format_version_is_rejected(self, table, tmp_path):
        import json

        path = table.save(tmp_path / "table.npz")
        with np.load(path, allow_pickle=False) as archive:
            arrays = {key: archive[key] for key in archive.files}

        metadata = json.loads(str(arrays["metadata"]))
        metadata["file_format_version"] = 999
        arrays["metadata"] = np.asarray(json.dumps(metadata))
        np.savez_compressed(path, **arrays)

        with pytest.raises(ValueError, match="unsupported file format version"):
            ThermoTable.load(path)


class TestSpeciesAccess:
    """Lookup helpers for the species axis."""

    def test_mole_fraction_of_returns_a_grid_shaped_field(self, table):
        assert table.mole_fraction_of("N2").shape == TEST_SHAPE

    def test_unknown_species_raises(self, table):
        with pytest.raises(KeyError):
            table.mole_fraction_of("CO2")
