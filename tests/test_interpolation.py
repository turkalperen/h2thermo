"""Tests for the interpolation module.

Accuracy is verified by comparing interpolated values against direct
equilibrium solves at randomly chosen states strictly inside the tabulated
envelope. Tolerances follow the measured accuracy reported in
``docs/validation.md``, with a margin.
"""

import numpy as np
import pytest

from h2thermo import (
    GridSpecification,
    ThermoInterpolator,
    ThermoTable,
    equilibrium_properties,
)

TEMPERATURE_RANGE = (500.0, 3000.0)
PRESSURE_RANGE = (1.0e5, 60.0e5)
EQUIVALENCE_RATIO_RANGE = (0.2, 1.0)

#: Tolerances on interpolated bulk properties, relative unless noted.
BULK_PROPERTY_TOLERANCE = 1.0e-3
DENSITY_TOLERANCE = 1.0e-3
ENTHALPY_TOLERANCE = 3.0e4  # J/kg, absolute

#: Number of random states used in the accuracy tests.
SAMPLE_COUNT = 60


@pytest.fixture(scope="module")
def table() -> ThermoTable:
    """Generate a table on a logarithmic pressure axis."""
    grid = GridSpecification(
        temperature=np.linspace(*TEMPERATURE_RANGE, 30),
        pressure=np.geomspace(*PRESSURE_RANGE, 12),
        equivalence_ratio=np.linspace(*EQUIVALENCE_RATIO_RANGE, 12),
    )
    return ThermoTable.generate(grid)


@pytest.fixture(scope="module")
def interpolator(table: ThermoTable) -> ThermoInterpolator:
    """Provide an interpolator over the shared table."""
    return ThermoInterpolator(table)


@pytest.fixture(scope="module")
def samples() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Draw random states strictly inside the tabulated envelope."""
    rng = np.random.default_rng(20260722)
    temperature = rng.uniform(600.0, 2900.0, SAMPLE_COUNT)
    pressure = np.exp(
        rng.uniform(np.log(1.1e5), np.log(55.0e5), SAMPLE_COUNT)
    )
    equivalence_ratio = rng.uniform(0.25, 0.98, SAMPLE_COUNT)
    return temperature, pressure, equivalence_ratio


class TestConstruction:
    """Guards on the tables an interpolator will accept."""

    def test_single_node_axis_is_rejected(self):
        grid = GridSpecification(
            temperature=np.array([1000.0]),
            pressure=np.array([1.0e5, 2.0e5]),
            equivalence_ratio=np.array([0.5, 1.0]),
        )
        with pytest.raises(ValueError, match="at least two nodes"):
            ThermoInterpolator(ThermoTable.generate(grid))

    def test_species_names_are_exposed(self, interpolator, table):
        assert interpolator.species_names == table.species_names


class TestAccuracy:
    """Agreement between interpolated and directly computed properties."""

    @staticmethod
    @pytest.fixture(scope="class")
    def comparison(interpolator, samples):
        temperature, pressure, equivalence_ratio = samples
        interpolated = interpolator.lookup(
            temperature, pressure, equivalence_ratio
        )
        exact = [
            equilibrium_properties(float(t), float(p), float(phi))
            for t, p, phi in zip(temperature, pressure, equivalence_ratio)
        ]
        return interpolated, exact

    @pytest.mark.parametrize(
        "name", ["cp", "cv", "gamma", "entropy", "mean_molecular_weight"]
    )
    def test_bulk_properties(self, comparison, name):
        interpolated, exact = comparison
        expected = np.array([getattr(state, name) for state in exact])
        actual = getattr(interpolated, name)
        relative_error = np.abs(actual - expected) / np.abs(expected)
        assert relative_error.max() < BULK_PROPERTY_TOLERANCE

    def test_enthalpy(self, comparison):
        interpolated, exact = comparison
        expected = np.array([state.enthalpy for state in exact])
        assert np.abs(interpolated.enthalpy - expected).max() < ENTHALPY_TOLERANCE

    def test_density_is_recovered_from_the_equation_of_state(self, comparison):
        """Density must beat what direct interpolation of the field achieves."""
        interpolated, exact = comparison
        expected = np.array([state.density for state in exact])
        relative_error = np.abs(interpolated.density - expected) / expected
        assert relative_error.max() < DENSITY_TOLERANCE

    def test_major_species_composition(self, interpolator, samples):
        """Species present in bulk must interpolate to better than one per cent."""
        temperature, pressure, equivalence_ratio = samples
        actual = interpolator.mole_fractions(
            temperature, pressure, equivalence_ratio
        )
        for index, (t, p, phi) in enumerate(
            zip(temperature, pressure, equivalence_ratio)
        ):
            exact = equilibrium_properties(float(t), float(p), float(phi))
            for species in ("H2O", "N2"):
                expected = exact.mole_fractions[species]
                if expected < 1.0e-2:
                    continue
                assert actual[species][index] == pytest.approx(
                    expected, rel=1.0e-2
                ), species


class TestQueryInterface:
    """Scalar and array handling."""

    def test_scalar_query_returns_scalars(self, interpolator):
        state = interpolator.lookup(1500.0, 20.0e5, 0.6)
        assert isinstance(state.cp, float)
        assert isinstance(state.density, float)

    def test_array_query_preserves_shape(self, interpolator):
        temperature = np.array([[1000.0, 1500.0], [2000.0, 2500.0]])
        state = interpolator.lookup(temperature, 10.0e5, 0.5)
        assert state.cp.shape == temperature.shape
        assert state.density.shape == temperature.shape

    def test_arguments_are_broadcast(self, interpolator):
        state = interpolator.lookup(
            np.array([1200.0, 1800.0]), 5.0e5, np.array([0.4, 0.9])
        )
        assert state.gamma.shape == (2,)

    def test_scalar_and_array_agree(self, interpolator):
        scalar = interpolator.lookup(1700.0, 15.0e5, 0.7)
        array = interpolator.lookup(
            np.array([1700.0]), np.array([15.0e5]), np.array([0.7])
        )
        assert array.cp[0] == pytest.approx(scalar.cp)

    def test_mole_fractions_sum_to_unity(self, interpolator):
        fractions = interpolator.mole_fractions(1600.0, 8.0e5, 0.55)
        assert sum(fractions.values()) == pytest.approx(1.0, rel=1.0e-6)


class TestBounds:
    """Behaviour outside the tabulated envelope."""

    def test_out_of_range_raises_by_default(self, interpolator):
        with pytest.raises(ValueError):
            interpolator.lookup(100.0, 20.0e5, 0.6)

    def test_out_of_range_returns_nan_when_permitted(self, table):
        permissive = ThermoInterpolator(table, bounds_error=False)
        assert np.isnan(permissive.lookup(100.0, 20.0e5, 0.6).cp)

    def test_non_positive_pressure_raises(self, interpolator):
        with pytest.raises(ValueError, match="pressure must be strictly positive"):
            interpolator.lookup(1500.0, 0.0, 0.6)


class TestPerformance:
    """The speed advantage that motivates tabulation.

    The comparison reuses a single Cantera solution object, which is the
    fastest way to call the solver directly. Comparing against the naive path,
    where a solution object is created per call, would overstate the benefit.
    """

    def test_batched_lookup_is_far_faster_than_solving(self, interpolator):
        import time

        from h2thermo import create_gas

        rng = np.random.default_rng(1)
        count = 20000
        temperature = rng.uniform(700.0, 2800.0, count)
        pressure = rng.uniform(2.0e5, 50.0e5, count)
        equivalence_ratio = rng.uniform(0.3, 0.95, count)

        start = time.perf_counter()
        interpolator.lookup(temperature, pressure, equivalence_ratio)
        batched_seconds_per_point = (time.perf_counter() - start) / count

        gas = create_gas()
        sample = 200
        start = time.perf_counter()
        for index in range(sample):
            equilibrium_properties(
                float(temperature[index]),
                float(pressure[index]),
                float(equivalence_ratio[index]),
                gas=gas,
            )
        solve_seconds_per_point = (time.perf_counter() - start) / sample

        assert batched_seconds_per_point < 0.05 * solve_seconds_per_point

    def test_scalar_lookup_beats_solving(self, interpolator):
        """The dedicated scalar path must be worth having.

        Cycle codes frequently query one state at a time, so a scalar lookup
        that were slower than a direct solve would defeat the purpose of
        tabulation for exactly the case that matters most.
        """
        import time

        from h2thermo import create_gas

        rng = np.random.default_rng(2)
        count = 500
        temperature = rng.uniform(700.0, 2800.0, count)
        pressure = rng.uniform(2.0e5, 50.0e5, count)
        equivalence_ratio = rng.uniform(0.3, 0.95, count)

        start = time.perf_counter()
        for index in range(count):
            interpolator.lookup(
                float(temperature[index]),
                float(pressure[index]),
                float(equivalence_ratio[index]),
            )
        lookup_seconds_per_point = (time.perf_counter() - start) / count

        gas = create_gas()
        start = time.perf_counter()
        for index in range(count):
            equilibrium_properties(
                float(temperature[index]),
                float(pressure[index]),
                float(equivalence_ratio[index]),
                gas=gas,
            )
        solve_seconds_per_point = (time.perf_counter() - start) / count

        assert lookup_seconds_per_point < solve_seconds_per_point
