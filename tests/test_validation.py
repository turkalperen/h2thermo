"""Validation of h2thermo against NASA CEA reference points.

The reference states in ``data/cea_reference_points.csv`` were produced by
``scripts/generate_cea_reference.py`` using the NASA CEA package, with the CEA
product species list restricted to those present in Cantera's ``h2o2.yaml``
mechanism. Comparing against the stored file rather than calling CEA directly
keeps this suite runnable in continuous integration without a compiled Fortran
dependency.

Tolerances are set from the observed agreement with a safety margin, and are
deliberately tiered: bulk thermodynamic properties agree far more closely than
radical mole fractions, whose equilibrium concentrations are sensitive to small
differences between the two thermodynamic databases.
"""

import csv
from pathlib import Path

import pytest

from h2thermo import equilibrium_properties

REFERENCE_PATH = (
    Path(__file__).resolve().parent.parent / "data" / "cea_reference_points.csv"
)

#: Relative tolerances on bulk properties.
MOLECULAR_WEIGHT_TOLERANCE = 1.0e-3
DENSITY_TOLERANCE = 1.0e-3
ENTROPY_TOLERANCE = 1.0e-3
FROZEN_CP_TOLERANCE = 5.0e-3

#: Absolute tolerance on specific enthalpy in J/kg. A relative tolerance is
#: unsuitable because the absolute enthalpy passes through zero within the
#: tabulated range. The value below corresponds to roughly ten kelvin of
#: equivalent temperature error at a typical specific heat.
ENTHALPY_TOLERANCE = 2.0e4

#: Species whose equilibrium concentration is set primarily by the overall
#: stoichiometry, and which therefore agree closely.
STABLE_SPECIES = ("H2O", "N2", "O2", "H2")

#: Radicals, whose concentrations depend exponentially on small differences in
#: the underlying Gibbs energy data.
RADICAL_SPECIES = ("OH", "H", "O")

STABLE_SPECIES_TOLERANCE = 5.0e-2
RADICAL_SPECIES_TOLERANCE = 2.0e-1

#: Mole fractions below this value are excluded, as trace species carry no
#: practical weight and their relative deviation is dominated by round off.
MINIMUM_REPORTED_MOLE_FRACTION = 1.0e-4


def load_reference_points() -> list[dict[str, str]]:
    """Read the stored CEA reference points.

    Returns
    -------
    list of dict
        One dictionary per reference state, keyed by column name.
    """
    if not REFERENCE_PATH.exists():
        pytest.skip(
            f"reference data not found at {REFERENCE_PATH}; regenerate it with "
            "scripts/generate_cea_reference.py"
        )
    with REFERENCE_PATH.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


REFERENCE_POINTS = load_reference_points()


def describe(point: dict[str, str]) -> str:
    """Return a short label identifying a reference state."""
    return (
        f"T={float(point['temperature_K']):.0f}K"
        f"_P={float(point['pressure_Pa']) / 1.0e5:.0f}bar"
        f"_phi={float(point['equivalence_ratio']):.1f}"
    )


@pytest.fixture(params=REFERENCE_POINTS, ids=describe)
def reference_point(request) -> dict[str, str]:
    """Provide each stored reference state in turn."""
    return request.param


@pytest.fixture
def computed(reference_point):
    """Evaluate h2thermo at the state of the current reference point."""
    return equilibrium_properties(
        float(reference_point["temperature_K"]),
        float(reference_point["pressure_Pa"]),
        float(reference_point["equivalence_ratio"]),
    )


def test_reference_data_is_present():
    """The reference file must contain a meaningful number of states."""
    assert len(REFERENCE_POINTS) >= 100


def test_mean_molecular_weight(reference_point, computed):
    expected = float(reference_point["mean_molecular_weight_kg_per_kmol"])
    assert computed.mean_molecular_weight == pytest.approx(
        expected, rel=MOLECULAR_WEIGHT_TOLERANCE
    )


def test_density(reference_point, computed):
    expected = float(reference_point["density_kg_per_m3"])
    assert computed.density == pytest.approx(expected, rel=DENSITY_TOLERANCE)


def test_entropy(reference_point, computed):
    expected = float(reference_point["entropy_J_per_kg_K"])
    assert computed.entropy == pytest.approx(expected, rel=ENTROPY_TOLERANCE)


def test_enthalpy(reference_point, computed):
    expected = float(reference_point["enthalpy_J_per_kg"])
    assert computed.enthalpy == pytest.approx(expected, abs=ENTHALPY_TOLERANCE)


def test_frozen_specific_heat(reference_point, computed):
    """The tabulated cp is a frozen-composition value, matching CEA's cp_fr."""
    expected = float(reference_point["cp_frozen_J_per_kg_K"])
    assert computed.cp == pytest.approx(expected, rel=FROZEN_CP_TOLERANCE)


def test_stable_species_composition(reference_point, computed):
    for species in STABLE_SPECIES:
        expected = float(reference_point[f"mole_fraction_{species}"])
        if expected <= MINIMUM_REPORTED_MOLE_FRACTION:
            continue
        assert computed.mole_fractions[species] == pytest.approx(
            expected, rel=STABLE_SPECIES_TOLERANCE
        ), species


def test_radical_species_composition(reference_point, computed):
    for species in RADICAL_SPECIES:
        expected = float(reference_point[f"mole_fraction_{species}"])
        if expected <= MINIMUM_REPORTED_MOLE_FRACTION:
            continue
        assert computed.mole_fractions[species] == pytest.approx(
            expected, rel=RADICAL_SPECIES_TOLERANCE
        ), species


def test_tabulated_cp_is_the_frozen_value_not_the_equilibrium_one():
    """Document the known limitation as an executable assertion.

    At high temperature and low pressure the equilibrium specific heat reported
    by CEA is several times the frozen value that h2thermo currently tabulates.
    This test pins that expectation so the difference cannot be forgotten, and
    will need updating when shifting specific heats are added.
    """
    hot_and_low_pressure = [
        point
        for point in REFERENCE_POINTS
        if float(point["temperature_K"]) >= 2600.0
        and float(point["pressure_Pa"]) <= 1.0e5
        and float(point["equivalence_ratio"]) >= 0.8
    ]
    assert hot_and_low_pressure, "expected high temperature reference points"

    for point in hot_and_low_pressure:
        frozen = float(point["cp_frozen_J_per_kg_K"])
        equilibrium = float(point["cp_equilibrium_J_per_kg_K"])
        assert equilibrium > 1.5 * frozen


def test_dissociation_contribution_to_cp_falls_with_pressure():
    """Raising the pressure suppresses dissociation and its heat capacity.

    Le Chatelier's principle requires the shifting contribution to the specific
    heat to shrink as pressure rises, because dissociation increases the number
    of moles. Reproducing this trend is evidence that the reference data
    describes the physics rather than an arbitrary dataset.
    """
    ratios = {}
    for point in REFERENCE_POINTS:
        if (
            float(point["temperature_K"]) != 2900.0
            or float(point["equivalence_ratio"]) != 1.0
        ):
            continue
        pressure = float(point["pressure_Pa"])
        ratios[pressure] = float(point["cp_equilibrium_J_per_kg_K"]) / float(
            point["cp_frozen_J_per_kg_K"]
        )

    assert len(ratios) >= 3
    pressures = sorted(ratios)
    values = [ratios[pressure] for pressure in pressures]
    assert values == sorted(values, reverse=True)
