"""Determine what pyCycle's tabular thermo format stores as gamma and Cp/Cv.

The tabular thermodynamics format used by pyCycle stores one set of specific
heats and a single ``gamma`` per grid point. Before h2thermo can write that
format, two questions have to be answered:

1. Is ``gamma`` the ratio of the stored specific heats, or an independently
   evaluated isentropic exponent?
2. Are the stored specific heats frozen-composition values, or do they include
   the contribution of shifting dissociation equilibria?

An earlier version of this script answered both by building a live pyCycle
model and evaluating its CEA thermodynamic path directly. That measures the
right physics, but it measures the wrong thing: what pyCycle's runtime
computes, not what its export path actually writes to the tabular file that
``export/pycycle.py`` has to reproduce. It also broke across pyCycle versions,
since it depended on wiring internal to pyCycle's ``ThermoAdd``/``Thermo``
components.

This version reads pyCycle's own shipped reference table directly --
``pycycle.constants.AIR_JETA_TAB_SPEC``, loaded from
``pycycle/thermo/tabular/air_jetA.pkl`` -- and evaluates the same three
discriminators against the values actually stored on disk. That is both more
robust (no live OpenMDAO model to break) and more directly relevant (it
measures the artifact h2thermo has to be compatible with).

Three independent discriminators are evaluated at each probe state:

* ``gamma`` against ``Cp / Cv``. Equality to machine precision means ``gamma``
  carries no information beyond the two specific heats.
* ``Cp`` against a finite difference of the tabulated specific enthalpy with
  respect to temperature, at fixed pressure and fuel-air ratio. Because the
  table already spans a full grid, this uses the table's own resolution
  rather than an arbitrarily chosen step.
* ``Cp - Cv`` against the specific gas constant ``R``. For an ideal gas of
  fixed composition the two are equal; once the composition is free to shift
  they are not.

At low temperature there is too little dissociation for the discriminators to
separate, so a low-temperature sanity check is included: if the three
quantities do not collapse to near-equality there, the interpretation below is
wrong. States near 2900 K are reported because that is the edge of the range
h2thermo compares against NASA CEA; the table's own hottest node, 3500 K, lies
outside that compared range and is reported separately for scale.

pyCycle is not a dependency of h2thermo and has to be installed separately:

    pip install om-pycycle

Run from the repository root:

    python scripts/probe_pycycle_definitions.py
"""

from __future__ import annotations

import numpy as np

try:
    from pycycle.constants import AIR_JETA_TAB_SPEC
except ImportError as error:  # pragma: no cover - developer tooling only
    raise SystemExit(
        "pyCycle is required to probe its tabular thermo definitions.\n"
        "Install it with: pip install om-pycycle"
    ) from error

#: Representative fuel-air ratio and pressure at which the temperature sweep
#: is evaluated. Chosen well inside the table's range so that the nearest grid
#: node is a close match rather than an extrapolation.
FUEL_AIR_RATIO = 0.03
PRESSURE = 101325.0

#: Temperatures at which the discriminators are reported. 300 K is a low
#: temperature sanity check, 2900 K is the edge of the range h2thermo compares
#: against NASA CEA, and 3500 K is the table's own hottest node, outside that
#: compared range.
REPORT_TEMPERATURES = (300.0, 1200.0, 1600.0, 2000.0, 2400.0, 2800.0, 2900.0, 3500.0)

#: Relative difference below which two quantities are treated as identical.
EQUALITY_TOLERANCE = 1.0e-2


def nearest_index(axis: np.ndarray, value: float) -> int:
    """Return the index of the grid node closest to ``value``."""
    return int(np.argmin(np.abs(axis - value)))


def enthalpy_temperature_derivative(
    enthalpy: np.ndarray, temperature: np.ndarray, index: int
) -> float:
    """Finite-difference dh/dT at a temperature-axis node.

    Uses a central difference when both neighbouring nodes exist, and a
    one-sided difference at the edges of the table.
    """
    if 0 < index < temperature.size - 1:
        return (enthalpy[index + 1] - enthalpy[index - 1]) / (
            temperature[index + 1] - temperature[index - 1]
        )
    if index == 0:
        return (enthalpy[index + 1] - enthalpy[index]) / (
            temperature[index + 1] - temperature[index]
        )
    return (enthalpy[index] - enthalpy[index - 1]) / (
        temperature[index] - temperature[index - 1]
    )


def relative_difference(value: float, reference: float) -> float:
    """Return the magnitude of the relative difference between two values."""
    return abs(value - reference) / abs(reference)


def probe() -> list[dict]:
    """Evaluate every discriminator at each of :data:`REPORT_TEMPERATURES`."""
    spec = AIR_JETA_TAB_SPEC
    temperature_axis = np.asarray(spec["T"])
    pressure_axis = np.asarray(spec["P"])
    far_axis = np.asarray(spec["FAR"])

    far_index = nearest_index(far_axis, FUEL_AIR_RATIO)
    pressure_index = nearest_index(pressure_axis, PRESSURE)

    enthalpy = np.asarray(spec["h"])[far_index, pressure_index, :]
    gamma_field = np.asarray(spec["gamma"])[far_index, pressure_index, :]
    cp_field = np.asarray(spec["Cp"])[far_index, pressure_index, :]
    cv_field = np.asarray(spec["Cv"])[far_index, pressure_index, :]
    r_field = np.asarray(spec["R"])[far_index, pressure_index, :]

    records = []
    for target_temperature in REPORT_TEMPERATURES:
        index = nearest_index(temperature_axis, target_temperature)
        cp = float(cp_field[index])
        cv = float(cv_field[index])
        records.append(
            {
                "target_T": target_temperature,
                "T": float(temperature_axis[index]),
                "gamma": float(gamma_field[index]),
                "cp_over_cv": cp / cv,
                "Cp": cp,
                "dh_dT": enthalpy_temperature_derivative(
                    enthalpy, temperature_axis, index
                ),
                "cp_minus_cv": cp - cv,
                "R": float(r_field[index]),
            }
        )

    return records


def report(records: list[dict]) -> None:
    """Print the three comparisons at every reported temperature."""
    header = (
        f"{'T [K]':>8}{'gamma':>12}{'Cp/Cv':>12}{'rel. diff.':>13}"
    )
    print("\nDiscriminator 1: gamma against the ratio of specific heats")
    print(header)
    for record in records:
        difference = relative_difference(record["gamma"], record["cp_over_cv"])
        print(
            f"{record['T']:>8.1f}{record['gamma']:>12.6f}"
            f"{record['cp_over_cv']:>12.6f}{difference:>13.2e}"
        )

    header = f"{'T [K]':>8}{'Cp':>12}{'dh/dT':>12}{'rel. diff.':>13}"
    print("\nDiscriminator 2: Cp against the finite-difference dh/dT")
    print(header)
    for record in records:
        difference = relative_difference(record["Cp"], record["dh_dT"])
        print(
            f"{record['T']:>8.1f}{record['Cp']:>12.2f}{record['dh_dT']:>12.2f}"
            f"{difference:>13.2e}"
        )

    header = f"{'T [K]':>8}{'Cp - Cv':>12}{'R':>12}{'rel. diff.':>13}"
    print("\nDiscriminator 3: Cp - Cv against the specific gas constant")
    print(header)
    for record in records:
        difference = relative_difference(record["cp_minus_cv"], record["R"])
        print(
            f"{record['T']:>8.1f}{record['cp_minus_cv']:>12.2f}{record['R']:>12.2f}"
            f"{difference:>13.2e}"
        )


def conclude(records: list[dict]) -> None:
    """State what the measurements imply, at the low and high ends of the sweep."""
    coldest = min(records, key=lambda record: record["T"])
    hottest = max(records, key=lambda record: record["T"])

    print(
        f"\nLow temperature sanity check, at {coldest['T']:.0f} K, where "
        "dissociation is negligible:"
    )
    gamma_ok = (
        relative_difference(coldest["gamma"], coldest["cp_over_cv"])
        < EQUALITY_TOLERANCE
    )
    r_ok = (
        relative_difference(coldest["cp_minus_cv"], coldest["R"])
        < EQUALITY_TOLERANCE
    )
    if gamma_ok and r_ok:
        print(
            "  All three quantities collapse to near-equality, as expected "
            "with no dissociation. The interpretation below is not "
            "contradicted at low temperature."
        )
    else:
        print(
            "  The quantities do NOT collapse at low temperature. Something "
            "in this script's interpretation, or in the table itself, needs "
            "re-checking before trusting the high temperature conclusion."
        )

    gamma_difference = relative_difference(hottest["gamma"], hottest["cp_over_cv"])
    specific_heat_difference = relative_difference(hottest["Cp"], hottest["dh_dT"])
    gas_constant_difference = relative_difference(hottest["cp_minus_cv"], hottest["R"])

    print(f"\nConclusion, drawn at {hottest['T']:.0f} K:")

    if gamma_difference < EQUALITY_TOLERANCE:
        print("  gamma is the ratio of the stored specific heats.")
    else:
        print(
            f"  gamma differs from Cp/Cv by {gamma_difference:.1%}, and is "
            f"below it ({hottest['gamma']:.4f} vs {hottest['cp_over_cv']:.4f}). "
            "A frozen-composition ratio would sit above the equilibrium "
            "Cp/Cv, not below it, since dissociation inflates Cp/Cv relative "
            "to its frozen value. A value below the equilibrium ratio is "
            "therefore inconsistent with gamma being any kind of frozen "
            "ratio, and is consistent with an independently evaluated "
            "isentropic exponent."
        )

    if specific_heat_difference < EQUALITY_TOLERANCE:
        print(
            "  Cp matches the finite-difference dh/dT, so the stored "
            "specific heats include shifting dissociation."
        )
    else:
        print(
            f"  Cp departs from dh/dT by {specific_heat_difference:.1%}, so "
            "the stored specific heats hold the composition fixed."
        )

    if gas_constant_difference < EQUALITY_TOLERANCE:
        print(
            "  Cp - Cv equals R, which is consistent with fixed composition."
        )
    else:
        print(
            f"  Cp - Cv departs from R by {gas_constant_difference:.0%}, "
            "which is only possible with a shifting composition."
        )


def main() -> None:
    records = probe()
    report(records)
    conclude(records)


if __name__ == "__main__":
    main()
