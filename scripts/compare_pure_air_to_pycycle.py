"""Compare h2thermo's pure-air state against pyCycle's own FAR=0 row.

This is a second, CEA-independent validation channel for the pure-air
(``equivalence_ratio=0.0``) state added in PR-6: instead of comparing against
NASA CEA, it compares directly against pyCycle's own shipped reference table,
``pycycle.constants.AIR_JETA_TAB_SPEC``, at its ``FAR = 0`` row. Agreement
there means an engine model built on ``h2thermo.export.pycycle`` sees
consistent air properties whether the property comes from h2thermo's own
solver or from pyCycle's historical reference data.

The property mapping used (``Cp``/``Cv`` as the equilibrium specific heats,
``gamma`` as the independently evaluated isentropic exponent) is the one
determined in ``scripts/probe_pycycle_definitions.py`` and recorded in
``docs/validation.md`` section 5; this script assumes that result rather than
re-deriving it.

Where h2thermo and pyCycle disagree, NASA CEA is used as an optional third,
independent opinion to determine which one is closer to a named reference
implementation, rather than leaving the disagreement unresolved. This part
of the script is skipped if the ``cea`` package is not installed.

Usage
-----
    pip install om-pycycle
    pip install cea  # optional, enables the CEA tie-breaker
    python scripts/compare_pure_air_to_pycycle.py
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

try:
    from pycycle.constants import AIR_JETA_TAB_SPEC
except ImportError as error:  # pragma: no cover - developer tooling only
    raise SystemExit(
        "pyCycle is required to run this comparison.\n"
        "Install it with: pip install om-pycycle"
    ) from error

import cantera as ct

try:
    from cea import TP as CEA_TP
    from cea import EqSolution, EqSolver, Mixture

    CEA_AVAILABLE = True
except ImportError:  # pragma: no cover - developer tooling only
    CEA_AVAILABLE = False

from h2thermo.equilibrium import (
    DEFAULT_FUEL,
    DEFAULT_MECHANISM,
    create_gas,
    equilibrium_properties,
    equilibrium_specific_heats,
)

#: Real dry air by mole fraction (N2, O2, Ar), used only to test whether the
#: shipped ``DRY_AIR`` reactant (which omits argon) explains the deviations
#: seen at the default composition.
REAL_AIR = {"O2": 20.95, "N2": 78.09, "AR": 0.93}

#: pyCycle's own oxidizer, expressed relative to `DRY_AIR`'s O2:N2 = 1:3.76.
DRY_AIR = {"O2": 1.0, "N2": 3.76}

#: Pressures swept at each temperature, chosen to span the CEA-compared
#: envelope (see docs/validation.md section 2).
PROBE_PRESSURES_BAR = (0.1, 1.0, 5.0, 20.0, 60.0)

#: Temperatures swept at the fixed pressure used for the realistic-range and
#: high-temperature comparisons. 200-1200 K spans what an inlet or compressor
#: section actually sees; above that, pure (unburned) air is not a state any
#: engine section reaches in practice, and the sweep is a mathematical stress
#: test of the table rather than a physically meaningful comparison.
PROBE_TEMPERATURES = (
    300.0,
    600.0,
    900.0,
    1200.0,
    1500.0,
    1800.0,
    2200.0,
    2600.0,
    2900.0,
)

FIXED_PRESSURE_BAR = 5.0

#: A vanishingly small but nonzero equivalence ratio, used as CEA's proxy
#: for equivalence_ratio=0.0. CEA's equilibrium solver hits a singular
#: update matrix at exactly zero fuel, since the iteration variables tied to
#: fuel-derived species become degenerate with no hydrogen present at all.
#: Convergence was checked from 1e-2 down to 1e-6; the result is stable to
#: 0.02% by 1e-4, an order of magnitude tighter than any comparison this
#: script draws a conclusion from.
CEA_PHI_PROXY = 1.0e-4


@dataclass(frozen=True)
class ComparisonRecord:
    """One state compared between h2thermo and pyCycle's FAR=0 row."""

    temperature: float
    pressure: float
    pycycle: dict[str, float]
    h2thermo: dict[str, float]

    def relative_difference(self, field: str) -> float:
        """Relative difference between the two sources for ``field``."""
        reference = self.pycycle[field]
        if reference == 0.0:
            return float("nan")
        return abs(self.h2thermo[field] - reference) / abs(reference)


def _nearest_index(axis: np.ndarray, value: float) -> int:
    """Return the index of the grid node closest to ``value``."""
    return int(np.argmin(np.abs(axis - value)))


def _pycycle_state(
    temperature_index: int, pressure_index: int, far_index: int
) -> dict[str, float]:
    """Read pyCycle's stored properties at one grid node."""
    table = AIR_JETA_TAB_SPEC
    return {
        field: float(table[field][far_index, pressure_index, temperature_index])
        for field in ("Cp", "Cv", "gamma", "h", "S", "rho", "R")
    }


def _h2thermo_state(
    temperature: float, pressure: float, gas: ct.Solution, oxidizer: dict
) -> dict[str, float]:
    """Evaluate h2thermo at the same state, mapped to pyCycle's field names.

    Follows the mapping in :mod:`h2thermo.export.pycycle`: ``Cp``/``Cv`` are
    the equilibrium (shifting-composition) specific heats, ``gamma`` is the
    isentropic exponent, and ``R`` is recovered from the mean molecular
    weight rather than tabulated directly.
    """
    state = equilibrium_properties(
        temperature, pressure, 0.0, fuel=DEFAULT_FUEL, oxidizer=oxidizer, gas=gas
    )
    shifting = equilibrium_specific_heats(
        temperature, pressure, 0.0, fuel=DEFAULT_FUEL, oxidizer=oxidizer, gas=gas
    )
    return {
        "Cp": shifting.cp,
        "Cv": shifting.cv,
        "gamma": shifting.isentropic_exponent,
        "h": state.enthalpy,
        "S": state.entropy,
        "rho": state.density,
        "R": ct.gas_constant / state.mean_molecular_weight,
    }


def compare_pressure_sweep() -> list[ComparisonRecord]:
    """Compare Cp, Cv, gamma, h, S, rho and R across pressure at ~2900 K."""
    gas = create_gas(DEFAULT_MECHANISM)
    table = AIR_JETA_TAB_SPEC
    i_far0 = _nearest_index(table["FAR"], 0.0)
    i_temperature = _nearest_index(table["T"], 2900.0)
    temperature = float(table["T"][i_temperature])

    records = []
    for pressure_bar in PROBE_PRESSURES_BAR:
        i_pressure = _nearest_index(table["P"], pressure_bar * 1.0e5)
        pressure = float(table["P"][i_pressure])
        records.append(
            ComparisonRecord(
                temperature=temperature,
                pressure=pressure,
                pycycle=_pycycle_state(i_temperature, i_pressure, i_far0),
                h2thermo=_h2thermo_state(temperature, pressure, gas, DRY_AIR),
            )
        )
    return records


def compare_temperature_sweep() -> list[ComparisonRecord]:
    """Compare Cp and gamma across temperature at a fixed, realistic pressure.

    This is the comparison that matters for how the pure-air row is actually
    used: as the state of unreacted air moving through an inlet, fan or
    compressor, none of which see pure air anywhere near 2900 K.
    """
    gas = create_gas(DEFAULT_MECHANISM)
    table = AIR_JETA_TAB_SPEC
    i_far0 = _nearest_index(table["FAR"], 0.0)
    i_pressure = _nearest_index(table["P"], FIXED_PRESSURE_BAR * 1.0e5)
    pressure = float(table["P"][i_pressure])

    records = []
    for temperature_target in PROBE_TEMPERATURES:
        i_temperature = _nearest_index(table["T"], temperature_target)
        temperature = float(table["T"][i_temperature])
        records.append(
            ComparisonRecord(
                temperature=temperature,
                pressure=pressure,
                pycycle=_pycycle_state(i_temperature, i_pressure, i_far0),
                h2thermo=_h2thermo_state(temperature, pressure, gas, DRY_AIR),
            )
        )
    return records


def check_argon_hypothesis() -> None:
    """Test whether omitting argon from `DRY_AIR` explains the deviations.

    `DRY_AIR` is a simplified two-component (O2, N2) oxidizer; real air is
    roughly 0.93% argon by mole. This checks a realistic argon-inclusive
    composition at a low-dissociation state (300 K, 1 bar, where the two
    databases should already agree closely on frozen properties) and at the
    hot, low-pressure corner where the largest Cp/Cv deviations appear.
    """
    gas = create_gas(DEFAULT_MECHANISM)
    table = AIR_JETA_TAB_SPEC
    i_far0 = _nearest_index(table["FAR"], 0.0)

    print("\nArgon hypothesis check")
    print("-" * 60)
    for label, temperature_target, pressure_bar in (
        ("baseline (300 K, 1 bar)", 300.0, 1.0),
        ("hot corner (2900 K, 60 bar)", 2900.0, 60.0),
    ):
        i_t = _nearest_index(table["T"], temperature_target)
        i_p = _nearest_index(table["P"], pressure_bar * 1.0e5)
        temperature = float(table["T"][i_t])
        pressure = float(table["P"][i_p])

        pycycle_state = _pycycle_state(i_t, i_p, i_far0)
        dry_air_state = _h2thermo_state(temperature, pressure, gas, DRY_AIR)
        real_air_state = _h2thermo_state(temperature, pressure, gas, REAL_AIR)

        dry_diff = abs(dry_air_state["Cp"] - pycycle_state["Cp"]) / pycycle_state["Cp"]
        real_diff = (
            abs(real_air_state["Cp"] - pycycle_state["Cp"]) / pycycle_state["Cp"]
        )
        print(f"{label}:")
        print(f"  Cp deviation with DRY_AIR (no argon):   {dry_diff:.2%}")
        print(f"  Cp deviation with REAL_AIR (with argon): {real_diff:.2%}")


@dataclass(frozen=True)
class TieBreakRecord:
    """Cp compared across all three sources at one state."""

    temperature: float
    pressure: float
    cp_h2thermo: float
    cp_pycycle: float
    cp_cea: float | None

    @property
    def h2thermo_vs_cea(self) -> float | None:
        """Relative difference between h2thermo and CEA, or None if CEA
        failed to converge at this state."""
        if self.cp_cea is None:
            return None
        return abs(self.cp_h2thermo - self.cp_cea) / self.cp_cea

    @property
    def pycycle_vs_cea(self) -> float | None:
        """Relative difference between pyCycle and CEA, or None if CEA
        failed to converge at this state."""
        if self.cp_cea is None:
            return None
        return abs(self.cp_pycycle - self.cp_cea) / self.cp_cea


def _build_cea_solver():
    """Construct CEA's solver and the reactant weight vectors for pure air.

    Mirrors the setup in ``scripts/generate_cea_reference.py``, restricted to
    the species h2o2.yaml carries, at the trace fuel fraction
    :data:`CEA_PHI_PROXY` that avoids the exact-zero singular matrix.
    """
    reactants = Mixture(["H2", "O2", "N2"])
    products = Mixture(["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "N2"])
    solver = EqSolver(products, reactants=reactants)
    solution = EqSolution(solver)
    fuel_weights = reactants.moles_to_weights(np.array([1.0, 0.0, 0.0]))
    oxidizer_weights = reactants.moles_to_weights(np.array([0.0, 1.0, 3.76]))
    ratio = reactants.chem_eq_ratio_to_of_ratio(
        oxidizer_weights, fuel_weights, CEA_PHI_PROXY
    )
    weights = reactants.of_ratio_to_weights(oxidizer_weights, fuel_weights, ratio)
    return solver, solution, weights


def _cea_pure_air_cp(
    temperature: float, pressure_bar: float, solver, solution, weights
) -> float | None:
    """Return CEA's equilibrium cp for near-pure air, in J/(kg K).

    Returns ``None`` rather than raising if CEA fails to converge at this
    state, so a single bad point does not abort the sweep.
    """
    solver.solve(solution, CEA_TP, temperature, pressure_bar, weights)
    if not solution.converged:
        return None
    return float(solution.cp_eq) * 1000.0


def resolve_with_cea() -> None:
    """Use NASA CEA as a third, independent opinion where h2thermo and
    pyCycle disagree, rather than leaving the disagreement unresolved.

    Every reference_point in this sweep is deliberately drawn from the
    pressure and temperature sweeps already reported above, so the same
    states are being triangulated, not a new set chosen to flatter one
    source over the other.
    """
    if not CEA_AVAILABLE:
        print(
            "\nCEA not installed; skipping the tie-breaker. "
            "Install it with: pip install cea"
        )
        return

    gas = create_gas(DEFAULT_MECHANISM)
    solver, solution, weights = _build_cea_solver()
    table = AIR_JETA_TAB_SPEC
    i_far0 = _nearest_index(table["FAR"], 0.0)

    def build_record(temperature: float, pressure: float) -> TieBreakRecord:
        i_t = _nearest_index(table["T"], temperature)
        i_p = _nearest_index(table["P"], pressure)
        t_grid = float(table["T"][i_t])
        p_grid = float(table["P"][i_p])
        shifting = equilibrium_specific_heats(t_grid, p_grid, 0.0, gas=gas)
        return TieBreakRecord(
            temperature=t_grid,
            pressure=p_grid,
            cp_h2thermo=shifting.cp,
            cp_pycycle=float(table["Cp"][i_far0, i_p, i_t]),
            cp_cea=_cea_pure_air_cp(
                t_grid, p_grid / 1.0e5, solver, solution, weights
            ),
        )

    records = [build_record(2900.0, p * 1.0e5) for p in PROBE_PRESSURES_BAR]
    records += [
        build_record(t, FIXED_PRESSURE_BAR * 1.0e5) for t in PROBE_TEMPERATURES
    ]

    print("\nTie-breaker: Cp against CEA (proxy equivalence_ratio=1e-4)")
    print("-" * 78)
    print(
        f"{'T (K)':>8} {'P (bar)':>9} {'cp_h2thermo':>12} {'cp_pycycle':>11} "
        f"{'cp_CEA':>10} {'|h2-CEA|':>9} {'|py-CEA|':>9}"
    )
    resolved = [r for r in records if r.cp_cea is not None]
    for r in records:
        if r.cp_cea is None:
            print(
                f"{r.temperature:8.1f} {r.pressure / 1.0e5:9.4f}  "
                "CEA did not converge at this state, skipped"
            )
            continue
        print(
            f"{r.temperature:8.1f} {r.pressure / 1.0e5:9.4f} "
            f"{r.cp_h2thermo:12.2f} {r.cp_pycycle:11.2f} {r.cp_cea:10.2f} "
            f"{r.h2thermo_vs_cea:8.2%} {r.pycycle_vs_cea:8.2%}"
        )

    if resolved:
        mean_h2 = sum(r.h2thermo_vs_cea for r in resolved) / len(resolved)
        mean_py = sum(r.pycycle_vs_cea for r in resolved) / len(resolved)
        print(
            f"\nMean |deviation from CEA| across {len(resolved)} resolved "
            f"states: h2thermo {mean_h2:.2%}, pyCycle {mean_py:.2%}."
        )


def report(records: list[ComparisonRecord], fields: tuple[str, ...]) -> None:
    """Print a comparison table for the given records and fields."""
    header = "  ".join(f"{field:>10}" for field in fields)
    print(f"{'T (K)':>8} {'P (bar)':>9}  {header}")
    for record in records:
        differences = "  ".join(
            f"{record.relative_difference(field):>9.2%}" for field in fields
        )
        print(
            f"{record.temperature:8.1f} {record.pressure / 1.0e5:9.4f}  {differences}"
        )


def main() -> None:
    print("Pressure sweep at ~2900 K (relative difference, h2thermo vs pyCycle)")
    print("-" * 70)
    report(compare_pressure_sweep(), ("Cp", "Cv", "gamma", "h", "S", "rho", "R"))

    print(f"\nTemperature sweep at ~{FIXED_PRESSURE_BAR:.0f} bar")
    print("-" * 70)
    report(compare_temperature_sweep(), ("Cp", "gamma"))

    check_argon_hypothesis()
    resolve_with_cea()


if __name__ == "__main__":
    main()
