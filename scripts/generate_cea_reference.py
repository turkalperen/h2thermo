"""Generate NASA CEA reference points for validating h2thermo.

This script is run infrequently and is the only part of the project that
requires the NASA CEA package. It writes a small comma separated file of
reference states that the test suite compares against, so that validation runs
in continuous integration without a compiled Fortran dependency.

The CEA product species list is deliberately restricted to the species present
in Cantera's ``h2o2.yaml`` mechanism. This isolates differences in the
thermodynamic databases and the solvers from differences in the species sets,
which would otherwise be conflated. The nitrogen chemistry that CEA offers and
``h2o2.yaml`` lacks is quantified separately in ``docs/validation.md``.

Usage
-----
    pip install cea
    python scripts/generate_cea_reference.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

try:
    from cea import TP, EqSolution, EqSolver, Mixture
except ImportError as error:  # pragma: no cover - developer tooling only
    raise SystemExit(
        "The NASA CEA package is required to regenerate reference data.\n"
        "Install it with: pip install cea"
    ) from error

#: Product species common to CEA and Cantera's h2o2.yaml mechanism.
MATCHED_PRODUCTS = (
    "H2",
    "H",
    "O",
    "O2",
    "OH",
    "H2O",
    "HO2",
    "H2O2",
    "N2",
)

#: Species whose mole fractions are recorded for composition validation.
REPORTED_SPECIES = ("H2O", "N2", "O2", "OH", "H", "O", "H2")

#: Reactant species. Air is expressed explicitly rather than using the CEA
#: built-in, whose element list includes argon and carbon and would therefore
#: not match the product element list.
REACTANTS = ("H2", "O2", "N2")

FUEL_MOLES = (1.0, 0.0, 0.0)
OXIDIZER_MOLES = (0.0, 1.0, 3.76)

REFERENCE_TEMPERATURES = (600.0, 1000.0, 1400.0, 1800.0, 2200.0, 2600.0, 2900.0)
REFERENCE_PRESSURES_BAR = (1.0, 5.0, 20.0, 60.0)
REFERENCE_EQUIVALENCE_RATIOS = (0.2, 0.4, 0.6, 0.8, 1.0)

OUTPUT_PATH = Path("data/cea_reference_points.csv")

BAR_TO_PASCAL = 1.0e5
KILO = 1.0e3


def build_solver() -> tuple[EqSolver, EqSolution, Mixture, np.ndarray, np.ndarray]:
    """Construct the CEA solver and the reactant weight vectors.

    Returns
    -------
    tuple
        The solver, a reusable solution object, the reactant mixture and the
        fuel and oxidizer mass fraction vectors.
    """
    reactants = Mixture(list(REACTANTS))
    products = Mixture(list(MATCHED_PRODUCTS))
    solver = EqSolver(products, reactants=reactants)
    solution = EqSolution(solver)

    fuel_weights = reactants.moles_to_weights(np.asarray(FUEL_MOLES, dtype=float))
    oxidizer_weights = reactants.moles_to_weights(
        np.asarray(OXIDIZER_MOLES, dtype=float)
    )
    return solver, solution, reactants, fuel_weights, oxidizer_weights


def solve_reference_point(
    solver: EqSolver,
    solution: EqSolution,
    reactants: Mixture,
    fuel_weights: np.ndarray,
    oxidizer_weights: np.ndarray,
    temperature: float,
    pressure_bar: float,
    equivalence_ratio: float,
) -> EqSolution | None:
    """Solve a single constant temperature and pressure equilibrium case.

    Returns
    -------
    EqSolution or None
        The converged solution, or ``None`` when CEA failed to converge.
    """
    oxidizer_fuel_ratio = reactants.chem_eq_ratio_to_of_ratio(
        oxidizer_weights, fuel_weights, equivalence_ratio
    )
    weights = reactants.of_ratio_to_weights(
        oxidizer_weights, fuel_weights, oxidizer_fuel_ratio
    )
    solver.solve(solution, TP, temperature, pressure_bar, weights)
    return solution if solution.converged else None


def main() -> None:
    """Sweep the reference grid and write the results to disk."""
    solver, solution, reactants, fuel_weights, oxidizer_weights = build_solver()

    header = [
        "temperature_K",
        "pressure_Pa",
        "equivalence_ratio",
        "mean_molecular_weight_kg_per_kmol",
        "density_kg_per_m3",
        "enthalpy_J_per_kg",
        "entropy_J_per_kg_K",
        "cp_frozen_J_per_kg_K",
        "cp_equilibrium_J_per_kg_K",
        "cv_frozen_J_per_kg_K",
        "cv_equilibrium_J_per_kg_K",
        "isentropic_exponent",
    ] + [f"mole_fraction_{species}" for species in REPORTED_SPECIES]

    rows = []
    failures = 0

    for temperature in REFERENCE_TEMPERATURES:
        for pressure_bar in REFERENCE_PRESSURES_BAR:
            for equivalence_ratio in REFERENCE_EQUIVALENCE_RATIOS:
                result = solve_reference_point(
                    solver,
                    solution,
                    reactants,
                    fuel_weights,
                    oxidizer_weights,
                    temperature,
                    pressure_bar,
                    equivalence_ratio,
                )
                if result is None:
                    failures += 1
                    continue

                mole_fractions = result.mole_fractions
                rows.append(
                    [
                        f"{temperature:.6g}",
                        f"{pressure_bar * BAR_TO_PASCAL:.6g}",
                        f"{equivalence_ratio:.6g}",
                        f"{result.MW:.8g}",
                        f"{result.density:.8g}",
                        f"{result.enthalpy * KILO:.8g}",
                        f"{result.entropy * KILO:.8g}",
                        f"{result.cp_fr * KILO:.8g}",
                        f"{result.cp_eq * KILO:.8g}",
                        f"{result.cv_fr * KILO:.8g}",
                        f"{result.cv_eq * KILO:.8g}",
                        f"{result.gamma_s:.8g}",
                    ]
                    + [
                        f"{float(mole_fractions[species]):.6e}"
                        for species in REPORTED_SPECIES
                    ]
                )

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)

    print(f"Wrote {len(rows)} reference points to {OUTPUT_PATH}")
    if failures:
        print(f"CEA failed to converge at {failures} points, which were skipped")


if __name__ == "__main__":
    main()
