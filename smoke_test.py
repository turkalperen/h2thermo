"""Single-point smoke test for the Cantera equilibrium toolchain.

Computes the adiabatic flame temperature of a stoichiometric hydrogen-air
mixture at ambient conditions and compares it against the accepted literature
value. This verifies that Cantera, the reaction mechanism and the equilibrium
solver are correctly installed before any library code is written.
"""

import cantera as ct

MECHANISM = "h2o2.yaml"
FUEL = "H2"
# Dry air, simplified to its two dominant components on a molar basis.
OXIDIZER = {"O2": 1.0, "N2": 3.76}

EQUIVALENCE_RATIO = 1.0
INLET_TEMPERATURE = 298.15  # K
INLET_PRESSURE = ct.one_atm  # Pa

# Accepted literature value for stoichiometric H2-air at ambient conditions.
EXPECTED_FLAME_TEMPERATURE = 2400.0  # K
TEMPERATURE_TOLERANCE = 100.0  # K

REPORTED_SPECIES = ("H2O", "N2", "O2", "OH", "H", "O", "H2")


def compute_adiabatic_flame_temperature() -> ct.Solution:
    """Equilibrate a stoichiometric hydrogen-air mixture at constant H and P.

    Returns
    -------
    ct.Solution
        The gas object in its equilibrium state, from which the adiabatic
        flame temperature and the product composition can be read.
    """
    gas = ct.Solution(MECHANISM)
    gas.set_equivalence_ratio(EQUIVALENCE_RATIO, FUEL, OXIDIZER)
    gas.TP = INLET_TEMPERATURE, INLET_PRESSURE

    # Constant enthalpy and pressure equilibrium defines the adiabatic
    # flame temperature.
    gas.equilibrate("HP")
    return gas


def main() -> None:
    """Run the smoke test and report the outcome."""
    gas = compute_adiabatic_flame_temperature()
    deviation = abs(gas.T - EXPECTED_FLAME_TEMPERATURE)

    print(f"Mechanism:                  {MECHANISM}")
    print(f"Equivalence ratio:          {EQUIVALENCE_RATIO:.2f}")
    print(f"Inlet state:                {INLET_TEMPERATURE:.2f} K, "
          f"{INLET_PRESSURE / 1e5:.3f} bar")
    print(f"Adiabatic flame temperature: {gas.T:.1f} K "
          f"(expected ~{EXPECTED_FLAME_TEMPERATURE:.0f} K)")

    print("\nEquilibrium mole fractions:")
    for species in REPORTED_SPECIES:
        print(f"  {species:<4s} {gas[species].X[0]:.4e}")

    if deviation <= TEMPERATURE_TOLERANCE:
        print("\nSmoke test PASSED: the toolchain is working as expected.")
    else:
        print(f"\nSmoke test FAILED: deviation of {deviation:.1f} K exceeds "
              f"the {TEMPERATURE_TOLERANCE:.0f} K tolerance.")


if __name__ == "__main__":
    main()
