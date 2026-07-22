"""Equilibrium thermodynamic properties for hydrogen-fuelled combustion.

``h2thermo`` generates validated thermodynamic property tables for hydrogen-air
combustion products, in a form suitable for gas turbine cycle analysis tools
such as pyCycle and T-MATS. The architecture is designed to extend to other
alternative fuels, including sustainable aviation fuels, ammonia and methane.
"""

from h2thermo.equilibrium import (
    DEFAULT_FUEL,
    DEFAULT_MECHANISM,
    DRY_AIR,
    EquilibriumProperties,
    adiabatic_flame_temperature,
    create_gas,
    equilibrium_properties,
)

__version__ = "0.1.0.dev0"

__all__ = [
    "DEFAULT_FUEL",
    "DEFAULT_MECHANISM",
    "DRY_AIR",
    "EquilibriumProperties",
    "__version__",
    "adiabatic_flame_temperature",
    "create_gas",
    "equilibrium_properties",
]
