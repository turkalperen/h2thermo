"""Equilibrium combustion calculations built on Cantera.

This module provides a thin, well-documented wrapper around Cantera's Gibbs
free energy minimisation solver. It exposes the two operations required to
build thermodynamic property tables for engine cycle analysis:

* :func:`equilibrium_properties` evaluates the equilibrium composition and the
  associated thermodynamic properties at a prescribed temperature, pressure
  and equivalence ratio.
* :func:`adiabatic_flame_temperature` determines the temperature reached by an
  adiabatic, constant pressure combustion process.

All quantities use SI units on a mass basis unless stated otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import cantera as ct
import numpy as np

__all__ = [
    "DEFAULT_MECHANISM",
    "DEFAULT_FUEL",
    "DEFAULT_RELATIVE_PRESSURE_STEP",
    "DEFAULT_TEMPERATURE_STEP",
    "DRY_AIR",
    "EquilibriumProperties",
    "EquilibriumSpecificHeats",
    "adiabatic_flame_temperature",
    "create_gas",
    "equilibrium_properties",
    "equilibrium_specific_heats",
]

#: Reaction mechanism shipped with Cantera, based on the hydrogen kinetics of
#: Burke et al. (2012). Only the species list and the thermodynamic polynomials
#: are relevant for equilibrium calculations; reaction rates are not used.
DEFAULT_MECHANISM = "h2o2.yaml"

#: Default fuel species.
DEFAULT_FUEL = "H2"

#: Dry air simplified to its two dominant components, on a molar basis.
DRY_AIR: Mapping[str, float] = {"O2": 1.0, "N2": 3.76}

#: Half width in K of the central difference used for temperature derivatives.
DEFAULT_TEMPERATURE_STEP = 1.0

#: Half width of the central difference used for pressure derivatives,
#: expressed as a fraction of the pressure.
DEFAULT_RELATIVE_PRESSURE_STEP = 1.0e-3


@dataclass(frozen=True)
class EquilibriumSpecificHeats:
    """Specific heats that account for a shifting chemical equilibrium.

    Attributes
    ----------
    cp : float
        Specific heat at constant pressure in J/(kg K), including the
        contribution of the shifting composition.
    cv : float
        Specific heat at constant volume in J/(kg K), likewise.
    isentropic_exponent : float
        The exponent relating pressure and specific volume along an isentrope.
        For a reacting mixture this differs from ``cp / cv`` and is the
        quantity required by compressible flow relations.
    """

    cp: float
    cv: float
    isentropic_exponent: float


@dataclass(frozen=True)
class EquilibriumProperties:
    """Thermodynamic state of a combustion product mixture at equilibrium.

    Attributes
    ----------
    temperature : float
        Static temperature in K.
    pressure : float
        Static pressure in Pa.
    equivalence_ratio : float
        Fuel-to-air equivalence ratio of the original reactant mixture.
    enthalpy : float
        Specific enthalpy in J/kg, referenced to the Cantera convention in
        which the enthalpy of formation is included.
    entropy : float
        Specific entropy in J/(kg K).
    cp : float
        Specific heat capacity at constant pressure in J/(kg K).
    cv : float
        Specific heat capacity at constant volume in J/(kg K).
    gamma : float
        Ratio of specific heats, ``cp / cv``.
    mean_molecular_weight : float
        Mean molecular weight of the mixture in kg/kmol.
    density : float
        Mixture density in kg/m^3.
    mole_fractions : dict of str to float
        Equilibrium mole fraction of every species in the mechanism.

    Notes
    -----
    The specific heats are frozen-composition values evaluated at the
    equilibrium composition. They therefore exclude the additional heat
    capacity contributed by shifting dissociation equilibria, which becomes
    significant above roughly 2000 K. Reporting equilibrium (shifting)
    specific heats is planned for a future release.
    """

    temperature: float
    pressure: float
    equivalence_ratio: float
    enthalpy: float
    entropy: float
    cp: float
    cv: float
    gamma: float
    mean_molecular_weight: float
    density: float
    mole_fractions: dict[str, float] = field(repr=False)


def create_gas(mechanism: str = DEFAULT_MECHANISM) -> ct.Solution:
    """Instantiate a Cantera solution object for the given mechanism.

    Parameters
    ----------
    mechanism : str, optional
        Name or path of the mechanism file.

    Returns
    -------
    ct.Solution
        A freshly created gas object.

    Notes
    -----
    Creating a :class:`~cantera.Solution` is comparatively expensive. Callers
    that evaluate many grid points should create the object once and reuse it
    by passing it through the ``gas`` argument of the functions in this module.
    """
    return ct.Solution(mechanism)


def _prepare_reactants(
    gas: ct.Solution,
    equivalence_ratio: float,
    fuel: str,
    oxidizer: Mapping[str, float],
) -> None:
    """Set the composition of ``gas`` to the specified reactant mixture."""
    if equivalence_ratio <= 0.0:
        raise ValueError(
            f"equivalence_ratio must be positive, got {equivalence_ratio}"
        )
    gas.set_equivalence_ratio(equivalence_ratio, fuel, dict(oxidizer))


def equilibrium_properties(
    temperature: float,
    pressure: float,
    equivalence_ratio: float,
    fuel: str = DEFAULT_FUEL,
    oxidizer: Mapping[str, float] = DRY_AIR,
    gas: ct.Solution | None = None,
    mechanism: str = DEFAULT_MECHANISM,
) -> EquilibriumProperties:
    """Evaluate equilibrium properties at a prescribed temperature and pressure.

    The reactant mixture defined by ``fuel``, ``oxidizer`` and
    ``equivalence_ratio`` is equilibrated at constant temperature and pressure.
    This is the operation performed at every node of a thermodynamic property
    table.

    Parameters
    ----------
    temperature : float
        Temperature in K.
    pressure : float
        Pressure in Pa.
    equivalence_ratio : float
        Fuel-to-air equivalence ratio. Values below unity are fuel lean.
    fuel : str, optional
        Fuel species name as defined in the mechanism.
    oxidizer : mapping of str to float, optional
        Oxidizer composition on a molar basis.
    gas : ct.Solution, optional
        Existing solution object to reuse. A new one is created when omitted.
    mechanism : str, optional
        Mechanism used when ``gas`` is not supplied.

    Returns
    -------
    EquilibriumProperties
        The converged equilibrium state.

    Raises
    ------
    ValueError
        If ``temperature``, ``pressure`` or ``equivalence_ratio`` is not
        strictly positive.
    """
    if temperature <= 0.0:
        raise ValueError(f"temperature must be positive, got {temperature}")
    if pressure <= 0.0:
        raise ValueError(f"pressure must be positive, got {pressure}")

    if gas is None:
        gas = create_gas(mechanism)

    _prepare_reactants(gas, equivalence_ratio, fuel, oxidizer)
    gas.TP = temperature, pressure
    gas.equilibrate("TP")

    return EquilibriumProperties(
        temperature=gas.T,
        pressure=gas.P,
        equivalence_ratio=equivalence_ratio,
        enthalpy=gas.enthalpy_mass,
        entropy=gas.entropy_mass,
        cp=gas.cp_mass,
        cv=gas.cv_mass,
        gamma=gas.cp_mass / gas.cv_mass,
        mean_molecular_weight=gas.mean_molecular_weight,
        density=gas.density,
        mole_fractions=dict(zip(gas.species_names, gas.X)),
    )


def _specific_volume(
    gas: ct.Solution,
    temperature: float,
    pressure: float,
    equivalence_ratio: float,
    fuel: str,
    oxidizer: Mapping[str, float],
) -> float:
    """Return the equilibrium specific volume at the given state, in m^3/kg."""
    _prepare_reactants(gas, equivalence_ratio, fuel, oxidizer)
    gas.TP = temperature, pressure
    gas.equilibrate("TP")
    return 1.0 / gas.density


def _specific_enthalpy(
    gas: ct.Solution,
    temperature: float,
    pressure: float,
    equivalence_ratio: float,
    fuel: str,
    oxidizer: Mapping[str, float],
) -> float:
    """Return the equilibrium specific enthalpy at the given state, in J/kg."""
    _prepare_reactants(gas, equivalence_ratio, fuel, oxidizer)
    gas.TP = temperature, pressure
    gas.equilibrate("TP")
    return gas.enthalpy_mass


def equilibrium_specific_heats(
    temperature: float,
    pressure: float,
    equivalence_ratio: float,
    fuel: str = DEFAULT_FUEL,
    oxidizer: Mapping[str, float] = DRY_AIR,
    gas: ct.Solution | None = None,
    mechanism: str = DEFAULT_MECHANISM,
    temperature_step: float = DEFAULT_TEMPERATURE_STEP,
    relative_pressure_step: float = DEFAULT_RELATIVE_PRESSURE_STEP,
) -> EquilibriumSpecificHeats:
    """Compute specific heats that account for shifting chemical equilibrium.

    The frozen specific heats returned by :func:`equilibrium_properties` hold
    the composition fixed. When a real mixture is heated, part of the energy
    instead breaks chemical bonds as the dissociation equilibrium shifts, so
    the effective specific heat is larger. Above roughly 2000 K the difference
    is substantial: at 2900 K and one atmosphere the equilibrium value exceeds
    the frozen one by more than a factor of three.

    Both quantities are physically meaningful limits. Frozen values apply when
    the flow is too fast for the chemistry to keep up, as in a turbine, while
    equilibrium values apply when the composition has time to adjust, as in a
    combustor. Real behaviour lies between them.

    Parameters
    ----------
    temperature : float
        Temperature in K.
    pressure : float
        Pressure in Pa.
    equivalence_ratio : float
        Fuel-to-air equivalence ratio.
    fuel : str, optional
        Fuel species name as defined in the mechanism.
    oxidizer : mapping of str to float, optional
        Oxidizer composition on a molar basis.
    gas : ct.Solution, optional
        Existing solution object to reuse. A new one is created when omitted.
    mechanism : str, optional
        Mechanism used when ``gas`` is not supplied.
    temperature_step : float, optional
        Half width of the temperature difference in K. Results are insensitive
        to this choice over at least two orders of magnitude.
    relative_pressure_step : float, optional
        Half width of the pressure difference, as a fraction of the pressure.

    Returns
    -------
    EquilibriumSpecificHeats
        The equilibrium specific heats and the isentropic exponent.

    Notes
    -----
    Derivatives are evaluated by central differences, which costs four
    additional equilibrium solves beyond the state itself. The specific heat at
    constant volume and the isentropic exponent follow from the standard
    thermodynamic relations used by NASA CEA:

    .. math::

        c_v = c_p + \\frac{p v}{T}
              \\frac{[(\\partial \\ln v / \\partial \\ln T)_p]^2}
                    {(\\partial \\ln v / \\partial \\ln p)_T}

    .. math::

        \\gamma_s = -\\frac{c_p / c_v}
                          {(\\partial \\ln v / \\partial \\ln p)_T}

    Because the isentropic exponent of a reacting mixture is not the ratio of
    its specific heats, it is reported separately rather than being derived
    from them.
    """
    if temperature <= temperature_step:
        raise ValueError(
            f"temperature must exceed the temperature step, got {temperature}"
        )
    if pressure <= 0.0:
        raise ValueError(f"pressure must be positive, got {pressure}")
    if not 0.0 < relative_pressure_step < 1.0:
        raise ValueError(
            "relative_pressure_step must lie between zero and one, got "
            f"{relative_pressure_step}"
        )

    if gas is None:
        gas = create_gas(mechanism)

    arguments = (equivalence_ratio, fuel, oxidizer)

    hot = temperature + temperature_step
    cold = temperature - temperature_step
    enthalpy_hot = _specific_enthalpy(gas, hot, pressure, *arguments)
    volume_hot = 1.0 / gas.density
    enthalpy_cold = _specific_enthalpy(gas, cold, pressure, *arguments)
    volume_cold = 1.0 / gas.density

    cp = (enthalpy_hot - enthalpy_cold) / (hot - cold)
    dlnv_dlnt = np.log(volume_hot / volume_cold) / np.log(hot / cold)

    pressure_step = pressure * relative_pressure_step
    high = pressure + pressure_step
    low = pressure - pressure_step
    volume_high = _specific_volume(gas, temperature, high, *arguments)
    volume_low = _specific_volume(gas, temperature, low, *arguments)
    dlnv_dlnp = np.log(volume_high / volume_low) / np.log(high / low)

    volume = _specific_volume(gas, temperature, pressure, *arguments)
    cv = cp + (pressure * volume / temperature) * dlnv_dlnt**2 / dlnv_dlnp
    isentropic_exponent = -(cp / cv) / dlnv_dlnp

    return EquilibriumSpecificHeats(
        cp=float(cp),
        cv=float(cv),
        isentropic_exponent=float(isentropic_exponent),
    )


def adiabatic_flame_temperature(
    inlet_temperature: float,
    pressure: float,
    equivalence_ratio: float,
    fuel: str = DEFAULT_FUEL,
    oxidizer: Mapping[str, float] = DRY_AIR,
    gas: ct.Solution | None = None,
    mechanism: str = DEFAULT_MECHANISM,
) -> float:
    """Compute the adiabatic flame temperature at constant pressure.

    The reactant mixture is equilibrated at constant enthalpy and pressure,
    which is the thermodynamic definition of adiabatic combustion.

    Parameters
    ----------
    inlet_temperature : float
        Temperature of the unburnt reactants in K.
    pressure : float
        Combustion pressure in Pa.
    equivalence_ratio : float
        Fuel-to-air equivalence ratio.
    fuel : str, optional
        Fuel species name as defined in the mechanism.
    oxidizer : mapping of str to float, optional
        Oxidizer composition on a molar basis.
    gas : ct.Solution, optional
        Existing solution object to reuse. A new one is created when omitted.
    mechanism : str, optional
        Mechanism used when ``gas`` is not supplied.

    Returns
    -------
    float
        Adiabatic flame temperature in K.

    Raises
    ------
    ValueError
        If ``inlet_temperature`` or ``pressure`` is not strictly positive.

    Examples
    --------
    >>> round(adiabatic_flame_temperature(298.15, 101325.0, 1.0))
    2387
    """
    if inlet_temperature <= 0.0:
        raise ValueError(
            f"inlet_temperature must be positive, got {inlet_temperature}"
        )
    if pressure <= 0.0:
        raise ValueError(f"pressure must be positive, got {pressure}")

    if gas is None:
        gas = create_gas(mechanism)

    _prepare_reactants(gas, equivalence_ratio, fuel, oxidizer)
    gas.TP = inlet_temperature, pressure
    gas.equilibrate("HP")

    return float(gas.T)
