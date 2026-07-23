"""Export a ThermoTable in the tabular thermo format read by pyCycle.

pyCycle's tabular thermo mode (``pycycle.thermo.tabular.tabular_thermo``)
reads a pickled dictionary of NumPy arrays: one-dimensional ``T``, ``P`` and
``FAR`` axes, and ``h``, ``S``, ``gamma``, ``Cp``, ``Cv``, ``rho`` and ``R``
each shaped ``[FAR, P, T]``. This was determined by reading pyCycle's own
shipped reference table (``pycycle.constants.AIR_JETA_TAB_SPEC``) rather than
assumed; see ``scripts/probe_pycycle_definitions.py`` and section 5 of
``docs/validation.md``.

That same measurement determined which of h2thermo's several specific-heat
and gamma definitions pyCycle's format actually stores: ``Cp`` and ``Cv`` are
equilibrium (shifting-composition) values, and ``gamma`` is an independently
evaluated isentropic exponent rather than their ratio. The mapping below
follows directly from that result.

Two conversions happen that are not simple renames:

* pyCycle indexes composition by fuel-air ratio (FAR), h2thermo by
  equivalence ratio. The two are related by a fixed multiplicative constant,
  ``FAR = equivalence_ratio * FAR_stoichiometric``, so relabelling the axis
  introduces no interpolation error. ``FAR_stoichiometric`` is computed from
  the mechanism rather than hard-coded, so this works for whichever fuel a
  table was generated with.
* pyCycle's format stores a specific gas constant ``R``, which h2thermo does
  not tabulate directly. It is recovered from the interpolated mean
  molecular weight through ``R = R_universal / M``, the same relation
  :mod:`h2thermo.interpolation` uses to recover density.

:class:`~h2thermo.table.GridSpecification` accepts an equivalence ratio of
zero, representing pure oxidizer with no fuel present. Including 0.0 as the
lowest node of the equivalence-ratio axis before calling
:meth:`~h2thermo.table.ThermoTable.generate` produces a ``FAR = 0`` row on
export, matching pyCycle's own air/Jet-A table and giving a full engine
model the pure-air state it needs for unburned sections such as an inlet or
compressor. A table generated without 0.0 on that axis will not have this
row; this is a property of the table passed in; :func:`write_pycycle_table`
does not add or require it.

pyCycle's tabular thermo format is a pickle file, which is not a safe format
to load from an untrusted source. Writing one is required for compatibility
with pyCycle, and the risk is pyCycle's to manage on load, not introduced by
writing it here.
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Mapping

import cantera as ct
import numpy as np

from h2thermo.equilibrium import DEFAULT_FUEL, DEFAULT_MECHANISM, DRY_AIR, create_gas
from h2thermo.table import ThermoTable

__all__ = [
    "PYCYCLE_PROPERTY_SOURCES",
    "stoichiometric_fuel_air_ratio",
    "write_pycycle_table",
]

#: pyCycle field name mapped to the h2thermo property it is filled from.
#: Determined by measurement rather than assumption; see the module
#: docstring and docs/validation.md section 5.
PYCYCLE_PROPERTY_SOURCES: dict[str, str] = {
    "h": "enthalpy",
    "S": "entropy",
    "Cp": "cp_equilibrium",
    "Cv": "cv_equilibrium",
    "gamma": "isentropic_exponent",
    "rho": "density",
}


def stoichiometric_fuel_air_ratio(
    fuel: str = DEFAULT_FUEL,
    oxidizer: Mapping[str, float] = DRY_AIR,
    mechanism: str = DEFAULT_MECHANISM,
) -> float:
    """Return the stoichiometric fuel-to-air mass ratio for a fuel/oxidizer pair.

    Parameters
    ----------
    fuel : str, optional
        Fuel species name as defined in the mechanism.
    oxidizer : mapping of str to float, optional
        Oxidizer composition on a molar basis.
    mechanism : str, optional
        Reaction mechanism file.

    Returns
    -------
    float
        Mass of fuel per unit mass of oxidizer at an equivalence ratio of
        one. Fuel-air ratio and equivalence ratio are related by
        ``FAR = equivalence_ratio * stoichiometric_fuel_air_ratio(...)``.

    Notes
    -----
    Computed from the mechanism rather than hard-coded, so it stays correct
    for whichever fuel a table was generated with.
    """
    gas = create_gas(mechanism)
    gas.set_equivalence_ratio(1.0, fuel, dict(oxidizer))
    fuel_mass_fraction = float(gas.Y[gas.species_index(fuel)])
    return fuel_mass_fraction / (1.0 - fuel_mass_fraction)


def _reorder_to_far_pressure_temperature(array: np.ndarray) -> np.ndarray:
    """Transpose a ``[T, P, FAR]``-shaped array to pyCycle's ``[FAR, P, T]``."""
    return np.transpose(array, axes=(2, 1, 0))


def write_pycycle_table(table: ThermoTable, path: str | Path) -> Path:
    """Write ``table`` in the pickle format read by pyCycle's tabular thermo.

    Parameters
    ----------
    table : ThermoTable
        Table to export. Fuel, oxidizer and mechanism are taken from
        ``table.metadata``.
    path : str or pathlib.Path
        Destination file. The ``.pkl`` suffix is appended when absent.

    Returns
    -------
    pathlib.Path
        The path actually written.

    Raises
    ------
    ValueError
        If the table contains any non-convergent (NaN) node. pyCycle's
        structured metamodel requires a complete grid.
    """
    if table.failed_node_count:
        raise ValueError(
            f"table contains {table.failed_node_count} non-convergent nodes "
            "and cannot be exported"
        )

    far_stoichiometric = stoichiometric_fuel_air_ratio(
        fuel=table.metadata["fuel"],
        oxidizer=table.metadata["oxidizer"],
        mechanism=table.metadata["mechanism"],
    )
    fuel_air_ratio = table.grid.equivalence_ratio * far_stoichiometric

    mean_molecular_weight = _reorder_to_far_pressure_temperature(
        table.properties["mean_molecular_weight"]
    )
    specific_gas_constant = ct.gas_constant / mean_molecular_weight

    payload = {
        "T": table.grid.temperature,
        "P": table.grid.pressure,
        "FAR": fuel_air_ratio,
        "R": specific_gas_constant,
    }
    for pycycle_name, h2thermo_name in PYCYCLE_PROPERTY_SOURCES.items():
        payload[pycycle_name] = _reorder_to_far_pressure_temperature(
            table.properties[h2thermo_name]
        )

    path = Path(path)
    if path.suffix != ".pkl":
        path = path.with_suffix(".pkl")
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("wb") as handle:
        pickle.dump(payload, handle)

    return path
