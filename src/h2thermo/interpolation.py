"""Fast interpolated lookup of tabulated equilibrium properties.

A :class:`ThermoInterpolator` wraps a :class:`~h2thermo.table.ThermoTable` and
returns properties at arbitrary states within the tabulated envelope, without
solving for chemical equilibrium. This is the operation that makes the tables
useful to engine cycle codes: a single equilibrium solve costs on the order of
a millisecond, whereas an interpolated lookup costs microseconds.

Two choices in this module are driven by measured accuracy rather than
convention, and are documented in ``docs/validation.md``:

* Pressure is interpolated on a logarithmic coordinate, because equilibrium
  composition responds to pressure through the logarithm of the equilibrium
  constant.
* Density is not interpolated. It is recomputed from the interpolated mean
  molecular weight through the ideal gas equation of state, which is roughly
  fifty times more accurate than interpolating the density field directly.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import cantera as ct
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from h2thermo.table import ThermoTable

__all__ = ["INTERPOLATED_PROPERTIES", "InterpolatedState", "ThermoInterpolator"]

#: Properties obtained directly by interpolation. Density is excluded because
#: it is recovered more accurately from the equation of state.
INTERPOLATED_PROPERTIES: tuple[str, ...] = (
    "enthalpy",
    "entropy",
    "cp",
    "cv",
    "gamma",
    "cp_equilibrium",
    "cv_equilibrium",
    "isentropic_exponent",
    "mean_molecular_weight",
)


@dataclass(frozen=True)
class InterpolatedState:
    """Thermodynamic properties returned by an interpolated lookup.

    Each attribute is a scalar when the query was scalar, and an array of the
    broadcast query shape otherwise. Units follow the rest of the library:
    SI on a mass basis.

    Attributes
    ----------
    temperature : float or numpy.ndarray
        Temperature in K.
    pressure : float or numpy.ndarray
        Pressure in Pa.
    equivalence_ratio : float or numpy.ndarray
        Fuel-to-air equivalence ratio.
    enthalpy : float or numpy.ndarray
        Specific enthalpy in J/kg.
    entropy : float or numpy.ndarray
        Specific entropy in J/(kg K).
    cp : float or numpy.ndarray
        Frozen-composition specific heat at constant pressure in J/(kg K).
    cv : float or numpy.ndarray
        Frozen-composition specific heat at constant volume in J/(kg K).
    gamma : float or numpy.ndarray
        Frozen-composition ratio of specific heats.
    cp_equilibrium : float or numpy.ndarray
        Specific heat at constant pressure including the shifting equilibrium
        contribution, in J/(kg K).
    cv_equilibrium : float or numpy.ndarray
        Specific heat at constant volume including the shifting equilibrium
        contribution, in J/(kg K).
    isentropic_exponent : float or numpy.ndarray
        Exponent relating pressure and specific volume along an isentrope of
        the reacting mixture.
    mean_molecular_weight : float or numpy.ndarray
        Mean molecular weight in kg/kmol.
    density : float or numpy.ndarray
        Density in kg/m^3, recovered from the equation of state.
    """

    temperature: float | np.ndarray
    pressure: float | np.ndarray
    equivalence_ratio: float | np.ndarray
    enthalpy: float | np.ndarray
    entropy: float | np.ndarray
    cp: float | np.ndarray
    cv: float | np.ndarray
    gamma: float | np.ndarray
    cp_equilibrium: float | np.ndarray
    cv_equilibrium: float | np.ndarray
    isentropic_exponent: float | np.ndarray
    mean_molecular_weight: float | np.ndarray
    density: float | np.ndarray


class ThermoInterpolator:
    """Interpolated access to a tabulated set of equilibrium properties.

    Parameters
    ----------
    table : ThermoTable
        Table to interpolate. Every axis must carry at least two nodes, and the
        table must not contain non-convergent nodes.
    bounds_error : bool, optional
        When true, queries outside the tabulated envelope raise. When false,
        such queries return NaN. Silent extrapolation is never performed,
        because the underlying data is strongly non-linear outside the
        sampled region.

    Raises
    ------
    ValueError
        If the table has an axis with fewer than two nodes, or contains nodes
        at which the equilibrium solver did not converge.

    Examples
    --------
    >>> import numpy as np
    >>> from h2thermo import GridSpecification, ThermoTable
    >>> grid = GridSpecification.linear(
    ...     (800.0, 2800.0), (1.0e5, 40.0e5), (0.4, 1.0), (6, 4, 4)
    ... )
    >>> interpolator = ThermoInterpolator(ThermoTable.generate(grid))
    >>> state = interpolator.lookup(1500.0, 20.0e5, 0.6)
    >>> bool(1000.0 < state.cp < 2000.0)
    True
    """

    def __init__(self, table: ThermoTable, bounds_error: bool = True) -> None:
        if min(table.grid.shape) < 2:
            raise ValueError(
                "interpolation requires at least two nodes on every axis, "
                f"got shape {table.grid.shape}"
            )
        if table.failed_node_count:
            raise ValueError(
                f"table contains {table.failed_node_count} non-convergent "
                "nodes and cannot be interpolated"
            )

        self._table = table
        self._bounds_error = bounds_error
        self._axes = (
            table.grid.temperature,
            np.log(table.grid.pressure),
            table.grid.equivalence_ratio,
        )

        stacked = np.stack(
            [table.properties[name] for name in INTERPOLATED_PROPERTIES], axis=-1
        )
        self._stacked_properties = stacked
        self._property_interpolator = RegularGridInterpolator(
            self._axes,
            stacked,
            method="linear",
            bounds_error=bounds_error,
            fill_value=np.nan,
        )
        self._composition_interpolator = RegularGridInterpolator(
            self._axes,
            table.mole_fractions,
            method="linear",
            bounds_error=bounds_error,
            fill_value=np.nan,
        )

    @property
    def table(self) -> ThermoTable:
        """The table being interpolated."""
        return self._table

    @property
    def species_names(self) -> tuple[str, ...]:
        """Species carried by the underlying table."""
        return self._table.species_names

    def _query_points(
        self,
        temperature: float | Iterable[float],
        pressure: float | Iterable[float],
        equivalence_ratio: float | Iterable[float],
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, bool]:
        """Broadcast the query arguments into an array of grid coordinates."""
        temperature_array, pressure_array, phi_array = np.broadcast_arrays(
            np.asarray(temperature, dtype=float),
            np.asarray(pressure, dtype=float),
            np.asarray(equivalence_ratio, dtype=float),
        )
        is_scalar = temperature_array.ndim == 0

        if np.any(pressure_array <= 0.0):
            raise ValueError("pressure must be strictly positive")

        points = np.column_stack(
            [
                temperature_array.ravel(),
                np.log(pressure_array.ravel()),
                phi_array.ravel(),
            ]
        )
        return temperature_array, pressure_array, phi_array, points, is_scalar

    def _interpolate_scalar(
        self, temperature: float, log_pressure: float, equivalence_ratio: float
    ) -> np.ndarray | None:
        """Trilinearly interpolate the property vector at a single state.

        A dedicated scalar path exists because the general purpose interpolator
        carries enough call overhead that a single lookup would otherwise be
        slower than solving for equilibrium outright. Returns ``None`` when the
        query lies outside the envelope and out of range queries are permitted.
        """
        coordinates = (temperature, log_pressure, equivalence_ratio)
        indices = []
        weights = []
        for axis, value in zip(self._axes, coordinates):
            if not axis[0] <= value <= axis[-1]:
                if self._bounds_error:
                    raise ValueError(
                        f"query value {value:g} lies outside the tabulated "
                        f"range [{axis[0]:g}, {axis[-1]:g}]"
                    )
                return None
            lower = min(max(int(np.searchsorted(axis, value)) - 1, 0), axis.size - 2)
            indices.append(lower)
            weights.append(
                (value - axis[lower]) / (axis[lower + 1] - axis[lower])
            )

        i, j, k = indices
        u, v, w = weights
        cell = self._stacked_properties[i:i + 2, j:j + 2, k:k + 2, :]
        along_temperature = cell[0] * (1.0 - u) + cell[1] * u
        along_pressure = (
            along_temperature[0] * (1.0 - v) + along_temperature[1] * v
        )
        return along_pressure[0] * (1.0 - w) + along_pressure[1] * w

    def _lookup_scalar(
        self, temperature: float, pressure: float, equivalence_ratio: float
    ) -> InterpolatedState:
        """Lean lookup path for a single state, avoiding array machinery."""
        if pressure <= 0.0:
            raise ValueError("pressure must be strictly positive")

        values = self._interpolate_scalar(
            temperature, float(np.log(pressure)), equivalence_ratio
        )
        if values is None:
            values = np.full(len(INTERPOLATED_PROPERTIES), np.nan)

        properties = dict(zip(INTERPOLATED_PROPERTIES, (float(v) for v in values)))
        density = (
            pressure
            * properties["mean_molecular_weight"]
            / (ct.gas_constant * temperature)
        )
        return InterpolatedState(
            temperature=temperature,
            pressure=pressure,
            equivalence_ratio=equivalence_ratio,
            density=density,
            **properties,
        )

    def lookup(
        self,
        temperature: float | Iterable[float],
        pressure: float | Iterable[float],
        equivalence_ratio: float | Iterable[float],
    ) -> InterpolatedState:
        """Return interpolated properties at one or many states.

        Parameters
        ----------
        temperature : float or array_like
            Temperature in K.
        pressure : float or array_like
            Pressure in Pa.
        equivalence_ratio : float or array_like
            Fuel-to-air equivalence ratio.

        Returns
        -------
        InterpolatedState
            Properties at the requested states. Arguments are broadcast against
            one another, so scalars may be mixed with arrays.

        Raises
        ------
        ValueError
            If a query lies outside the tabulated envelope and the
            interpolator was constructed with ``bounds_error=True``.
        """
        if (
            np.isscalar(temperature)
            and np.isscalar(pressure)
            and np.isscalar(equivalence_ratio)
        ):
            return self._lookup_scalar(
                float(temperature), float(pressure), float(equivalence_ratio)
            )

        (
            temperature_array,
            pressure_array,
            phi_array,
            points,
            is_scalar,
        ) = self._query_points(temperature, pressure, equivalence_ratio)

        shape = temperature_array.shape
        if is_scalar:
            scalar_values = self._interpolate_scalar(
                float(temperature_array),
                float(points[0, 1]),
                float(phi_array),
            )
            values = (
                np.full((1, len(INTERPOLATED_PROPERTIES)), np.nan)
                if scalar_values is None
                else scalar_values.reshape(1, -1)
            )
        else:
            values = self._property_interpolator(points)
        properties = {
            name: values[:, index].reshape(shape)
            for index, name in enumerate(INTERPOLATED_PROPERTIES)
        }

        # The equation of state recovers density far more accurately than
        # interpolating the density field, which varies almost linearly with
        # pressure and is therefore poorly represented on a logarithmic axis.
        density = (
            pressure_array
            * properties["mean_molecular_weight"]
            / (ct.gas_constant * temperature_array)
        )

        def finish(array: np.ndarray) -> float | np.ndarray:
            return float(array) if is_scalar else array

        return InterpolatedState(
            temperature=finish(temperature_array),
            pressure=finish(pressure_array),
            equivalence_ratio=finish(phi_array),
            enthalpy=finish(properties["enthalpy"]),
            entropy=finish(properties["entropy"]),
            cp=finish(properties["cp"]),
            cv=finish(properties["cv"]),
            gamma=finish(properties["gamma"]),
            cp_equilibrium=finish(properties["cp_equilibrium"]),
            cv_equilibrium=finish(properties["cv_equilibrium"]),
            isentropic_exponent=finish(properties["isentropic_exponent"]),
            mean_molecular_weight=finish(properties["mean_molecular_weight"]),
            density=finish(density),
        )

    def mole_fractions(
        self,
        temperature: float | Iterable[float],
        pressure: float | Iterable[float],
        equivalence_ratio: float | Iterable[float],
    ) -> dict[str, float | np.ndarray]:
        """Return interpolated mole fractions of every species.

        Parameters
        ----------
        temperature : float or array_like
            Temperature in K.
        pressure : float or array_like
            Pressure in Pa.
        equivalence_ratio : float or array_like
            Fuel-to-air equivalence ratio.

        Returns
        -------
        dict of str to float or numpy.ndarray
            Mole fraction of each species in the mechanism.

        Notes
        -----
        Interpolated composition is considerably less accurate than the bulk
        properties, and the error grows as the mole fraction falls. Species
        present above a mole fraction of about 0.01 are reproduced to better
        than one per cent, whereas trace radicals below 0.0001 can deviate by
        tens of per cent. Call the equilibrium solver directly when trace
        composition matters.
        """
        temperature_array, _, _, points, is_scalar = self._query_points(
            temperature, pressure, equivalence_ratio
        )
        values = self._composition_interpolator(points)
        shape = temperature_array.shape
        return {
            species: (
                float(values[:, index].reshape(shape))
                if is_scalar
                else values[:, index].reshape(shape)
            )
            for index, species in enumerate(self.species_names)
        }
