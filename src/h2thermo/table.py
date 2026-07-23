"""Generation and storage of thermodynamic property tables.

A :class:`ThermoTable` holds equilibrium combustion product properties sampled
on a structured grid of temperature, pressure and equivalence ratio. Tables are
generated once by repeated calls to the equilibrium solver and then persisted,
so that downstream cycle analysis can retrieve properties without paying the
cost of a Gibbs minimisation at every operating point.

Equivalence ratio is used as the mixture coordinate throughout. Conversion to
the fuel-air ratio expected by engine cycle codes is the responsibility of the
export layer.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Mapping

import cantera as ct
import numpy as np

from h2thermo.equilibrium import (
    DEFAULT_FUEL,
    DEFAULT_MECHANISM,
    DRY_AIR,
    create_gas,
    equilibrium_properties,
    equilibrium_specific_heats,
)

__all__ = [
    "EQUILIBRIUM_PROPERTY_SOURCES",
    "GridSpecification",
    "PROPERTY_NAMES",
    "ThermoTable",
    "report_progress",
]

#: Scalar properties stored at every grid node, in the order used internally.
PROPERTY_NAMES: tuple[str, ...] = (
    "enthalpy",
    "entropy",
    "cp",
    "cv",
    "gamma",
    "cp_equilibrium",
    "cv_equilibrium",
    "isentropic_exponent",
    "mean_molecular_weight",
    "density",
)

#: Properties obtained from the shifting equilibrium calculation rather than
#: from the state itself, mapped to their attribute on the returned object.
EQUILIBRIUM_PROPERTY_SOURCES: dict[str, str] = {
    "cp_equilibrium": "cp",
    "cv_equilibrium": "cv",
    "isentropic_exponent": "isentropic_exponent",
}

#: Format version of the persisted file, incremented on breaking changes.
FILE_FORMAT_VERSION = 2


@dataclass(frozen=True)
class GridSpecification:
    """Structured sampling grid in temperature, pressure and equivalence ratio.

    Attributes
    ----------
    temperature : numpy.ndarray
        Strictly increasing, strictly positive temperature nodes in K.
    pressure : numpy.ndarray
        Strictly increasing, strictly positive pressure nodes in Pa.
    equivalence_ratio : numpy.ndarray
        Strictly increasing equivalence ratio nodes, dimensionless. May start
        at exactly zero, representing pure oxidizer with no fuel present;
        this is the row a full engine model needs for unburned sections such
        as an inlet or compressor.
    """

    temperature: np.ndarray
    pressure: np.ndarray
    equivalence_ratio: np.ndarray

    def __post_init__(self) -> None:
        for name in ("temperature", "pressure", "equivalence_ratio"):
            axis = np.asarray(getattr(self, name), dtype=float)
            object.__setattr__(self, name, axis)
            if axis.ndim != 1 or axis.size == 0:
                raise ValueError(f"{name} must be a non-empty one-dimensional array")
            if np.any(np.diff(axis) <= 0.0):
                raise ValueError(f"{name} values must be strictly increasing")
        if np.any(self.temperature <= 0.0):
            raise ValueError("temperature values must be strictly positive")
        if np.any(self.pressure <= 0.0):
            raise ValueError("pressure values must be strictly positive")
        if np.any(self.equivalence_ratio < 0.0):
            raise ValueError("equivalence_ratio values must be non-negative")

    @classmethod
    def linear(
        cls,
        temperature_range: tuple[float, float],
        pressure_range: tuple[float, float],
        equivalence_ratio_range: tuple[float, float],
        shape: tuple[int, int, int],
    ) -> "GridSpecification":
        """Build a grid with uniformly spaced nodes along each axis.

        Parameters
        ----------
        temperature_range : tuple of float
            Lower and upper temperature bounds in K.
        pressure_range : tuple of float
            Lower and upper pressure bounds in Pa.
        equivalence_ratio_range : tuple of float
            Lower and upper equivalence ratio bounds. The lower bound may be
            zero, for a pure-oxidizer (no fuel) row.
        shape : tuple of int
            Number of nodes along the temperature, pressure and equivalence
            ratio axes respectively.

        Returns
        -------
        GridSpecification
            The resulting grid.
        """
        n_temperature, n_pressure, n_equivalence_ratio = shape
        return cls(
            temperature=np.linspace(*temperature_range, n_temperature),
            pressure=np.linspace(*pressure_range, n_pressure),
            equivalence_ratio=np.linspace(
                *equivalence_ratio_range, n_equivalence_ratio
            ),
        )

    @property
    def shape(self) -> tuple[int, int, int]:
        """Number of nodes along each axis."""
        return (
            self.temperature.size,
            self.pressure.size,
            self.equivalence_ratio.size,
        )

    @property
    def size(self) -> int:
        """Total number of grid nodes."""
        return int(np.prod(self.shape))


@dataclass(frozen=True)
class ThermoTable:
    """Tabulated equilibrium properties of combustion products.

    Attributes
    ----------
    grid : GridSpecification
        The sampling grid the table was evaluated on.
    properties : dict of str to numpy.ndarray
        Scalar property arrays of shape ``grid.shape``, keyed by the names in
        :data:`PROPERTY_NAMES`. Nodes at which the solver failed hold NaN.
    species_names : tuple of str
        Species present in the mechanism, in mechanism order.
    mole_fractions : numpy.ndarray
        Equilibrium mole fractions of shape ``grid.shape + (n_species,)``.
    metadata : dict
        Provenance information such as fuel, oxidizer, mechanism and the time
        of generation.
    """

    grid: GridSpecification
    properties: dict[str, np.ndarray]
    species_names: tuple[str, ...]
    mole_fractions: np.ndarray
    metadata: dict

    @property
    def failed_node_count(self) -> int:
        """Number of grid nodes at which the equilibrium solver did not converge."""
        return int(np.count_nonzero(np.isnan(self.properties["cp"])))

    @classmethod
    def generate(
        cls,
        grid: GridSpecification,
        fuel: str = DEFAULT_FUEL,
        oxidizer: Mapping[str, float] = DRY_AIR,
        mechanism: str = DEFAULT_MECHANISM,
        progress: Callable[[int, int], None] | None = None,
    ) -> "ThermoTable":
        """Evaluate equilibrium properties at every node of ``grid``.

        Parameters
        ----------
        grid : GridSpecification
            Grid to sample.
        fuel : str, optional
            Fuel species name as defined in the mechanism.
        oxidizer : mapping of str to float, optional
            Oxidizer composition on a molar basis.
        mechanism : str, optional
            Reaction mechanism file.
        progress : callable, optional
            Called as ``progress(completed, total)`` after each node. Useful
            for reporting on large grids, which may take several minutes.

        Returns
        -------
        ThermoTable
            The populated table.

        Notes
        -----
        A single :class:`~cantera.Solution` object is reused across all nodes,
        because instantiating one is considerably more expensive than solving
        for equilibrium. Nodes at which the solver raises are recorded as NaN
        rather than aborting the run, so that a single pathological point does
        not discard hours of computation.
        """
        gas = create_gas(mechanism)
        species_names = tuple(gas.species_names)

        shape = grid.shape
        properties = {
            name: np.full(shape, np.nan, dtype=float) for name in PROPERTY_NAMES
        }
        mole_fractions = np.full(shape + (len(species_names),), np.nan, dtype=float)

        total = grid.size
        completed = 0

        for i, temperature in enumerate(grid.temperature):
            for j, pressure in enumerate(grid.pressure):
                for k, phi in enumerate(grid.equivalence_ratio):
                    try:
                        state = equilibrium_properties(
                            float(temperature),
                            float(pressure),
                            float(phi),
                            fuel=fuel,
                            oxidizer=oxidizer,
                            gas=gas,
                        )
                        shifting = equilibrium_specific_heats(
                            float(temperature),
                            float(pressure),
                            float(phi),
                            fuel=fuel,
                            oxidizer=oxidizer,
                            gas=gas,
                        )
                    except (ct.CanteraError, ValueError):
                        # Leave this node as NaN and continue; the count is
                        # reported through `failed_node_count`.
                        pass
                    else:
                        for name in PROPERTY_NAMES:
                            source = EQUILIBRIUM_PROPERTY_SOURCES.get(name)
                            properties[name][i, j, k] = (
                                getattr(shifting, source)
                                if source is not None
                                else getattr(state, name)
                            )
                        mole_fractions[i, j, k, :] = [
                            state.mole_fractions[species]
                            for species in species_names
                        ]

                    completed += 1
                    if progress is not None:
                        progress(completed, total)

        metadata = {
            "fuel": fuel,
            "oxidizer": dict(oxidizer),
            "mechanism": mechanism,
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "file_format_version": FILE_FORMAT_VERSION,
        }

        return cls(
            grid=grid,
            properties=properties,
            species_names=species_names,
            mole_fractions=mole_fractions,
            metadata=metadata,
        )

    def save(self, path: str | Path) -> Path:
        """Write the table to a compressed NumPy archive.

        Parameters
        ----------
        path : str or pathlib.Path
            Destination file. The ``.npz`` suffix is appended when absent.

        Returns
        -------
        pathlib.Path
            The path actually written.
        """
        path = Path(path)
        if path.suffix != ".npz":
            path = path.with_suffix(".npz")
        path.parent.mkdir(parents=True, exist_ok=True)

        arrays = {
            "axis_temperature": self.grid.temperature,
            "axis_pressure": self.grid.pressure,
            "axis_equivalence_ratio": self.grid.equivalence_ratio,
            "species_names": np.asarray(self.species_names),
            "mole_fractions": self.mole_fractions,
            "metadata": np.asarray(json.dumps(self.metadata)),
        }
        arrays.update(
            {f"property_{name}": array for name, array in self.properties.items()}
        )

        np.savez_compressed(path, **arrays)
        return path

    @classmethod
    def load(cls, path: str | Path) -> "ThermoTable":
        """Read a table previously written by :meth:`save`.

        Parameters
        ----------
        path : str or pathlib.Path
            Source file.

        Returns
        -------
        ThermoTable
            The reconstructed table.

        Raises
        ------
        ValueError
            If the file was written by an incompatible format version.
        """
        with np.load(Path(path), allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata"]))
            stored_version = metadata.get("file_format_version")
            if stored_version != FILE_FORMAT_VERSION:
                raise ValueError(
                    f"unsupported file format version {stored_version}, "
                    f"expected {FILE_FORMAT_VERSION}"
                )

            grid = GridSpecification(
                temperature=archive["axis_temperature"],
                pressure=archive["axis_pressure"],
                equivalence_ratio=archive["axis_equivalence_ratio"],
            )
            properties = {
                name: archive[f"property_{name}"] for name in PROPERTY_NAMES
            }
            species_names = tuple(str(name) for name in archive["species_names"])
            mole_fractions = archive["mole_fractions"]

        return cls(
            grid=grid,
            properties=properties,
            species_names=species_names,
            mole_fractions=mole_fractions,
            metadata=metadata,
        )

    def species_index(self, species: str) -> int:
        """Return the position of ``species`` within :attr:`species_names`.

        Parameters
        ----------
        species : str
            Species name as defined in the mechanism.

        Returns
        -------
        int
            Index into the last axis of :attr:`mole_fractions`.

        Raises
        ------
        KeyError
            If the species is not present in the mechanism.
        """
        try:
            return self.species_names.index(species)
        except ValueError as error:
            raise KeyError(
                f"species {species!r} is not present in the mechanism"
            ) from error

    def mole_fraction_of(self, species: str) -> np.ndarray:
        """Return the mole fraction field of a single species.

        Parameters
        ----------
        species : str
            Species name as defined in the mechanism.

        Returns
        -------
        numpy.ndarray
            Array of shape ``grid.shape``.
        """
        return self.mole_fractions[..., self.species_index(species)]


def report_progress(completed: int, total: int, step: int = 500) -> None:
    """Print a coarse progress line, intended as a default ``progress`` hook.

    Parameters
    ----------
    completed : int
        Number of nodes evaluated so far.
    total : int
        Total number of nodes.
    step : int, optional
        Reporting interval in nodes.
    """
    if completed % step == 0 or completed == total:
        percentage = 100.0 * completed / total
        print(f"  {completed:>7d} / {total} nodes ({percentage:5.1f} %)")
