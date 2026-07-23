"""Adapters that write ThermoTable data in formats external cycle codes read.

Each target format gets its own module, since the on-disk schema, axis
convention and property mapping are specific to the tool that reads them and
have to be established by inspecting or measuring against that tool rather
than assumed. See :mod:`h2thermo.export.pycycle` for the first of these.
"""

from __future__ import annotations

__all__: list[str] = []
