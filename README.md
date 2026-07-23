# h2thermo

[![tests](https://github.com/turkalperen/h2thermo/actions/workflows/tests.yml/badge.svg)](https://github.com/turkalperen/h2thermo/actions/workflows/tests.yml)

Equilibrium thermodynamic property tables for hydrogen-fuelled gas turbine
cycle analysis. The architecture is designed to extend to other alternative
fuels, including sustainable aviation fuels, ammonia and methane.

> **Status: early development.** Table generation, interpolated lookup,
> validation against NASA CEA and a pyCycle export adapter are complete for
> hydrogen-air. The exported table does not yet include a pure-air (zero
> equivalence ratio) row; T-MATS export, further fuels and a stable public
> API are not yet implemented.

## Motivation

Engine cycle analysis tools such as [pyCycle](https://github.com/OpenMDAO/pyCycle)
and T-MATS offer a fast tabular thermodynamics path, but the tables they ship
with are prepared for conventional kerosene. Anyone modelling a
hydrogen-burning engine has to generate an equivalent dataset first. Both tools
provide example scripts for doing so, so the work is not blocked; it is instead
repeated from scratch in each thesis and paper, and the resulting tables are
generally used without a published statement of how far they depart from a
reference solution.

`h2thermo` aims to close that gap: a documented, openly licensed generator whose
output is validated against NASA CEA across the operating envelope, with the
residual error measured rather than assumed.

## Installation

The package depends on [Cantera](https://cantera.org), which is most reliably
installed through Conda:

```bash
conda env create -f environment.yml
conda activate ct-env
pip install -e ".[dev]"
```

## Quick start

```python
from h2thermo import adiabatic_flame_temperature, equilibrium_properties

# Adiabatic flame temperature of stoichiometric hydrogen-air at 1 atm.
flame_temperature = adiabatic_flame_temperature(
    inlet_temperature=298.15, pressure=101325.0, equivalence_ratio=1.0
)
print(f"{flame_temperature:.1f} K")

# Product properties at a prescribed combustor state.
state = equilibrium_properties(
    temperature=1800.0, pressure=20.0e5, equivalence_ratio=0.6
)
print(f"cp = {state.cp:.1f} J/(kg K), gamma = {state.gamma:.4f}")
```

Generating a table and querying it:

```python
import numpy as np

from h2thermo import GridSpecification, ThermoInterpolator, ThermoTable

grid = GridSpecification(
    temperature=np.linspace(500.0, 3000.0, 30),
    pressure=np.geomspace(1.0e5, 60.0e5, 12),
    equivalence_ratio=np.linspace(0.2, 1.0, 12),
)
table = ThermoTable.generate(grid)
table.save("data/generated/h2_air_table.npz")

interpolator = ThermoInterpolator(table)

# A single state, or many at once.
state = interpolator.lookup(1837.0, 13.7e5, 0.63)
states = interpolator.lookup(
    np.linspace(1000.0, 2500.0, 1000), 20.0e5, 0.6
)
```

## Scope

| Input | Supported | Compared against CEA |
| --- | --- | --- |
| Fuel | Hydrogen (further fuels planned) | Hydrogen |
| Temperature | 200 to 3000 K | 600 to 2900 K |
| Pressure | 1 to 60 bar | 1 to 60 bar |
| Equivalence ratio | 0.2 to 1.0 | 0.2 to 1.0 |

The solver returns results across the whole supported envelope. Outside the
compared range the agreement with CEA reported below has not been established,
so the temperature extremes should be treated as unverified.

Outputs are the equilibrium composition and the corresponding specific
enthalpy, entropy, frozen and equilibrium specific heats, ratio of specific
heats, isentropic exponent, mean molecular weight and density, all in SI units
on a mass basis.

## Method

Equilibrium compositions are obtained by Gibbs free energy minimisation using
Cantera's solver, which follows the same thermodynamic principle as NASA CEA.
The `h2o2.yaml` mechanism distributed with Cantera is used; it is based on the
hydrogen kinetics of Burke et al. (2012) and includes the dissociation species
that matter at high temperature.

### Frozen and equilibrium specific heats

Both definitions are provided. `cp` and `cv` hold the composition fixed, while
`cp_equilibrium` and `cv_equilibrium` include the heat capacity contributed by
shifting dissociation. Below 2000 K they agree to within one per cent; at
2900 K and 1 bar the equilibrium value is more than three times the frozen one.
Frozen values suit fast flows such as a turbine, equilibrium values suit a
combustor, and real behaviour lies between them.

Since the isentropic exponent of a reacting mixture is not the ratio of its
specific heats, it is reported separately as `isentropic_exponent`.

### Known limitation

`GridSpecification` requires a strictly positive equivalence ratio, so an
exported table never has pyCycle's `FAR = 0` (pure air) row. A full engine
model that includes unburned sections such as an inlet or compressor needs
that row. See [`h2thermo.export.pycycle`](src/h2thermo/export/pycycle.py) for
details.

## Validation

Properties are validated against NASA CEA over 140 reference states spanning
600 to 2900 K, 1 to 60 bar and equivalence ratios from 0.2 to 1.0. Mean
molecular weight and density agree to within 0.06 per cent, entropy to within
0.04 per cent and frozen specific heat to within 0.15 per cent. Internal
consistency checks on element conservation and the equation of state hold to
machine precision.

Interpolation adds well under 0.05 per cent on top of that, an order of
magnitude below the agreement with CEA. Full results, including the measured
effect of the inert-nitrogen treatment, are in [docs/validation.md](docs/validation.md).

```bash
pytest
```

## Roadmap

1. Hydrogen-air property tables with validation against NASA CEA (complete)
2. Interpolated lookup (complete)
3. Frozen and equilibrium specific heats (complete)
4. Export adapter for pyCycle (complete); T-MATS and generic CSV/JSON output remaining
5. Additional fuels: sustainable aviation fuels, ammonia, methane
6. Packaging and documentation

## References

- Burke, M. P., Chaos, M., Ju, Y., Dryer, F. L., Klippenstein, S. J. (2012).
  Comprehensive H2/O2 kinetic model for high-pressure combustion.
  *International Journal of Chemical Kinetics*, 44(7), 444-474.
- Goodwin, D. G., Moffat, H. K., Speth, R. L. et al. Cantera: An object-oriented
  software toolkit for chemical kinetics, thermodynamics, and transport
  processes. https://cantera.org

## License

Released under the MIT License. See [LICENSE](LICENSE).
