# h2thermo

Equilibrium thermodynamic property tables for hydrogen-fuelled gas turbine
cycle analysis. The architecture is designed to extend to other alternative
fuels, including sustainable aviation fuels, ammonia and methane.

> **Status: early development.** The equilibrium core is implemented and
> validated at a single reference point. Table generation, interpolation and
> cycle-code export adapters are in progress. The public API is not yet stable.

## Motivation

Engine cycle analysis tools such as [pyCycle](https://github.com/OpenMDAO/pyCycle)
and T-MATS offer a fast tabular thermodynamics path, but the supplied tables are
prepared for conventional kerosene. Anyone modelling a hydrogen-burning engine
must generate an equivalent dataset themselves. This has so far been done
repeatedly as one-off work inside individual theses and papers, with no
reusable implementation published.

`h2thermo` aims to provide that missing component: a documented, tested and
openly licensed generator of combustion product properties.

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

| Input | Range |
| --- | --- |
| Fuel | Hydrogen (further fuels planned) |
| Temperature | 200 to 3000 K |
| Pressure | 1 to 60 bar |
| Equivalence ratio | 0.2 to 1.0 |

Outputs are the equilibrium composition and the corresponding specific
enthalpy, entropy, specific heats, ratio of specific heats, mean molecular
weight and density, all in SI units on a mass basis.

## Method

Equilibrium compositions are obtained by Gibbs free energy minimisation using
Cantera's solver, which follows the same thermodynamic principle as NASA CEA.
The `h2o2.yaml` mechanism distributed with Cantera is used; it is based on the
hydrogen kinetics of Burke et al. (2012) and includes the dissociation species
that matter at high temperature.

### Known limitations

Specific heats are currently frozen-composition values evaluated at the
equilibrium composition. They exclude the heat capacity contributed by shifting
dissociation equilibria. Below 2000 K the difference is under one per cent;
at 2900 K and 1 bar the equilibrium specific heat is more than three times the
frozen value. Shifting specific heats are the next priority.

The `h2o2.yaml` mechanism treats nitrogen as inert, so nitric oxide is absent.
The effect on bulk properties is below 0.04 per cent and therefore negligible
here, but the library is not suitable for emissions work as it stands.

Both limitations are quantified in [docs/validation.md](docs/validation.md).

## Validation

Properties are validated against NASA CEA over 140 reference states spanning
600 to 2900 K, 1 to 60 bar and equivalence ratios from 0.2 to 1.0. Mean
molecular weight and density agree to within 0.06 per cent, entropy to within
0.04 per cent and frozen specific heat to within 0.15 per cent. Internal
consistency checks on element conservation and the equation of state hold to
machine precision.

Interpolation adds well under 0.05 per cent on top of that, an order of
magnitude below the agreement with CEA. Full results, including the measured
cost of the two known limitations, are in
[docs/validation.md](docs/validation.md).

```bash
pytest
```

## Roadmap

1. Hydrogen-air property tables with validation against NASA CEA (complete)
2. Interpolated lookup (complete)
3. Shifting specific heats
4. Additional fuels: sustainable aviation fuels, ammonia, methane
5. Export adapters for pyCycle and T-MATS, plus generic CSV and JSON output
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
