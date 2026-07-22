# Validation

This document records the checks that establish the accuracy of the
equilibrium properties produced by `h2thermo`, together with the quantitative
results. Every figure below is reproducible from the scripts and tests in this
repository.

Validation proceeds in three layers: internal consistency checks that require
no external reference, comparison against NASA CEA, and quantification of the
known limitations.

## 1. Internal consistency

These checks confirm that the tabulated data does not contradict itself. They
are cheap, require no external software, and run on every commit.

| Check | Maximum relative error |
| --- | --- |
| Element conservation (H, O and N atom ratios against the reactants) | 4.8 x 10<sup>-10</sup> |
| Ideal gas equation of state, rho = p M / (R T) | 5.6 x 10<sup>-10</sup> |
| Ratio of specific heats, gamma = cp / cv | 0 |
| Mole fractions summing to unity | < 10<sup>-9</sup> |

All four agree to machine precision, evaluated over a 50 x 20 x 20 grid
spanning 200 to 3000 K, 1 to 60 bar and equivalence ratios from 0.2 to 1.0.

## 2. Comparison against NASA CEA

NASA CEA is the reference implementation for equilibrium combustion
thermodynamics and is the origin of the thermodynamic data used by engine cycle
codes such as pyCycle. Agreement with CEA is therefore the central accuracy
claim of this project.

### Method

Reference states were generated with the NASA CEA Python package and are stored
in [`data/cea_reference_points.csv`](../data/cea_reference_points.csv). The
grid covers seven temperatures from 600 to 2900 K, four pressures from 1 to
60 bar and five equivalence ratios from 0.2 to 1.0, giving 140 states.

The CEA product species list was restricted to the nine species present in
Cantera's `h2o2.yaml` mechanism. Without this restriction, differences in the
species sets would be conflated with differences in the thermodynamic databases
and the solvers. The consequence of the restriction is quantified in section 3.

Storing the reference data rather than calling CEA from the test suite keeps
validation runnable in continuous integration without a compiled Fortran
dependency. Regenerate it with:

```bash
pip install cea
python scripts/generate_cea_reference.py
```

### Bulk properties

| Property | Maximum relative deviation | Mean |
| --- | --- | --- |
| Mean molecular weight | 0.052 % | 0.011 % |
| Density | 0.052 % | 0.012 % |
| Specific entropy | 0.032 % | 0.017 % |
| Frozen specific heat | 0.147 % | 0.060 % |
| Specific enthalpy | 11.4 kJ/kg absolute | 2.2 kJ/kg |

Enthalpy is reported as an absolute deviation because the absolute specific
enthalpy passes through zero within the tabulated range, which makes a relative
measure meaningless. For scale, 11.4 kJ/kg corresponds to roughly seven kelvin
of equivalent temperature error at a representative specific heat.

### Composition

| Species group | Maximum relative deviation | Mean |
| --- | --- | --- |
| Stable species (H2O, N2, O2, H2) | 3.0 % | 0.30 % |
| Radicals (OH, H, O) | 12.6 % | 4.8 % |

States with a mole fraction below 10<sup>-4</sup> are excluded, since their
relative deviation is dominated by round off and they carry no practical
weight.

Radicals deviate an order of magnitude more than stable species. This is
expected rather than alarming: radical concentrations depend exponentially on
Gibbs energies, so small differences between the two thermodynamic databases
are amplified. The stable species that determine the bulk properties agree far
more closely, which is consistent with the bulk property agreement above.

### Adiabatic flame temperature

Stoichiometric hydrogen-air ignited from ambient conditions at 1 atm reaches
2386.7 K, against an accepted literature value of approximately 2400 K.

## 3. Known limitations

### Frozen rather than shifting specific heats

The tabulated `cp` is a frozen-composition value evaluated at the equilibrium
composition. It excludes the heat capacity contributed by the shifting
dissociation equilibrium. CEA reports both quantities, which allows the
omission to be measured directly.

Ratio of equilibrium to frozen specific heat at stoichiometric conditions:

| Temperature | 1 bar | 5 bar | 20 bar | 60 bar |
| --- | --- | --- | --- | --- |
| 2600 K | 2.07 | 1.57 | 1.34 | 1.23 |
| 2900 K | 3.33 | 2.19 | 1.69 | 1.45 |

Below roughly 2000 K the two differ by less than one per cent, so the current
tables are appropriate there. Above that the difference cannot be ignored, and
users performing combustor exit calculations should be aware of it.

The pressure dependence is itself a useful check on the data. Dissociation
increases the number of moles, so raising the pressure suppresses it, and the
shifting contribution shrinks accordingly. The reference data reproduces this
trend, which the test suite asserts explicitly.

Adding shifting specific heats is the highest priority item on the roadmap.

### Absent nitrogen chemistry

Cantera's `h2o2.yaml` treats nitrogen as inert, whereas CEA includes nitric
oxide and related species. The consequence was measured by solving with CEA
twice, once with the matched species list and once with nitrogen chemistry
included, at an equivalence ratio of 0.6 and 20 bar:

| Temperature | NO mole fraction | Effect on mean molecular weight | Effect on frozen cp |
| --- | --- | --- | --- |
| 1500 K | 6.9 x 10<sup>-4</sup> | 0.0005 % | 0.000 % |
| 2000 K | 4.2 x 10<sup>-3</sup> | 0.0013 % | 0.002 % |
| 2400 K | 1.0 x 10<sup>-2</sup> | 0.0039 % | 0.012 % |
| 2800 K | 1.9 x 10<sup>-2</sup> | 0.0117 % | 0.035 % |

The thermodynamic consequence is negligible, two orders of magnitude below the
agreement already achieved on bulk properties. The omission is therefore
acceptable for property tables. It would not be acceptable for emissions work,
where nitric oxide is the quantity of interest, and a mechanism carrying
nitrogen chemistry would be required.

## 4. Interpolation accuracy

Interpolated lookup is what makes the tables useful, and its error is separate
from the accuracy of the underlying equilibrium data. The figures below were
measured on a 30 x 12 x 12 grid spanning 500 to 3000 K, 1 to 60 bar and
equivalence ratios from 0.2 to 1.0, by comparing interpolated values against
direct equilibrium solves at randomly chosen states inside the envelope.

### Bulk properties

| Property | Maximum relative error | Mean |
| --- | --- | --- |
| Ratio of specific heats | 0.016 % | 0.004 % |
| Frozen specific heat | 0.039 % | 0.007 % |
| Specific entropy | 0.046 % | 0.008 % |
| Mean molecular weight | 0.047 % | 0.005 % |
| Density | 0.043 % | 0.004 % |

Interpolation error is therefore an order of magnitude below the agreement
with CEA, so the total error remains dominated by the underlying
thermodynamic data rather than by tabulation.

### Two choices driven by measurement

**Density is derived, not interpolated.** Interpolating the density field gives
a maximum error of 2.21 per cent, because density varies almost linearly with
pressure and is poorly represented on a logarithmic axis. Recomputing it from
the interpolated mean molecular weight through the ideal gas equation of state
gives 0.043 per cent, an improvement of roughly fifty times.

**Pressure is interpolated logarithmically.** The benefit is modest for bulk
properties, reducing the maximum specific heat error from 0.022 to 0.015 per
cent on a 50 x 20 x 20 grid, but the coordinate is the physically natural one
and costs nothing.

### Composition

Interpolated mole fractions are far less accurate than bulk properties, and the
error grows as the species becomes rarer. For the hydroxyl radical:

| Mole fraction range | Maximum error | Mean |
| --- | --- | --- |
| above 0.01 | 0.53 % | 0.16 % |
| 0.001 to 0.01 | 9.4 % | 0.86 % |
| 0.0001 to 0.001 | 13.7 % | 2.4 % |
| below 0.0001 | 41.9 % | 8.5 % |

Species matter to the bulk properties in proportion to their abundance, so
this pattern is benign for property tables. It is not acceptable when trace
composition is itself the quantity of interest, and the equilibrium solver
should be called directly in that case.

### Grid resolution

Maximum specific heat interpolation error against grid size:

| Grid | Nodes | Maximum error |
| --- | --- | --- |
| 20 x 8 x 8 | 1280 | 0.097 % |
| 30 x 12 x 12 | 4320 | 0.040 % |
| 50 x 20 x 20 | 20000 | 0.015 % |

Even the coarsest grid stays well inside the accuracy of the underlying data,
so grid resolution can be chosen for file size and generation time rather than
for accuracy.

### Speed

Measured against a direct solve that reuses a single Cantera solution object,
which is the fastest way to call the solver:

| Operation | Time per state | Speed-up |
| --- | --- | --- |
| Equilibrium solve | 99 us | 1x |
| Scalar lookup | 34 us | 2.9x |
| Batched lookup | 0.87 us | 115x |

The scalar path is a hand written trilinear interpolation. Routing scalar
queries through the general purpose interpolator instead costs roughly 215 us,
which would make a single lookup slower than solving outright. The batched path
is where the benefit is decisive, so callers with many states should pass
arrays rather than looping.


## Reproducing these results

```bash
conda env create -f environment.yml
conda activate ct-env
pip install -e ".[dev]"
pytest
```

The validation suite compares every stored reference state against the library
and reports any deviation beyond the tolerances listed above.
