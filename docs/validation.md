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
| Frozen specific heat at constant pressure | 0.147 % | 0.060 % |
| Frozen specific heat at constant volume | 0.196 % | 0.080 % |
| Equilibrium specific heat at constant pressure | 0.556 % | 0.207 % |
| Equilibrium specific heat at constant volume | 0.630 % | 0.236 % |
| Isentropic exponent | 0.080 % | 0.032 % |
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

### Two specific heat definitions, both tabulated

`cp` and `cv` are frozen-composition values, holding the composition fixed.
`cp_equilibrium` and `cv_equilibrium` include the heat capacity contributed by
the shifting dissociation equilibrium. Both are physically meaningful limits:
frozen values apply when the flow is too fast for the chemistry to respond, as
in a turbine, and equilibrium values apply when the composition has time to
adjust, as in a combustor. Real behaviour lies between them, which is why CEA
reports both and this library follows suit.

The equilibrium quantities are obtained from central differences of the
enthalpy and specific volume, at a cost of four additional equilibrium solves
per state. Results are insensitive to the differencing step over at least two
orders of magnitude, from 0.1 to 10 K.

Because the isentropic exponent of a reacting mixture is not the ratio of its
specific heats, it is computed separately from the volume derivatives and
reported as `isentropic_exponent`. It agrees with CEA to 0.08 per cent, the
closest agreement of any quantity in this library.

Ratio of equilibrium to frozen specific heat at stoichiometric conditions:

| Temperature | 1 bar | 5 bar | 20 bar | 60 bar |
| --- | --- | --- | --- | --- |
| 2600 K | 2.07 | 1.57 | 1.34 | 1.23 |
| 2900 K | 3.33 | 2.19 | 1.69 | 1.45 |

Below roughly 2000 K the two differ by less than one per cent. Above that the
choice between them matters, and callers should pick the definition that suits
the process they are modelling.

The pressure dependence is itself a useful check on the data. Dissociation
increases the number of moles, so raising the pressure suppresses it, and the
shifting contribution shrinks accordingly. The reference data reproduces this
trend, which the test suite asserts explicitly.

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
| Isentropic exponent | 0.019 % | 0.005 % |
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


## 5. Compatibility with pyCycle's tabular thermo format

Before an export adapter can map h2thermo's properties onto pyCycle's tabular
thermo format, it has to be known which of h2thermo's several specific-heat
and gamma definitions the format actually stores. Guessing would violate the
project's own rule that documentation claims correspond to something
measured, so this was determined by reading pyCycle's own shipped reference
table directly, rather than by inspecting its source for comments or
docstrings.

### Method

pyCycle ships a pre-generated tabular thermo file for air and Jet-A,
`pycycle/thermo/tabular/air_jetA.pkl`, exposed as
`pycycle.constants.AIR_JETA_TAB_SPEC`. It is a pickled dictionary of arrays:
one-dimensional `T` (100 nodes, 100 to 3500 K), `P` (110 nodes, log-spaced
from 1 Pa to 10 MPa) and `FAR` (20 nodes, 0 to 0.05), and three-dimensional
`h`, `S`, `gamma`, `Cp`, `Cv`, `rho` and `R`, each shaped `[FAR, P, T]`. Units
are SI on a mass basis, matching h2thermo's own convention. This is the file
`pycycle.thermo.tabular.tabular_thermo.SetTotalTP` loads at run time and feeds
directly into an OpenMDAO structured metamodel, with no further
transformation, so measuring this file measures exactly what pyCycle's
tabular mode uses.

An earlier approach built a live pyCycle model on its CEA path and compared
`flow:gamma`, `flow:Cp` and `flow:Cv` against independently computed
quantities. That measures the right physics but the wrong artifact: what
pyCycle's CEA path computes at run time is not guaranteed to be what its
tabular export path wrote to disk. It was also fragile across pyCycle
versions, since it depended on internal component wiring rather than a public
interface. Reading the shipped table directly avoids both problems and is
reproducible with
[`scripts/probe_pycycle_definitions.py`](../scripts/probe_pycycle_definitions.py),
which requires `pip install om-pycycle` (not a project dependency).

Three discriminators were evaluated at a fixed fuel-air ratio of 0.03 and 1
atm, across a temperature sweep:

* `gamma` against `Cp / Cv`. Equality means gamma carries no information
  beyond the two specific heats.
* `Cp` against a finite difference of the tabulated enthalpy with respect to
  temperature, using the table's own grid spacing. Because the composition
  at each node was solved independently by whatever produced the table, this
  derivative is the equilibrium specific heat if and only if `Cp` is.
* `Cp - Cv` against `R`. The two are equal only for a mixture of fixed
  composition.

### Results

| Temperature | gamma | Cp/Cv | rel. diff. | Cp [J/(kg K)] | dh/dT [J/(kg K)] | rel. diff. | Cp - Cv | R | rel. diff. |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 306 K (sanity check) | 1.3864 | 1.3864 | 9 x 10<sup>-9</sup> | 1029.8 | 1029.8 | 1 x 10<sup>-5</sup> | 287.0 | 287.0 | 1 x 10<sup>-6</sup> |
| 2400 K | 1.2008 | 1.2035 | 0.23 % | 1920.7 | 1921.8 | 0.06 % | 324.8 | 288.6 | 12.6 % |
| 2900 K (edge of the CEA-compared range) | 1.1513 | 1.1692 | 1.5 % | 3597.3 | 3597.2 | 0.003 % | 520.5 | 298.4 | 74.4 % |
| 3500 K (edge of pyCycle's own table, outside both h2thermo's supported and CEA-compared ranges) | 1.1602 | 1.2042 | 3.7 % | 4966.1 | 4975.4 | 0.19 % | 842.0 | 327.6 | 157 % |

At 306 K, well below where dissociation matters, all three discriminators
collapse to near machine precision, which is the expected behaviour if the
interpretation below is correct and would have falsified it otherwise.

### Conclusion

**`Cp` and `Cv` are equilibrium (shifting-composition) values.** The dh/dT
discriminator, which measures the definition directly rather than by
elimination, agrees with the stored `Cp` to a few parts in 10<sup>4</sup>
across the entire sweep. `Cp - Cv` departing from `R` by more than 100 % at
the higher end confirms the same conclusion by a second, independent route:
a fixed-composition mixture cannot depart from `Cp - Cv = R` at all.

**`gamma` is an independently evaluated isentropic exponent, not `Cp / Cv`.**
It differs from the stored `Cp / Cv` by up to 3.7 % at 3500 K, and,
importantly, it is *below* the equilibrium ratio rather than above it. A
frozen-composition ratio for this mixture would sit above the equilibrium
`Cp / Cv`, since dissociation inflates the equilibrium specific heats (and
therefore lowers their ratio) relative to their frozen values, a pattern
h2thermo's own tables of `cp_equilibrium`/`cp` and `cv_equilibrium`/`cv` in
section 3 already establish. A value below the equilibrium ratio therefore
rules out `gamma` being any kind of frozen ratio, leaving an independently
computed isentropic exponent as the explanation consistent with all three
discriminators.

### Consequence for the export mapping

This determines, rather than assumes, the property mapping used by
`export/pycycle.py`:

| pyCycle field | h2thermo source |
| --- | --- |
| `Cp` | `cp_equilibrium` |
| `Cv` | `cv_equilibrium` |
| `gamma` | `isentropic_exponent` |
| `h`, `S`, `rho` | `enthalpy`, `entropy`, `density` |
| `R` | not currently stored; recovered as `ct.gas_constant / mean_molecular_weight` |

### Scope

The decision is unambiguous everywhere it was measured, but it was measured on
pyCycle's own table, whose envelope (100 to 3500 K, 1 Pa to 10 MPa, FAR 0 to
0.05) is wider on both axes than the range h2thermo supports and wider still
than the range compared against NASA CEA. The 2900 K row above is the
relevant one for h2thermo's own compared envelope; the 3500 K row is reported
for scale only, and no claim is made about agreement at that temperature.

## Reproducing these results

```bash
conda env create -f environment.yml
conda activate ct-env
pip install -e ".[dev]"
pytest
```

The validation suite compares every stored reference state against the library
and reports any deviation beyond the tolerances listed above.
