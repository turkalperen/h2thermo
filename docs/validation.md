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
| Element conservation (H, O and N atom ratios against the reactants) | 1.2 x 10<sup>-9</sup> |
| Ideal gas equation of state, rho = p M / (R T) | 5.6 x 10<sup>-10</sup> |
| Ratio of specific heats, gamma = cp / cv | 0 |
| Mole fractions summing to unity | < 10<sup>-9</sup> |

All four agree to machine precision, evaluated over a 50 x 20 x 20 grid
spanning 200 to 3000 K, 10 kPa to 60 bar and equivalence ratios from 0 to 1.0,
including the pure-air (zero equivalence ratio) row.

## 2. Comparison against NASA CEA

NASA CEA is the reference implementation for equilibrium combustion
thermodynamics and is the origin of the thermodynamic data used by engine cycle
codes such as pyCycle. Agreement with CEA is therefore the central accuracy
claim of this project.

### Method

Reference states were generated with the NASA CEA Python package and are stored
in [`data/cea_reference_points.csv`](../data/cea_reference_points.csv). The
grid covers eight temperatures from 300 to 2900 K, five pressures from 0.1 to
60 bar and five equivalence ratios from 0.2 to 1.0, giving 200 states. The
0.1 bar (10 kPa) tier is close to the ~22.6 kPa static pressure at 11 km
altitude, the low end of the range a compressor or combustor inlet would see
in practice.

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
| Density | 0.052 % | 0.011 % |
| Specific entropy | 0.032 % | 0.017 % |
| Frozen specific heat at constant pressure | 0.147 % | 0.067 % |
| Frozen specific heat at constant volume | 0.196 % | 0.091 % |
| Equilibrium specific heat at constant pressure | 1.176 % | 0.221 % |
| Equilibrium specific heat at constant volume | 1.167 % | 0.253 % |
| Isentropic exponent | 0.084 % | 0.034 % |
| Specific enthalpy | 11.3 kJ/kg absolute | 2.1 kJ/kg |

Enthalpy is reported as an absolute deviation because the absolute specific
enthalpy passes through zero within the tabulated range, which makes a relative
measure meaningless. For scale, 11.3 kJ/kg corresponds to roughly seven kelvin
of equivalent temperature error at a representative specific heat.

The equilibrium specific heats agree more than twice as loosely as they did
before the 0.1 bar tier was added: the worst points are all at 0.1 bar and
2200-2900 K, where dissociation, and therefore the shifting-composition
contribution to `cp`/`cv`, is largest. This is the same effect the ratio table
in section 3 quantifies directly, not a new source of error: the two
thermodynamic databases agree closely on the underlying frozen properties, so
their disagreement is amplified wherever the equilibrium and frozen values
diverge most from each other. The tolerance in the test suite already carried
enough margin to absorb this without a change.

### Composition

| Species group | Maximum relative deviation | Mean |
| --- | --- | --- |
| Stable species (H2O, N2, O2, H2) | 3.0 % | 0.29 % |
| Radicals (OH, H, O) | 12.6 % | 4.7 % |

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

| Temperature | 0.1 bar | 1 bar | 5 bar | 20 bar | 60 bar |
| --- | --- | --- | --- | --- | --- |
| 2600 K | 3.85 | 2.07 | 1.57 | 1.34 | 1.23 |
| 2900 K | 7.13 | 3.33 | 2.19 | 1.69 | 1.45 |

At 0.1 bar, the altitude-relevant pressure added for this range, the
equilibrium value reaches more than seven times the frozen one at 2900 K.
This is also why the CEA agreement on the equilibrium specific heats in
section 2 loosens at that corner: the larger the shifting contribution, the
more the result depends on radical thermodynamics the two databases
reproduce least closely.

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

## 6. Cross-validation of the pure-air state against pyCycle

Section 2 validates hydrogen-air combustion products against NASA CEA. The
pure-air state added in PR-6 (`equivalence_ratio=0.0`) burns no fuel, so it
is not covered by that reference set. This section checks it a different
way: directly against pyCycle's own shipped `FAR = 0` row in
`pycycle.constants.AIR_JETA_TAB_SPEC`, the same file used in section 5. This
is a second channel entirely independent of CEA, at the cost of comparing
against a table whose own provenance is unknown rather than a named
reference implementation.

### Method

At each probed grid node, h2thermo's `equilibrium_properties` and
`equilibrium_specific_heats` are evaluated at `equivalence_ratio=0.0` and
mapped to pyCycle's field names using the same convention established in
section 5 (`Cp`/`Cv` as the equilibrium specific heats, `gamma` as the
isentropic exponent, `R` recovered from mean molecular weight). Where the two
disagree, NASA CEA is brought in as a third, independent opinion rather than
leaving the disagreement unresolved; CEA's solver hits a singular update
matrix at exactly zero fuel, so a proxy equivalence ratio of 10<sup>-4</sup>
stands in for pure air, checked for convergence down to 10<sup>-6</sup> with
the result stable to 0.02%. Reproducible with
[`scripts/compare_pure_air_to_pycycle.py`](../scripts/compare_pure_air_to_pycycle.py),
which requires `pip install om-pycycle` and, for the tie-breaker,
`pip install cea`.

### Realistic operating range: close agreement

An unburned-air row exists for sections that never see fuel: an inlet, a
fan, a compressor, cooling bypass air. None of those reach anywhere near
2900 K; a high-pressure-ratio compressor discharge is a more representative
upper bound. Swept at 5 bar:

| Temperature | Cp deviation | Gamma deviation |
| --- | --- | --- |
| 300 K | 0.55 % | 0.06 % |
| 600 K | 0.58 % | 0.07 % |
| 900 K | 0.75 % | 0.12 % |
| 1200 K | 0.27 % | 0.04 % |

Across the range a real air-breathing section actually operates in, h2thermo
agrees with pyCycle's own reference data about as closely as section 2 shows
it agreeing with CEA on combustion products.

### Bulk properties across the full pressure range, at 2900 K

Pushed to the hottest node shared by both tables, to see where agreement
degrades:

| Property | 0.1 bar | 1 bar | 5 bar | 20 bar | 60 bar |
| --- | --- | --- | --- | --- | --- |
| Density | 0.65 % | 0.48 % | 0.44 % | 0.42 % | 0.41 % |
| Specific gas constant | 0.65 % | 0.48 % | 0.44 % | 0.42 % | 0.42 % |
| Specific entropy | 0.12 % | 0.07 % | 0.14 % | 0.19 % | 0.22 % |
| Specific enthalpy | 0.70 % | 2.13 % | 2.60 % | 2.80 % | 2.90 % |
| Isentropic exponent | 0.09 % | 0.72 % | 1.39 % | 1.92 % | 2.28 % |

Density and R sit at a roughly constant ~0.4-0.65%, and the cause is
identifiable: `DRY_AIR`, h2thermo's default oxidizer, is a simplified
two-component (O2, N2) mixture that omits argon, at 0.93% of real air by
mole. At 300 K and 1 bar, a low-dissociation state where the two tables
should already agree closely on frozen properties, pyCycle's implied mean
molecular weight is 28.965 kg/kmol, matching the handbook value for real dry
air (28.9647) to four figures, against 28.851 kg/kmol for h2thermo's
argon-free `DRY_AIR`. Substituting a realistic argon-inclusive oxidizer
closes the gap: the Cp deviation at that state falls from 0.55% to 0.09%,
confirmed with `check_argon_hypothesis()` in the comparison script.

### Resolved: pyCycle's table, not h2thermo, is the less accurate one above 1500 K

The equilibrium `Cp`/`Cv` deviation between h2thermo and pyCycle is small
(under 1%) up to about 1200 K and grows to 5-10% by 2200-2900 K. Two
explanations were checked and ruled out rather than assumed. The argon
composition difference does not explain it: substituting a realistic
argon-inclusive oxidizer at the hottest, highest-pressure node makes the
deviation slightly worse (8.00% to 8.69%), not better. Dissociation strength
does not explain it either: at 100 bar, the highest pressure pyCycle's table
stores, O2 dissociation is almost fully suppressed (h2thermo's own
atomic-oxygen mole fraction there is 0.0039), yet pyCycle's `Cp` still
exceeds h2thermo's frozen `cp` by 19.5%. A gap that large with essentially no
dissociation on either side rules out a shifting-composition artifact.

With both plausible causes eliminated, the disagreement was resolved rather
than left open, by bringing in NASA CEA as a third, independent opinion at
every state in both sweeps:

| State | cp h2thermo | cp pyCycle | cp CEA | \|h2thermo - CEA\| | \|pyCycle - CEA\| |
| --- | --- | --- | --- | --- | --- |
| 2900 K, 0.1 bar | 4404.05 | 4242.62 | 4390.09 | 0.32 % | 3.36 % |
| 2900 K, 1 bar | 2453.85 | 2488.46 | 2447.01 | 0.28 % | 1.69 % |
| 2900 K, 5 bar | 1835.19 | 1927.86 | 1831.81 | 0.18 % | 5.24 % |
| 2900 K, 20 bar | 1580.66 | 1696.43 | 1578.83 | 0.12 % | 7.45 % |
| 2900 K, 60 bar | 1457.52 | 1584.27 | 1456.49 | 0.07 % | 8.77 % |
| 300-2900 K, 5 bar (9 points) | -- | -- | -- | 0.02-0.18 % | 0.34-7.23 % |

Across all 14 states where CEA converged, h2thermo's mean deviation from CEA
is 0.13%, essentially the same level of agreement section 2 already
established for combustion products. pyCycle's mean deviation from the same
reference is 3.62%, and reaches 8.8% at the hottest, highest-pressure node.
h2thermo is not the source of the disagreement documented above; pyCycle's
own shipped air table is the less accurate of the two at this corner, most
likely because it was built and tuned around realistic combustion
conditions and was never validated at temperatures no unburned air actually
reaches.

This does not change the practical conclusion. No inlet, fan or compressor
carries unburned air anywhere close to 1500 K in an actual engine; air only
reaches that temperature after combustion, at which point
`equivalence_ratio` is no longer zero and the state is already covered by
section 2's CEA comparison. What changes is which library a caller should
trust if they do probe that corner: h2thermo, not pyCycle's own reference
table.

## Reproducing these results

```bash
conda env create -f environment.yml
conda activate ct-env
pip install -e ".[dev]"
pytest
```

The validation suite compares every stored reference state against the library
and reports any deviation beyond the tolerances listed above.
