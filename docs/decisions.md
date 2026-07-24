# PR-8

Cross-validates the pure-air state added in PR-6 against pyCycle's own
`FAR = 0` row, a channel independent of the hydrogen-combustion NASA CEA
reference set in section 2, which has nothing to say about a fuel-free
mixture. Where h2thermo and pyCycle disagreed, NASA CEA was brought back in
as a third opinion to settle which one was actually accurate, rather than
leaving two disagreeing tables undifferentiated.

## What is included

- `scripts/compare_pure_air_to_pycycle.py`: reads
  `pycycle.constants.AIR_JETA_TAB_SPEC` directly (the same artifact
  `scripts/probe_pycycle_definitions.py` uses) and compares it against
  h2thermo at matching grid nodes, using the property mapping established
  in PR-5. Includes a temperature sweep at a realistic pressure, a pressure
  sweep at the hottest shared temperature, a targeted check of whether
  `DRY_AIR`'s missing argon explains the deviations found, and an optional
  NASA CEA tie-breaker (`pip install cea`) run at every probed state.
- `docs/validation.md` section 6 and a pointer in `README.md`'s Validation
  section, with the results below.

## Result: good agreement where the row is actually used

Below about 1200 K, at 5 bar, Cp agrees with pyCycle's own table to within
0.3-0.9% and gamma to within 0.04-0.12% -- comparable to the CEA agreement
on combustion products in section 2. This is the range that matters: an
inlet, fan or compressor never carries unburned air anywhere near 2900 K.

## Above ~1500 K: two hypotheses ruled out, then resolved with a third opinion

The Cp/Cv deviation between h2thermo and pyCycle grows above roughly 1500 K,
reaching 5-10% by 2200-2900 K. Two explanations were checked and ruled out:

- **Argon.** `DRY_AIR` omits it; real air is 0.93% argon by mole. At 300 K,
  1 bar this fully explains a smaller deviation seen there: pyCycle's
  implied molecular weight (28.965 kg/kmol) matches the handbook value for
  real dry air almost exactly, and substituting a realistic argon-inclusive
  oxidizer drops the Cp deviation from 0.55% to 0.09%. At the hot,
  high-pressure corner it does not help -- the deviation gets marginally
  worse (8.00% to 8.69%) with argon included.
- **Dissociation strength.** At 100 bar, the highest pressure pyCycle's
  table stores, O2 dissociation is nearly fully suppressed (h2thermo's own
  atomic-oxygen mole fraction is 0.0039), yet pyCycle's Cp still exceeds
  h2thermo's frozen cp by 19.5%. A gap that large with essentially no
  dissociation on either side cannot be a dissociation-modelling artifact.

With neither hypothesis holding, the question of which table is actually
right was answerable, not just askable: CEA can solve pure air too, given a
proxy equivalence ratio of 10<sup>-4</sup> to route around a singular update
matrix at exactly zero fuel (stable to 0.02% down to 10<sup>-6</sup>, checked
before trusting it). Run at all 14 states in both sweeps, h2thermo's mean
deviation from CEA is 0.13%; pyCycle's is 3.62%, reaching 8.8% at the
hottest, highest-pressure node. **h2thermo is not the source of the
disagreement; pyCycle's own shipped air table is the less accurate one at
this corner.**

## Why the practical conclusion is unchanged even though the diagnosis is

No inlet, fan or compressor carries unburned air anywhere close to 1500 K in
an actual engine; air only reaches that temperature after combustion, at
which point equivalence_ratio is no longer zero and the state is already
covered by section 2's CEA comparison. The high-temperature pure-air corner
remains a mathematically valid but physically unreachable input. What the
CEA tie-breaker changes is not whether this matters operationally, but who
should be trusted if someone does probe that corner: h2thermo, not pyCycle's
own reference table.

## Test plan

- [x] `pip install om-pycycle && python scripts/compare_pure_air_to_pycycle.py`
      -- reproduces every figure above
- [x] `pip install cea` additionally installed -- tie-breaker section runs
      and reproduces the 0.13%/3.62% means
- [x] `flake8 src tests scripts examples` -- clean

# PR-7

Extends the NASA CEA validation envelope down to 10 kPa (0.1 bar) and 300 K,
closing the gap left after PR-5: the previously validated envelope (1-60 bar,
600-2900 K) sat entirely above the ~22.6 kPa static pressure at 11 km
altitude that a compressor or inlet section would actually see, and left the
library's own claimed 200 K lower temperature bound entirely unvalidated
below 600 K.

## What is included

- `scripts/generate_cea_reference.py`: added a 0.1 bar pressure tier and a
  300 K temperature tier to the reference sweep, taking it from 140 to 200
  stored states.
- `data/cea_reference_points.csv`: regenerated with `pip install cea` at the
  extended envelope; all 200 points converged in both CEA and h2thermo, none
  were dropped.
- `examples/generate_table.py`: the full-envelope example table now spans
  10 kPa to 60 bar and equivalence ratio 0 to 1.0 (previously 1-60 bar and
  0.2-1.0), matching what the library actually supports after PR-6. All
  20,000 nodes converge.
- `docs/validation.md` and `README.md` updated throughout with the figures
  measured on the extended set, replacing the PR-2/PR-4 figures rather than
  appending alongside them.
- Four new tests in `tests/test_equilibrium.py` pinning convergence and the
  expected dissociation behaviour at 10 kPa.

## What the extension changed, and what it didn't

Bulk properties (molecular weight, density, entropy, frozen specific heats,
isentropic exponent) are essentially unchanged: all within a few thousandths
of a percent of the previous figures. This is expected; these depend on the
overall stoichiometry and the major species, which the two databases already
agreed on closely.

The equilibrium (shifting-composition) specific heats are the exception:
maximum deviation from CEA roughly doubled, from 0.556% / 0.630% to 1.176% /
1.167% for cp/cv respectively. Every point driving that maximum sits at the
new 0.1 bar tier and 2200 K or above. This is not a new problem; it is the
one PR-2 already identified and quantified, observed for the first time at a
pressure low enough to make it visible at this magnitude: radical
concentrations, and therefore the shifting contribution to the specific
heats, depend exponentially on Gibbs energies the two databases do not
reproduce identically. The ratio of equilibrium to frozen cp at 2900 K and
0.1 bar is 7.13, the largest anywhere in the tabulated envelope, against
3.33 at 1 bar; the CEA deviation follows the same shape. The fixed test
tolerance (2%, set with margin in PR-2) already absorbed this without being
loosened.

## Why 300 K and not lower

200 K remains the library's stated lower temperature bound, but no CEA
points were added there. The 200-300 K decade is a small extrapolation
beyond the now-validated range, and CEA's own equilibrium solver already
prints a low-temperature range warning close to that boundary in an
adjacent code path (see PR-6's note on `adiabatic_flame_temperature`), so
points there would need closer scrutiny than a routine extension warranted.
The gap is stated in the Scope table rather than left implicit.

## Test plan

- [x] `pip install cea && python scripts/generate_cea_reference.py` -- 200
      points written, 0 CEA convergence failures
- [x] `python examples/generate_table.py` -- 20,000 nodes, 0 non-convergent
- [x] `pytest` -- 2291 passed
- [x] `flake8 src tests scripts examples` -- clean

# PR-6

Allows `equivalence_ratio = 0.0` (pure oxidizer, no fuel present) throughout
the equilibrium and table layers, closing the limitation recorded in PR-5:
`GridSpecification` rejected a zero equivalence ratio outright, so an
exported table could never carry pyCycle's `FAR = 0` row, which a full
engine model needs for unburned sections such as an inlet or compressor.

## What is included

- `equilibrium._prepare_reactants` now rejects only negative equivalence
  ratios; zero is accepted and sets the mixture to pure oxidizer.
- `GridSpecification` accepts an equivalence-ratio axis starting at 0.0,
  while temperature and pressure remain strictly positive as before.
- `h2thermo.export.pycycle`'s module docstring is corrected: the `FAR = 0`
  row is available whenever the caller's grid includes 0.0, rather than
  being unconditionally absent.
- 13 new tests: composition and specific heat of the pure-oxidizer state,
  its high-temperature dissociation, `equilibrium_specific_heats` and
  `adiabatic_flame_temperature` at zero, rejection of negative ratios,
  `GridSpecification` boundary tests, table generation over a grid that
  includes the zero row, and export of that row's `FAR = 0` entry.

## Why zero and not some small positive floor

An equivalence ratio of exactly zero has an unambiguous physical meaning,
pure oxidizer, whereas a small positive floor would be an arbitrary
approximation to the same state. Cantera's `set_equivalence_ratio` accepts
zero directly and returns the oxidizer composition exactly, so there was no
numerical reason to avoid the boundary. Negative values remain rejected;
they correspond to nothing physical.

## Verified rather than assumed

Accepting zero without raising was the easy part; what needed checking was
whether the result is physically correct.

- At 300 K the pure-oxidizer mole fractions match the `DRY_AIR` molar ratio
  (1 : 3.76) to within floating point precision, and the specific heat is
  1010 J/(kg K) against a handbook value of 1005 J/(kg K) for dry air, 0.5%
  high. The difference is expected: `DRY_AIR` is a simplified two-component
  mixture with no argon.
- At 2900 K the same state shows measurable O2 dissociation into atomic O,
  the same qualitative high-temperature behaviour the library already
  reports for fuel-air mixtures, evidence that the equilibrium solver is
  treating oxygen chemistry consistently regardless of whether fuel is
  present.
- `adiabatic_flame_temperature` at zero correctly returns the inlet
  temperature unchanged, since there is no fuel to react. Below 300 K inlet,
  Cantera's `HP` solver emits a range warning, because with no reaction the
  state never moves away from the sub-300 K starting point. This is a
  property of that one function's iterative path, not of the library's
  primary path: `ThermoTable.generate` always calls the direct `TP` solver
  at the requested grid temperature, which was checked separately down to
  200 K, the bottom of the currently supported envelope, and does not warn
  there for any equivalence ratio, zero included.

## What was deliberately left alone

`h2thermo.interpolation` required no changes. The equivalence-ratio axis is
handled the same way as temperature and pressure throughout that module, so
a grid with a zero-valued first node interpolates correctly without any
special case. This was confirmed by reasoning about the code path rather
than by a new test, since the existing interpolation test suite already
exercises arbitrary axis values.

No new CEA reference points were added for the pure-air state. The stored
140-point reference set validates hydrogen combustion products, which is
not what a fuel-free mixture is; comparing it against a handbook air
property instead of CEA is the more honest reference for what is actually
being checked.

## Test plan

- [x] `pytest` -- 1625 passed
- [x] `flake8 src tests scripts examples` -- clean

# PR-5
## Summary
- Determine, by reading pyCycle's own shipped reference table (`pycycle.constants.AIR_JETA_TAB_SPEC`) rather than assuming, that its tabular thermo format stores equilibrium (shifting-composition) `Cp`/`Cv` and an independently evaluated isentropic exponent as `gamma` -- not the frozen values or their ratio. Measured with `scripts/probe_pycycle_definitions.py` and recorded in `docs/validation.md` section 5, replacing an earlier version of that script that drove a live, version-fragile OpenMDAO model instead of the shipped artifact.
- Add `h2thermo.export.pycycle.write_pycycle_table`, which writes a `ThermoTable` in the pickle format pyCycle's tabular thermo mode reads: `Cp`/`Cv`/`gamma` mapped per the measured result above, the equivalence-ratio axis converted to fuel-air ratio via a mechanism-derived stoichiometric constant, and the specific gas constant recovered from mean molecular weight.
- Known limitation, documented in the module docstring: `GridSpecification` requires equivalence ratio > 0, so an exported table never has pyCycle's `FAR = 0` (pure air) row.

## Test plan
- [x] `pytest` -- 1612 passed
- [x] `flake8 src tests scripts examples` -- clean
- [x] `python scripts/probe_pycycle_definitions.py` (requires `pip install om-pycycle`) -- reproduces the measured discriminator values in `docs/validation.md` section 5

# PR-4
Adds specific heats that account for shifting chemical equilibrium, closing the
largest outstanding gap in the physics. What was previously documented as a
limitation is now a deliberate choice between two tabulated definitions.

## What is included

- `equilibrium_specific_heats()`: shifting `cp`, `cv` and the isentropic
  exponent at a single state
- Three new tabulated fields: `cp_equilibrium`, `cv_equilibrium`,
  `isentropic_exponent`, available through both the table and the interpolator
- The CEA reference file gains columns for frozen `cv`, equilibrium `cv` and
  the isentropic exponent
- 572 new tests, giving 1603 in total
- File format version raised to 2

## Why both definitions are kept

Frozen and equilibrium specific heats are both physically meaningful limits
rather than an approximation and a correction. Frozen values apply when the
flow moves faster than the chemistry can respond, as in a turbine. Equilibrium
values apply when the composition has time to adjust, as in a combustor. Real
behaviour lies between them, which is why CEA reports both, and this library
now follows the same convention.

The previous framing, which treated the frozen value as a shortcoming, was
misleading and has been corrected in the README and in `docs/validation.md`.

## Agreement with CEA

Measured across all 140 stored reference points:

| Quantity | Maximum relative deviation | Mean |
| --- | --- | --- |
| Frozen cp | 0.147 % | 0.060 % |
| Frozen cv | 0.196 % | 0.080 % |
| Equilibrium cp | 0.556 % | 0.207 % |
| Equilibrium cv | 0.630 % | 0.236 % |
| Isentropic exponent | 0.080 % | 0.032 % |

The equilibrium quantities agree less closely than the frozen ones, which is
expected: they depend on how fast the composition shifts, and the radical
concentrations that drive that shift are the part of the composition the two
thermodynamic databases reproduce least closely. Tolerances in the test suite
are tiered accordingly.

## Design decisions

**The isentropic exponent is computed separately, not derived from cp and cv.**
For a reacting mixture the exponent relating pressure to specific volume along
an isentrope is not the ratio of the specific heats, because the composition
shifts during the process. Deriving it that way would have produced a plausible
looking number that is wrong by several per cent in exactly the conditions
where it matters. It is obtained instead from the volume derivatives, following
the relations used by CEA, and agrees to 0.08 per cent, the closest agreement
of any quantity in this library.

**Derivatives are evaluated by central differences, and the step size was
verified rather than assumed.** Varying the temperature step from 0.1 to 10 K
changes the result only in the fourth significant figure. A test pins this,
because a derivative that moved with the step size would be unreliable
regardless of how well it happened to agree with a reference.

**The calculation is unconditional rather than opt-in.** It costs four extra
equilibrium solves per node, raising generation time by a factor of 4.1, so a
full 50 x 20 x 20 grid now takes about two minutes rather than six seconds. An
opt-in flag was considered and rejected: conditional fields would complicate
the file format, the interpolator and the tests, and earlier measurements
showed that even coarse grids are accurate enough that generation time is not
a binding constraint. Two minutes is paid once.

**The file format version guard earned its place.** Raising it to 2 means
tables written by earlier versions now fail loudly instead of being read with
three missing fields.

## Tests worth noting

Beyond agreement with CEA, several tests assert physical behaviour that would
survive a change of reference data:

- The equilibrium specific heat can never fall below the frozen one, since
  shifting equilibrium only ever absorbs additional energy.
- The two coincide below the onset of dissociation, at 800 K.
- The gap narrows as pressure rises, because dissociation increases the number
  of moles and is therefore suppressed by pressure.

# PR-3

Adds the interpolation layer, which is what makes the generated tables usable:
properties at arbitrary states within the tabulated envelope, without solving
for chemical equilibrium.

## What is included

- `ThermoInterpolator`: interpolated lookup over a `ThermoTable`
- `InterpolatedState`: the returned property set
- Scalar and array queries, with arguments broadcast against one another
- `mole_fractions()` for interpolated composition
- 19 new tests covering accuracy, the query interface, bounds and performance
- `docs/validation.md` gains a section quantifying interpolation error
- scipy added as a runtime dependency

## Accuracy

Measured on a 30 x 12 x 12 grid against direct equilibrium solves at random
states inside the envelope:

| Property | Maximum relative error |
| --- | --- |
| Ratio of specific heats | 0.016 % |
| Frozen specific heat | 0.039 % |
| Specific entropy | 0.046 % |
| Mean molecular weight | 0.047 % |
| Density | 0.043 % |

Interpolation error is an order of magnitude below the agreement with CEA, so
the total error stays dominated by the underlying thermodynamic data rather
than by tabulation.

## Design decisions

**Density is derived from the equation of state, not interpolated.**
Interpolating the density field directly gives a maximum error of 2.21 per
cent, because density varies almost linearly with pressure and is therefore
poorly represented on a logarithmic axis. Recomputing it from the interpolated
mean molecular weight gives 0.043 per cent, an improvement of roughly fifty
times for one line of arithmetic.

**Scalar queries take a dedicated trilinear path.** Routing a single lookup
through the general purpose interpolator costs about 215 us, which is slower
than solving for equilibrium outright at 99 us. Since cycle codes frequently
query one state at a time, that would have defeated the purpose of tabulation
for the case that matters most. A hand written trilinear path brings this to
34 us.

**Out of range queries raise rather than extrapolate.** The underlying data is
strongly non-linear outside the sampled region, so silent extrapolation would
return plausible looking but meaningless numbers. Callers who prefer NaN can
construct the interpolator with `bounds_error=False`.

**Pressure is interpolated on a logarithmic coordinate.** This was expected to
matter more than it does: the maximum specific heat error falls from 0.022 to
0.015 per cent on a 50 x 20 x 20 grid. It is retained because it is the
physically natural coordinate and costs nothing, but it is not the improvement
it was anticipated to be.

## Performance

Measured against a direct solve that reuses a single Cantera solution object,
which is the fastest way to call the solver. Comparing against the naive path,
where a solution object is created per call, would overstate the benefit by
roughly a factor of fifty.

| Operation | Time per state | Speed-up |
| --- | --- | --- |
| Equilibrium solve | 99 us | 1x |
| Scalar lookup | 34 us | 2.9x |
| Batched lookup | 0.87 us | 115x |

The batched path is where the advantage is decisive, so callers with many
states should pass arrays rather than looping. This is documented in the
README.

## Composition accuracy is tiered, and documented as such

Interpolated mole fractions degrade as the species becomes rarer. For the
hydroxyl radical, error rises from 0.5 per cent above a mole fraction of 0.01
to 42 per cent below 0.0001. Species influence the bulk properties in
proportion to their abundance, so this is benign for property tables, but it
is not acceptable when trace composition is itself the quantity of interest.
The docstring directs those callers to the equilibrium solver.

## Grid resolution

Maximum specific heat interpolation error is 0.097 per cent on a 20 x 8 x 8
grid and 0.015 per cent on 50 x 20 x 20. Even the coarsest grid stays well
inside the accuracy of the underlying data, so resolution can be chosen for
file size and generation time rather than for accuracy.

# PR-2

Adds a validation layer establishing the accuracy of the equilibrium
properties against NASA CEA, together with documentation of the results.

## What is included

- `scripts/generate_cea_reference.py`: produces reference states using the
  NASA CEA Python package
- `data/cea_reference_points.csv`: 140 stored reference states spanning
  600-2900 K, 1-60 bar and equivalence ratios from 0.2 to 1.0
- `tests/test_validation.py`: compares every stored state against the library
- `docs/validation.md`: full results, including the measured cost of both
  known limitations
- README updated to report the headline agreement and link to the details

## Results

| Property | Maximum deviation from CEA |
| --- | --- |
| Mean molecular weight | 0.052 % |
| Density | 0.052 % |
| Specific entropy | 0.032 % |
| Frozen specific heat | 0.147 % |
| Specific enthalpy | 11.4 kJ/kg absolute |

Internal consistency checks on element conservation, the ideal gas equation of
state and the ratio of specific heats hold to machine precision.

## Design decisions

**The CEA product species list is restricted to match `h2o2.yaml`.** This is
the most consequential choice in the change. Comparing against CEA's full
product set would conflate two independent sources of disagreement: the
different species available to each solver, and the different thermodynamic
databases behind them. Matching the species lists isolates the second, which is
what the comparison is meant to measure. The cost of the restriction is then
quantified separately in `docs/validation.md`, so nothing is hidden by it.

**Reference data is stored rather than generated during testing.** CEA is a
compiled Fortran package and would make the test suite fragile in continuous
integration. Generating the reference file is a separate, infrequent step, and
the resulting CSV is small enough to version. The committed file also serves as
durable evidence of the accuracy claim.

**Composition tolerances are tiered by species type.** Stable species agree to
0.30 % on average while radicals agree to 4.8 %. This is expected rather than
concerning: radical concentrations depend exponentially on Gibbs energies, so
small database differences are amplified. A single tolerance would be either
too loose for the stable species or produce false failures on radicals.

**All tolerances were measured before being set, then given margin.** No
tolerance in this change was chosen to make a test pass.

## Limitations, now measured rather than asserted

The two known limitations were previously documented as qualitative caveats.
They are now quantified.

**Frozen specific heats.** The ratio of equilibrium to frozen `cp` at
stoichiometric conditions:

| Temperature | 1 bar | 5 bar | 20 bar | 60 bar |
| --- | --- | --- | --- | --- |
| 2600 K | 2.07 | 1.57 | 1.34 | 1.23 |
| 2900 K | 3.33 | 2.19 | 1.69 | 1.45 |

Below 2000 K the difference is under one per cent. Above it the omission is
significant, and a test now pins this expectation so it cannot be forgotten;
that test will need updating when shifting specific heats are added.

A second test asserts that the ratio falls with pressure. Dissociation
increases the number of moles, so higher pressure suppresses it and the
shifting contribution shrinks. Reproducing this trend is evidence that the
reference data carries the physics rather than arbitrary numbers.

**Absent nitrogen chemistry.** `h2o2.yaml` treats nitrogen as inert. Measured
against CEA with nitrogen chemistry enabled, at phi = 0.6 and 20 bar, nitric
oxide reaches a mole fraction of 1.9 % at 2800 K while affecting the mean
molecular weight by 0.012 % and the frozen `cp` by 0.035 %. The thermodynamic
consequence is two orders of magnitude below the agreement already achieved, so
the omission is acceptable for property tables. It would not be acceptable for
emissions work.

## Notes

The NASA CEA reimplementation is now installable from PyPI with `pip install
cea`, which is what made programmatic validation across the full envelope
practical. Earlier plans assumed reference points would be transcribed by hand
from the CEARUN web interface.

# PR-1
Adds the first working version of the thermodynamic table layer.

`ThermoTable` samples equilibrium combustion product properties on a
structured grid of temperature, pressure and equivalence ratio, and persists
the result to a compressed NumPy archive.

## What is included

- `GridSpecification`: validated sampling grid with a `linear` constructor
- `ThermoTable.generate`: evaluates equilibrium at every node, reusing a
  single Cantera `Solution` object
- `ThermoTable.save` / `ThermoTable.load`: round trip through `.npz`
- Species mole fraction fields alongside the scalar properties
- 19 new tests covering validation, generation, physical trends and
  persistence

## Design decisions

**Equivalence ratio as the mixture coordinate.** The library speaks the
language of combustion internally. Conversion to the fuel-air ratio expected
by engine cycle codes is deferred to the export layer, keeping the two
concerns separate.

**Non-convergent nodes are recorded as NaN rather than raising.** A single
pathological point should not discard a long generation run. The count is
exposed through `failed_node_count` so failures stay visible.

**Species mole fractions are stored, not only scalar properties.** The size
cost is small for the current mechanism and the composition data is needed
for CEA validation and for any later emissions work.

**The file format carries a version number.** Future format changes will fail
loudly on old files instead of being misread silently.

## Performance

A 50 x 20 x 20 grid spanning 200-3000 K, 1-60 bar and phi 0.2-1.0 generates
20,000 nodes in roughly 9 seconds with no convergence failures, producing a
2.2 MB archive.

## Not included

Interpolation and cycle code export adapters are deliberately left for
separate changes.
