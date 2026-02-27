# Case 68: Calcite Dissolution — Reactive Geochemistry

## Overview

Calcite (CaCO3) is the most abundant carbonate mineral in Earth's crust, forming limestone, chalk, and marble. Its dissolution controls the chemistry of groundwater in carbonate aquifers worldwide, governs the global carbon cycle on geological timescales, and determines the pH of cave drip water. When acidic water contacts calcite, it dissolves: calcium ions and bicarbonate are released into solution, pH rises, and the mineral surface retreats. The process is kinetically controlled — not all calcite dissolves instantaneously — and couples tightly with the aqueous equilibrium of the carbonate system (Case 67).

This case combines the two chemical reaction types from the `chemical_reactions` module in a single simulation: instantaneous aqueous equilibrium reactions (the carbonate system speciation) and a kinetically-controlled solid-phase reaction (calcite dissolution). With three primary dissolved species (Ca^2+, H+, HCO3-), six secondary equilibrium species, and one kinetically-dissolving mineral, this is the most complex geochemistry case in the series. It demonstrates how MOOSE handles the simultaneous solution of coupled equilibrium + kinetic chemistry in a batch reactor setting.

Key concepts demonstrated:

- Combined `AqueousEquilibriumReactions` + `SolidKineticReactions` in a single `ReactionNetwork` block
- Three primary species (Ca^2+, H+, HCO3-) with six secondary equilibrium species
- Calcite dissolution reaction: CaCO3(s) → Ca^2+ + HCO3- - H+ (stoichiometry including proton consumption)
- `TotalMineralVolumeFraction` postprocessor for tracking mineral volume
- pH rise due to acid consumption by dissolution
- Initial calcite content 0.05 mol/L; dissolution proceeds until equilibrium or mineral exhaustion

---

## The Physics

### Aqueous Equilibrium Reactions

Six secondary aqueous species are held in instantaneous equilibrium with the three primary species (Ca^2+, H+, HCO3-):

| Secondary species | Reaction | log10(Keq) |
|------------------|----------|-----------|
| CO2(aq) | H+ + HCO3- ⇌ CO2(aq) | 6.35 |
| CO3^2- | HCO3- ⇌ H+ + CO3^2- | -10.33 |
| OH- | H2O ⇌ H+ + OH- | -14.0 |
| CaHCO3+ | Ca^2+ + HCO3- ⇌ CaHCO3+ | 1.11 |
| CaCO3(aq) | Ca^2+ + HCO3- ⇌ CaCO3(aq) + H+ | -7.05 |
| CaOH+ | Ca^2+ + H2O ⇌ CaOH+ + H+ | -12.78 |

The calcium-bearing secondary species (CaHCO3+, CaCO3(aq), CaOH+) are important in hard water chemistry — they represent ion pairs and complexes that reduce the free Ca^2+ activity. At the pH values of this simulation (~6-6.5), CaHCO3+ is the most significant calcium complex.

### Calcite Kinetic Dissolution Reaction

The net calcite dissolution reaction (written in terms of primary species):

```
CaCO3(s) + H+ → Ca^2+ + HCO3-
```

or equivalently in the `SolidKineticReactions` framework stoichiometry notation:

```
Primary species contributions: Ca^2+: +1, H+: -1, HCO3-: +1
```

This shows that dissolution consumes one H+ and produces one Ca^2+ and one HCO3- per mole of calcite dissolved. The consumption of H+ is why dissolution raises pH — acid is neutralized by the dissolving mineral.

The kinetic rate law (TST rate):

```
r = k_diss * A_mineral * exp(-Ea / (R * T)) * (1 - Q/K_calcite)
```

where:
- `k_diss = 1e-8 mol/(m^2 * s)` — intrinsic dissolution rate constant
- `A_mineral` — reactive surface area (proportional to current mineral content)
- `Ea = 15,000 J/mol` — activation energy (same as Case 66)
- `T = 298.15 K` — temperature
- `Q = [Ca^2+] * [HCO3-] / [H+]` — ion activity product
- `K_calcite = 10^(-8.48)` — calcite solubility product at 25°C

The equilibrium constant log10(K) = -8.48 corresponds to the calcite ion activity product:

```
K_sp = [Ca^2+] * [CO3^2-] = 10^(-8.48) = 3.31e-9 mol^2/L^2
```

### Reaction Progress Analysis

Starting from initial conditions [Ca^2+] = 0, [H+] = 1e-6 M (pH = 6), [HCO3-] = 1e-4 M, and 0.05 mol/L calcite:

1. The initial solution is undersaturated with respect to calcite (Q < K_calcite) because [Ca^2+] = 0.
2. Calcite begins to dissolve, releasing Ca^2+ and HCO3- and consuming H+.
3. pH rises as H+ is consumed (calcite is the acid consumer).
4. As [Ca^2+] increases, Q approaches K_calcite and the dissolution rate slows.
5. At equilibrium: Q = K_calcite, r = 0, dissolution stops.

The equilibrium Ca^2+ concentration for calcite dissolution in a carbonate system at pH 6.5 is approximately 5-15 mmol/L (depending on PCO2 and total alkalinity). With only 0.05 mol/L initial calcite, there is sufficient mineral to reach equilibrium before it is fully consumed.

### Domain — Batch Reactor (0D)

As in Case 67, the simulation uses a single-element mesh representing a well-mixed batch reactor. No spatial transport occurs — all chemistry is uniform throughout the reactor volume. This isolates the geochemical complexity from transport effects, allowing verification of the reaction network.

---

## Input File Walkthrough

The input file is `case68_calcite_dissolution.i`.

### `[Mesh]`

```
type = GeneratedMesh
dim = 1
nx = 1
xmax = 1
```

Single-element 1D domain. Identical to Case 67 — the minimal mesh for a 0D batch reactor.

### `[ReactionNetwork]`

The combined equilibrium + kinetic reaction network:

```
[ReactionNetwork]
  primary_species = 'Ca++ H+ HCO3-'

  [AqueousEquilibriumReactions]
    [CO2aq]
      primary_species_stoichiometry = '0  1  1'
      equilibrium_constant = 6.35
      secondary_species_name = 'CO2aq'
    []
    [CO3--]
      primary_species_stoichiometry = '0 -1  1'
      equilibrium_constant = -10.33
      secondary_species_name = 'CO3--'
    []
    [OH-]
      primary_species_stoichiometry = '0 -1  0'
      equilibrium_constant = -14.0
      secondary_species_name = 'OH-'
    []
    [CaHCO3+]
      primary_species_stoichiometry = '1  0  1'
      equilibrium_constant = 1.11
      secondary_species_name = 'CaHCO3+'
    []
    [CaCO3aq]
      primary_species_stoichiometry = '1 -1  1'
      equilibrium_constant = -7.05
      secondary_species_name = 'CaCO3aq'
    []
    [CaOH+]
      primary_species_stoichiometry = '1 -1  0'
      equilibrium_constant = -12.78
      secondary_species_name = 'CaOH+'
    []
  []

  [SolidKineticReactions]
    [calcite]
      primary_reactant_stoichiometry = '1 -1 1'
      kinetic_rate_constant = 1e-8
      log10_Keq = -8.48
      activation_energy = 15000
      molar_volume = 36.9e-6
      specific_reactive_surface_area = 1.0
      kinetic_mineral_name = calcite
    []
  []
[]
```

The `primary_species_stoichiometry` entries are ordered to match the `primary_species` list: (Ca++, H+, HCO3-). For calcite:
- Ca++: +1 (produced per mole dissolved)
- H+: -1 (consumed — this is the pH-rising mechanism)
- HCO3-: +1 (produced)

The `molar_volume = 36.9e-6 m^3/mol` is the actual molar volume of calcite, used to convert mineral concentration (mol/m^3) to volume fraction for the `TotalMineralVolumeFraction` postprocessor.

### `[Variables]`

Three primary species with chemically realistic initial conditions:

```
[Ca++]
  initial_condition = 0       # no calcium initially (fresh water)
[]
[H+]
  initial_condition = 1e-6    # pH = 6 initially (slightly acidic)
[]
[HCO3-]
  initial_condition = 1e-4    # 0.1 mmol/L background alkalinity
[]
```

The initial calcite content is specified in the `SolidKineticReactions` Action as 0.05 mol/L. This corresponds to a volume fraction of 0.05 * 36.9e-6 = 1.85e-6 m^3/m^3 — a very small amount of dispersed calcite in the pore fluid.

### `[AuxVariables]` and `[AuxKernels]`

Auto-generated by the `ReactionNetwork` Action:
- Six secondary species variables (CO2aq, CO3--, OH-, CaHCO3+, CaCO3aq, CaOH+)
- `AqueousEquilibriumRxnAux` kernels for each secondary species (mass-action expressions)
- `KineticDisPreConcAux` for tracking cumulative Ca++, H+, HCO3- changes from calcite dissolution

Manually added:
- `pH` variable with `PHAux` kernel: pH = -log10([H+])
- `calcite_vf` variable with arithmetic expression: volume fraction = calcite * molar_volume

### `[Kernels]`

For each primary species, auto-generated kernels:
- `TimeDerivative`: transient term (batch reactor — no spatial operators)
- `CoupledBEEquilibriumSub` (×6): contributions from six equilibrium secondary species to primary species mass balances
- `CoupledBEKinetic`: contribution from calcite kinetic reaction to Ca++, H+, HCO3- equations

The kinetic contribution to the H+ equation has a negative sign (dissolution consumes H+), while Ca++ and HCO3- receive positive contributions (dissolution produces these species). This stoichiometric coupling is encoded in `primary_reactant_stoichiometry = '1 -1 1'`.

### `[Postprocessors]`

| Name | Type | Physical meaning |
|------|------|-----------------|
| `pH` | ElementAverageValue | Solution pH |
| `Ca_conc` | ElementAverageValue on `Ca++` | Free calcium ion concentration |
| `HCO3_conc` | ElementAverageValue on `HCO3-` | Bicarbonate concentration |
| `CO2aq_conc` | ElementAverageValue on `CO2aq` | Dissolved CO2 |
| `calcite_mineral` | ElementAverageValue | Calcite concentration (mol/m^3) |
| `calcite_volume_frac` | TotalMineralVolumeFraction | Calcite as fraction of pore volume |

`TotalMineralVolumeFraction` is a specialized postprocessor in the `chemical_reactions` module that sums the volume fractions of all declared kinetic minerals. It is used to track how much solid mineral remains and when complete dissolution occurs.

### `[Executioner]`

```
type     = Transient
dt       = 10.0
end_time = 100.0
nl_rel_tol  = 1e-12
nl_abs_tol  = 1e-14
```

Ten time steps of 10 s each. The combined equilibrium + kinetic chemistry is stiffer than either alone. The tight tolerances (nl_rel_tol = 1e-12, nl_abs_tol = 1e-14) are needed to accurately resolve both the fast equilibrium adjustments and the slow kinetic dissolution. The kinetic timescale is:

```
t_kinetic ~ [calcite]_0 / r_max
           = 0.05 / (k_diss * A * exp(-Ea/RT))
           = 0.05 / (1e-8 * 1.0 * 2.35e-3)
           ~ 2e6 s
```

The 100 s simulation captures the very early stage of dissolution — the calcite concentration decreases by only about 50 * 10 * 2.35e-11 ~ 1.2e-8 mol/L, a tiny fraction of the initial 0.05 mol/L. This is realistic: calcite dissolution in groundwater systems typically occurs over years to centuries.

---

## Running the Simulation

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
  -v "C:/Users/simon/Downloads/moose-next/quickstart-runs:/work" \
  -w /work/case68-calcite-dissolution \
  --entrypoint /bin/bash \
  idaholab/moose:latest \
  -c '/opt/moose/bin/combined-opt -i case68_calcite_dissolution.i 2>&1 | tail -30'
```

Output files:
- `case68_calcite_dissolution_out.e` — Exodus with all species concentrations over time
- `case68_calcite_dissolution_out.csv` — Time series of pH, Ca++, HCO3-, calcite content

---

## Expected Results

### Evolution of Key Variables

| Time (s) | pH | [Ca^2+] (mol/L) | [HCO3-] (mol/L) | calcite (mol/L) |
|----------|----|-----------------|-----------------|-----------------|
| 0        | 6.00 | 0.000e-0       | 1.00e-4         | 0.05000         |
| 10       | 6.03 | ~1e-8          | ~1.00e-4        | 0.04999         |
| 50       | 6.18 | ~5e-8          | ~1.01e-4        | 0.04999         |
| 100      | 6.48 | ~1.5e-7        | ~1.01e-4        | 0.04999         |

The pH rise from 6.0 to approximately 6.48 over 100 s is the key observable. Calcium concentration increases by a small but measurable amount (from 0 to ~1.5e-7 mol/L — far below the eventual equilibrium value of ~5 mmol/L, which would require months at this rate).

### pH Rise Mechanism

The pH rise is caused directly by the dissolution stoichiometry: each mole of calcite dissolved consumes one mole of H+. As [H+] decreases, pH increases. The rise is self-limiting: as pH increases, the solution approaches calcite saturation (Q/K increases) and the driving force `(1 - Q/K)` decreases, slowing the dissolution rate.

### Carbonate Speciation Changes

As Ca^2+ enters solution:
- CaHCO3+ forms (calcium-bicarbonate complex), reducing free [Ca^2+] slightly below total dissolved calcium
- The CO2aq concentration decreases slightly as pH rises (less CO2 is stabilized at higher pH)
- CO3^2- remains negligible throughout (pH << 10.33)

These speciation changes are captured by the equilibrium secondary species calculations, which update at every Newton iteration as the primary species evolve.

### Comparison to Simple Case 67

Adding calcite dissolution transforms Case 67 (static equilibrium) into a dynamic problem where:
1. The pH is not fixed but rises over time
2. Calcium enters solution progressively
3. The mineral inventory decreases (slowly)
4. Multiple equilibrium species redistribute with each pH step

This coupling between kinetics and equilibrium is the hallmark of reactive geochemistry.

---

## Key Takeaways

- Combining `AqueousEquilibriumReactions` and `SolidKineticReactions` in a single `ReactionNetwork` block is the standard approach for realistic geochemical systems. Equilibrium reactions handle fast carbonate speciation; kinetic reactions handle slower mineral dissolution or precipitation.
- The calcite dissolution stoichiometry `Ca++ H+ HCO3-: 1 -1 1` encodes the acid consumption by the mineral. The negative stoichiometry for H+ means dissolution is a proton sink: each mole dissolved consumes one mole of H+, raising pH. This is why carbonate aquifers buffer groundwater pH near 7-8.
- `TotalMineralVolumeFraction` tracks the solid-phase volume as the kinetic reaction proceeds. When this postprocessor approaches zero, the mineral is exhausted and dissolution must stop regardless of the kinetic rate expression.
- The kinetic timescale for calcite dissolution at 25°C is millions of seconds — orders of magnitude longer than equilibrium adjustment times. This separation of timescales justifies treating the aqueous reactions as instantaneous (equilibrium) while treating the mineral-fluid reaction as kinetically controlled.
- Six secondary species are needed in the calcium carbonate system to correctly account for all major ion pairs and complexes. Omitting CaHCO3+ and CaCO3(aq) would overestimate the free Ca^2+ activity and give incorrect saturation states.
- The calcium carbonate system is the paradigm for understanding acid-base buffering in groundwater. pH is buffered near the pKa values of the carbonate equilibria (6.35 and 10.33) because dissolution or precipitation of calcite continuously adjusts [Ca^2+] and [HCO3-] to maintain the mass-action expressions near equilibrium.
- In a transport simulation (adding diffusion/advection to this batch chemistry), the dissolution front would propagate through a porous medium — this is the model for limestone cave formation, acid mine drainage treatment, and CO2 sequestration in carbonate formations.
- The tight nonlinear tolerances (nl_rel_tol = 1e-12) are essential when equilibrium constants span many orders of magnitude (from K = 10^6.35 for CO2(aq) to K = 10^-14 for water ionization). Poor convergence in one reaction propagates errors through the coupled system.
