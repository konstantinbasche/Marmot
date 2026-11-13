# Changelog entries (moved from the Marmot CHANGELOG, kept with the standalone material modules)

### Added
- `materials/CohesiveFrictionalCrack` (`COHESIVEFRICTIONALCRACK`): cohesive-frictional traction-separation law with
  crack dilatancy for concrete (Gudžulić & Meschke 2025). For implicit solvers (`computeTraction`) it is integrated by
  the implicit-explicit (Impl-Ex) scheme of Oliver et al. (2008), as used by Gudžulić & Meschke: damage, frictional
  slip and the tangential jump of the dilatancy are extrapolated from the two previous increments, and the contact
  indicator is frozen at the last converged increment, which gives a symmetric, positive semi-definite tangent in
  softening (analytical). For explicit time integration
  (`computeTractionExplicit`) the traction is that of the implicitly integrated internal variables for the current
  jump, without extrapolation, with the contact indicator of the current jump and without a tangent.
- Mass per unit area of the interface as an optional material parameter (`materialProperties[11]`), returned by
  `getMassPerArea()` for the inertia and the critical time step of explicit time integration; like the density of
  `LINEARELASTIC`, it throws if it was not given, which implicit static analyses never ask for.

### Changed
- `COHESIVEFRICTIONALCRACKIMPLEX` is merged into `COHESIVEFRICTIONALCRACK`, which is now the Impl-Ex law; the fully
  implicit update is no longer available for implicit solvers (it did not converge in structural simulations), but is
  the explicit update. Input decks of `COHESIVEFRICTIONALCRACKIMPLEX` change the material name and drop the
  parameters 11-16; the mass per area moves from `materialProperties[17]` to `[11]`.
- The mass per area of `COHESIVEFRICTIONALCRACK` moves from `materialProperties[12]` to `[11]`; decks with more than 12
  properties (the former layout) are rejected.

### Removed
- The contact-scheme parameter of the Impl-Ex law: the contact indicator is always frozen at the last converged
  increment (former scheme 1), the only one that worked in implicit simulations.
- The a-priori time-step control of the Impl-Ex law (tolerance, minimum controlled time increment, reference opening,
  increment budget; state variables `implexError`, `implexErrorCount`).
- The tangent by automatic differentiation (`materialProperties[11]`, tangent method); the analytical tangent is the
  only one.
