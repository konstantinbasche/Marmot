# Changelog

All notable changes to Marmot are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Changes that break existing input or code are marked
**Breaking**; changes that alter numerical results are marked **Results**.

## [Unreleased]

### Added
- The meshfree layer in-tree, imported from the BOKU-CMP-MEC-MAT repositories with their full history
  (`git subtree`): `core/MarmotMeshfreeCore` (reproducing kernel approximations, kernel functions, cells, particle
  domains, factories), `materialpoints/DisplacementMaterialPoint`,
  `materialpoints/GradientEnhancedFiniteStrainMaterialPoint`, `particles/DisplacementParticle`,
  `particles/GradientEnhancedFiniteStrainParticle` and `cells/GradientEnhancedFiniteStrainCell`. The separate
  module repositories are superseded.
- `cells/DisplacementCell`, the MPM cell for `DisplacementMaterialPoint` (Lagrangian and B-spline geometries).
- `BODYFORCE` for all displacement and gradient-enhanced finite-strain particles (it was advertised but empty, or
  threw for the SDI particles).
- `CWFCORRECTION`: `load[0]` selects the corrected traction components as a bit mask (1 = x, 2 = y, 4 = z;
  0 or no load = all), for faces where only the normal displacement is prescribed.
- Tests for all registered meshfree particles, cells and material points, and for the meshfree core.
- `elements/CohesiveInterfaceElement` with the cohesive interfaces `CILINE2` (2D, 4 nodes), `CITRI3` (3D, 6 nodes)
  and `CIQUAD4` (3D, 8 nodes): zero-thickness interface for discrete crack models, local frame, total jump, nodal
  integration; `CILINE2G`, `CITRI3G` and `CIQUAD4G` with full Gauss integration.
- `MarmotMaterialCohesive` and `MarmotMaterialCohesiveFactory` in `MarmotMechanicsCore`: base class and factory of
  cohesive traction-separation laws, `computeTraction( Response&, Tangents&, const Increment& )` with the recoverable
  and dissipated energy per unit area in the response; the element integrates them to its internal energy.
- `Tri3` shape functions and triangle Gauss rules, nodal (Newton-Cotes / Lobatto) quadrature rules for `Bar2`,
  `Quad4` and `Tri3` (`IntegrationTypes::NodalIntegration`), `getQuadraturePointInfo( shape, integrationType )`
  serving Gauss and nodal rules alike, and the triangular interface geometry
  `MarmotGeometryInterfaceElement<3, 6>`.
- Linear triangles `CPS3` (plane stress) and `CPE3` (plane strain) in `elements/DisplacementFiniteElement`, with one
  Gauss point, and the `Tri3` specialisations of `MarmotGeometryElement` (EnSight shape `tria3`).
- Distributed loads (pressure, surface traction) on `Tetra4` and `Tri3` elements (`C3D4`, `CPS3`, `CPE3`): face
  definitions `Tetra4::getBoundaryElementIndices` (Abaqus faces S1–S4, previously a stub that always threw) and
  `Tri3::getBoundaryElementIndices`, and `BoundaryElement` support for both parent shapes.

### Changed
- **Breaking:** a rank-deficient moment matrix of the reproducing kernel approximations now throws
  (`factorizeMomentMatrix`) instead of silently returning shape functions without partition of unity. This
  includes every 3D RKPM model with a single layer of particles through the thickness (all kernel centres
  coplanar): use at least two layers of kernels, larger supports, or a lower completeness order.
- **Breaking:** `DisplacementMaterialPoint::setInitialCondition` throws for every condition, including
  `geostaticstress`, which was silently ignored.
- **Results:** the gradients of reproducing kernel shape functions of completeness order ≥ 2
  (`computeMonomialBasisGradient`), the VCI basis in 3D, the weak-form correction `CWFCORRECTION` (boundary term in
  the intermediate configuration, `1/J_Y`), and the nonlocal tangent `dL/dN` of the gradient-enhanced
  finite-strain particles.
- **Results:** the subdomains of the SDI particles are deformed about the centroid of the parent particle, so that
  they tile the deformed particle; the pressure load of `DisplacementParticleSQCNIxSDI` evaluates its test
  function at the face centre of the geometry, like the load vector and the VCI boundary term.
- The displacement NSNI particles are named like all other particles: `DisplacementSQCNIxNSNI`,
  `DisplacementSNNIxNSNI`, `DisplacementSQCNI_RxNSNI` and `DisplacementSQCNI_RUxNSNI` (2D and 3D). The former names
  `Displacement/SQCNIxNSNI`, `Displacement/SNNIxNSNI`, `Displacement/R-SNNIxNSNI` and `Displacement/RS-SNNIxNSNI`
  remain registered as aliases.
- Point particles have no faces and no longer advertise `PRESSURE` / `CWFCORRECTION`.
- The cell, particle and material point registries are function-local statics; the public `register*` functions
  are unchanged.

### Fixed
- Result-file keyword of 2D interface elements (`ILINE2`): `MarmotGeometryInterfaceElement` named the four-node
  element `bar2`, a two-node cell, so the EnSight geometry was invalid; it is now `quad4`, the cell spanning both
  faces, as `hexa8` and `penta6` for the 3D interfaces.
- **Results:** sign of the hoop term in the axisymmetric finite-strain UL elements (`CX8UL`, `CX8RUL`, residual and
  stiffness): it was subtracted, so any state with nonzero hoop stress was wrong (e.g. homogeneous radial expansion
  off by a factor of roughly 5–60). Uniaxial loading with zero hoop stress was unaffected.
- Uninitialized intermediate second moments of the SQCNIxNSNI particles in the first increment.
- The inverse isoparametric map of distorted Lagrangian cells (it threw) and their point location test.
- Stack overflows of the 64-node B-spline cells on Windows (1 MB default stack).
