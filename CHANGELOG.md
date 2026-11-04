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
- `Marmot::Microplane` in `MarmotMechanicsCore`: the microplane quadrature schemes (N21 symmetric/asymmetric, N28,
  N37, N61) and the volumetric-deviatoric (`VDSplit`) and normal-tangential (`NTSplit`) projection tensors with
  their homogenization helpers, shared by all microplane models.
- `materials/MicroplaneM7` (`MICROPLANEM7`), the microplane model M7 for concrete (Caner and Bažant, 2013), and
  `materials/GM7` (`GM7`), its gradient-enhanced nonlocal extension. The former standalone module repositories are
  superseded.
- `elements/CohesiveInterfaceElement` with the cohesive interfaces `CILINE2` (2D, 4 nodes), `CITRI3` (3D, 6 nodes)
  and `CIQUAD4` (3D, 8 nodes): zero-thickness interface for discrete crack models, local frame, total jump, nodal
  integration; `CILINE2G`, `CITRI3G` and `CIQUAD4G` with full Gauss integration.
- `Marmot::Materials::MaterialCohesive` and `Marmot::Registry::MaterialCohesiveFactory` in `MarmotMechanicsCore`:
  base class and factory of cohesive traction-separation laws, `computeTraction( Response&, Tangents&,
  const Increment& )` with the recoverable and dissipated energy per unit area in the response; the element
  integrates them to its internal energy.
- Explicit time integration for `elements/CohesiveInterfaceElement`: `computeKernelsExplicit`, lumped and
  consistent inertia of an artificial interface mass (half on each face), and
  `computeCriticalTimeStepForExplicitDynamics` from that mass and `MaterialCohesive::getMaximumStiffness()`. Two new
  virtuals of the cohesive material base:
  - `MaterialCohesive::getMaximumStiffness()`, the counterpart of `MaterialHypoElastic::getMaximumWaveSpeed()`: by
    default the largest diagonal entry of the tangent at the current state, evaluated on a copy of it;
  - `MaterialCohesive::getMassPerArea( stateVars )` (pure virtual), the counterpart of
    `MaterialHypoElastic::getDensity()`: the mass per unit area of the interface. As for the bulk elements, the mass
    comes from the material, not from an element property, and only the inertia asks for it, so implicit static
    analyses need none.
- `Tri3` shape functions and triangle Gauss rules, nodal (Newton-Cotes / Lobatto) quadrature rules for `Bar2`,
  `Quad4` and `Tri3` (`IntegrationTypes::NodalIntegration`), `getQuadraturePointInfo( shape, integrationType )`
  serving Gauss and nodal rules alike, and the triangular interface geometry
  `GeometryInterfaceElement<3, 6>`.
- Linear triangles `CPS3` (plane stress) and `CPE3` (plane strain) in `elements/DisplacementFiniteElement`, with one
  Gauss point, and the `Tri3` specialisations of `GeometryElement` (EnSight shape `tria3`).
- Distributed loads (pressure, surface traction) on `Tetra4` and `Tri3` elements (`C3D4`, `CPS3`, `CPE3`): face
  definitions `Tetra4::getBoundaryElementIndices` (Abaqus faces S1–S4, previously a stub that always threw) and
  `Tri3::getBoundaryElementIndices`, and `BoundaryElement` support for both parent shapes.

### Changed
- **Breaking:** complete namespace rework. All code now lives in the global `Marmot` namespace with a closed set
  of second-level namespaces (`Utilities`, `Math`, `ContinuumMechanics`, `Materials`, `Discretizations`,
  `Registry`, `Testing`); classes dropped the redundant `Marmot` name prefix; the top-level `MarmotLibrary`
  namespace is gone. There are **no compatibility aliases** — consumers migrate with this table.
  Header file names are unchanged, with three exceptions: `Marmot/MarmotUtility.h` (empty, deleted — include
  `Marmot/MarmotTypedefs.h` instead) and `Marmot/VonMisesConstants.h` / `Marmot/ADVonMisesConstants.h`
  (dissolved into the model headers, see table):

  | Old | New |
  |---|---|
  | `MarmotLibrary::Marmot<X>Factory` | `Marmot::Registry::<X>Factory` (e.g. `ElementFactory`, `MaterialHypoElasticFactory`, `MaterialFiniteStrainFactory`, `MaterialGradientEnhancedFiniteStrainFactory`, `MaterialGradientEnhancedHypoElasticFactory`, `ParticleFactory`, `MaterialPointFactory`, `CellFactory`, `CellElementFactory`) |
  | `::MarmotMaterialHypoElastic(/AD)`, `::MarmotInterfaceMaterialHypoElastic` | `Marmot::Materials::MaterialHypoElastic(/AD)`, `::InterfaceMaterialHypoElastic` |
  | `::MarmotMaterialFiniteStrain(/AD)`, `::MarmotMaterialGradientEnhancedFiniteStrain` | `Marmot::Materials::MaterialFiniteStrain(/AD)`, `::MaterialGradientEnhancedFiniteStrain` |
  | `::MarmotMaterialGeneralGradientEnhancedHypoElastic<>` | `Marmot::Materials::MaterialGradientEnhancedHypoElastic<>` |
  | `::MarmotElement`, `::MarmotGeometryElement<>`, `::MarmotElementSpatialWrapper` | `Marmot::Discretizations::Elements::{Element, GeometryElement<>, ElementSpatialWrapper}` |
  | `::MarmotMaterialSection`, `::ElementProperties` | `Marmot::Discretizations::{MaterialSection, ElementProperties}` |
  | `::MarmotMaterialPoint`, `::MarmotCell`, `::MarmotCellElement` | `Marmot::Discretizations::{MaterialPoints::MaterialPoint, Cells::Cell, Cells::CellElement}` |
  | `::MarmotJournal`, `::MakeString`, `::StateView`, `::MarmotStateVarVectorManager` | `Marmot::Utilities::{Journal, MakeString, StateView, StateVarVectorManager}` |
  | `Marmot::Elements` / `Marmot::Meshfree` / `Marmot::MaterialPoints` / `Marmot::Cells` | `Marmot::Discretizations::{Elements, Meshfree, MaterialPoints, Cells}` (meshfree classes dropped prefixes: `Particle`, `Approximation`, `KernelFunction*`) |
  | `Marmot::FiniteElement::*` (shape functions, quadrature, EAS, …) | `Marmot::Discretizations::Elements::*` (`Spatial2D::axisymmetric` → `Axisymmetric`) |
  | `Marmot::Constants`, `Marmot::EigenTensors` + `ContinuumMechanics::{CommonTensors, TensorUtility}` | `Marmot::Math::{Constants, Tensors}` |
  | `Marmot::FastorStandardTensors`, `Marmot::FastorIndices`, `ContinuumMechanics::TensorUtility::TensorExponential`, plain-`Marmot` Fastor helpers (`mapEigenToFastor`, `makeDual`/`makeReal`, `deviatoric`, `expandTo3D`/`reduceTo2D`, …) | `Marmot::Math::Tensors::{FastorTensors, FastorIndices, TensorExponential}` and functions directly in `Marmot::Math::Tensors` |
  | `Marmot::AutomaticDifferentiation`, `Marmot::NumericalAlgorithms::Differentiation` | `Marmot::Math::Differentiation::{Automatic, Numerical}` |
  | `Marmot::NumericalAlgorithms::Integration` | `Marmot::Math::Integration` |
  | `Marmot::NumericalAlgorithms` (substeppers, `HughesWinget`, `YieldSurfaceCombinationManager`) | `Marmot::ContinuumMechanics::Algorithms` |
  | `Marmot::TimeIntegration`, `Marmot::GeostaticStress` | under `Marmot::ContinuumMechanics` |
  | `Marmot::GradientDamage`, `Marmot::PhaseField` | `Marmot::ContinuumMechanics::Degradation::{GradientDamage, PhaseField}` |
  | `Marmot::Solvers::MarmotMaterialPointSolver{HypoElastic, FiniteStrain}` | `Marmot::Materials::MaterialPointSolver{HypoElastic, FiniteStrain}` |
  | `ContinuumMechanics::FiniteStrain::*` (namespace dissolved) | `ContinuumMechanics::Kinematics::DeformationMeasures`, `::Elasticity::EnergyDensityFunctions`, `::StressMeasures`, `::Plasticity`, `::Viscoelasticity` (merged with the small-strain one) |
  | `ContinuumMechanics::CommonConstitutiveModels` (namespace dissolved) | `ContinuumMechanics::Plasticity::MenetreyWillam`, `ContinuumMechanics::Viscoplasticity::DuvautLionsViscosity` |
  | `Marmot::Materials::{KelvinChain, Wiechert, SolidificationTheory}` | `Marmot::ContinuumMechanics::Viscoelasticity::{KelvinChain, Wiechert, SolidificationTheory}` |
  | `ContinuumMechanics::VoigtNotation` (conversions, `Transformations`) | `ContinuumMechanics::Voigt` |
  | `ContinuumMechanics::VoigtNotation::{Invariants, Derivatives}`, `ContinuumMechanics::HaighWestergaard` | `ContinuumMechanics::Invariants`, `::Invariants::Derivatives`, `::Invariants::HaighWestergaard` |
  | derivative namespaces `*OrderDerived` spellings | unified `FirstOrder` / `SecondOrder` / `ThirdOrder` |
  | `Marmot::Materials::Shrinkage::B4` | `Marmot::Materials::B4Shrinkage` |
  | `Marmot::Materials::VonMisesConstants` / `::ADVonMisesConstants` (headers `VonMisesConstants.h`, `ADVonMisesConstants.h`) | class members `VonMisesModel::Constants` / `ADVonMises::Constants` (in `VonMises.h`, `ADVonMises.h`) |
  | `Marmot::Microplane` (quadrature schemes, `VDSplit`, `NTSplit`) | `Marmot::ContinuumMechanics::Microplane` |
  | plain-`Marmot` typedefs (`Vector6d`, `Matrix6d`, `mVector6d`, …), general math functions and B-spline/monomial bases | `Marmot::Math::*` (root) |

  Factory registration **string keys are unchanged**, so name-based material/element creation through the
  factories behaves identically. Registration translation units now use unnamed namespaces. The convention and
  its rules are documented in `AGENTS.md`; `scripts/checkNamespaceConventions.sh` checks them.
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
- The reference wave speed of the artificial bulk viscosity is cached also when it is zero (it was asked for again on
  every increment), in `DisplacementFiniteElement` and `GeneralGradientEnhancedDisplacementFiniteElement`. A material
  meant to run with the damage degradation should override `getMaximumWaveSpeed()` with its closed form; the new
  `Marmot::Testing::checkMaximumWaveSpeedOverrideMatchesDefault()` (`MarmotWaveSpeedContract.h`) checks an override
  against the default.
- `GeneralGradientEnhancedDisplacementFiniteElement` no longer stores the strain-displacement operator per quadrature
  point: `QuadraturePoint::B` and the unused `B_K` are removed (code reading them must be rebuilt and use
  `localGeometryElement.B( qp.dNdX )`), and `B` is built from `dNdX` where the kernels use it. 16.7 of 20.7 GB on a
  107,508-element `GC3D20` mesh; results are bitwise identical.
- **Results:** the kernels of `GeneralGradientEnhancedDisplacementFiniteElement` (`GC3D8`, `GC3D20`, `GCPS4`, ... and
  the mixed-order `M` variants) compute the strain increment and the internal force from `dNdX` through
  `Discretizations::Elements::SmallStrain::strainIncrement` and `addInternalForce` (new, in
  `MarmotSmallStrainKinematics.h`)
  instead of the dense B-matrix, and the gradient term of the non-local residual and tangent in O(nNodes) per
  quadrature point. Residuals and tangents change in the last digits (about 1e-16 relative to the residual scale);
  the explicit edge-breakout run is 25 % faster. The quadrature points are smaller (about 5 KB to under 1 KB per
  point for `GC3D20R`); equal-order elements store the non-local shape functions once.

### Fixed
- Lumped inertia and critical time step of the linear simplices `CPS3`, `CPE3` and `C3D4`: the manifold-based mass
  lumping blended in the shape functions of a quadrilateral or hexahedron and wrote past the end of the nodal
  vectors, which corrupted the heap of an explicit analysis. A linear simplex is now its own linear counterpart
  (`MassLumping::nLinearNodes()`, also used by the finite-strain and the gradient-enhanced elements), with the
  mass / nNodes at each node. **Results**
- **Results:** elastic strain energy of `LinearElastic`, `ADLinearElastic` and every material on
  `MaterialHypoElasticAD`. The fixes:
  - `LinearElastic` did not compute the elastic energy density at all, so the internal energy of the elements was
    zero.
  - `ADLinearElastic` accumulated ½ Δε : Δσ without the stress at the start of the increment, so it was too small
    by about the number of increments.
  - `MaterialHypoElasticAD::computeStress` returned only the stress of `computeStressAD` and dropped its energy and
    dissipation.
  - Both linear elastic laws now integrate W by the trapezoidal rule, which is exact for linear elasticity.
  - Tests compare W with ½ ε : C : ε after a non-proportional strain path in 3D and in plane stress.
- Result-file keyword of 2D interface elements (`ILINE2`): `GeometryInterfaceElement` named the four-node
  element `bar2`, a two-node cell, so the EnSight geometry was invalid; it is now `quad4`, the cell spanning both
  faces, as `hexa8` and `penta6` for the 3D interfaces.
- **Results:** sign of the hoop term in the axisymmetric finite-strain UL elements (`CX8UL`, `CX8RUL`, residual and
  stiffness): it was subtracted, so any state with nonzero hoop stress was wrong (e.g. homogeneous radial expansion
  off by a factor of roughly 5–60). Uniaxial loading with zero hoop stress was unaffected.
- Uninitialized intermediate second moments of the SQCNIxNSNI particles in the first increment.
- The inverse isoparametric map of distorted Lagrangian cells (it threw) and their point location test.
- Stack overflows of the 64-node B-spline cells on Windows (1 MB default stack).
