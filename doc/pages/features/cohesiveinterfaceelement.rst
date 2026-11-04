.. _cohesiveinterfaceelement:

Cohesive Interface Element
==========================

Theory
------

The cohesive interface element is a **zero-thickness interface** for discrete crack models: it is inserted between
the faces of two bulk elements, and all inelasticity of the model is described by the
cohesive material of the interface, while the bulk material can remain linear
elastic. The element has :math:`2 n` nodes: the first :math:`n` nodes form the bottom face and the last :math:`n`
nodes the top face, in the same order. Both faces coincide in the reference configuration.

The displacement jump across the interface,

.. math::

   [\![\mathbf{u}]\!] = \mathbf{u}^{+} - \mathbf{u}^{-} = \mathbf{N}_{\mathrm{jump}}\, \mathbf{q}, \qquad
   \mathbf{N}_{\mathrm{jump}} = \big[\, -\mathbf{N} \;\; +\mathbf{N} \,\big],

is rotated into the local frame of the interface, :math:`[\![\mathbf{u}]\!]' = \mathbf{R}\, [\![\mathbf{u}]\!]`.
The rows of :math:`\mathbf{R}` are the local axes:

- :math:`x'`: the unit normal of the bottom face; for a correctly oriented element it points towards the top face,
  so that a positive normal jump opens the interface. In 3D it follows the right-hand rule of the node order of the
  bottom face; in 2D it is the tangent of the bottom face, from its first to its second node, rotated by
  :math:`+90^\circ`, i.e. the top face lies on the left of the bottom face,
- :math:`y'`: the direction of the first edge of the face, :math:`\partial \mathbf{x} / \partial \xi_1`,
- :math:`z' = x' \times y'` (3D only).

The element works with :math:`n_{dim}` local components. The cohesive material always has three: an element
embedded in 2D passes a zero out-of-plane jump :math:`[\![u]\!]_{z'} = 0`, as for plane strain, and keeps the
in-plane part of the traction and the tangent. A line element with thickness :math:`b` therefore responds exactly
like a quadrilateral element that extrudes it by :math:`b`.

The cohesive material returns the local traction :math:`\mathbf{t}'` and its tangent, and the element assembles

.. math::

   \mathbf{P}_e = \sum_{qp} \mathbf{N}_{\mathrm{jump}}^\mathsf{T} \mathbf{R}^\mathsf{T}\, \mathbf{t}'\, J_0\, w,
   \qquad
   \mathbf{K}_e = \sum_{qp} \mathbf{N}_{\mathrm{jump}}^\mathsf{T} \mathbf{R}^\mathsf{T}\,
   \frac{\partial \mathbf{t}'}{\partial [\![\mathbf{u}]\!]'}\, \mathbf{R}\, \mathbf{N}_{\mathrm{jump}}\, J_0\, w,

with the surface Jacobian determinant :math:`J_0 = \sqrt{\det \mathbf{G}}` and the quadrature weight :math:`w`. The
kinematics are geometrically linear and the local frame is fixed in the reference configuration.

The material also returns the recoverable energy :math:`\psi` and the dissipated energy :math:`D` per unit area. The
internal energy of the element (``computeInternalEnergy``), :math:`\sum_{qp} ( \psi + D )\, J_0\, w`, is the work
done on the interface, e.g. for an energy balance of the model.

The integration rule of the surface is chosen at registration, through ``getQuadraturePointInfo( shape,
integrationType )``. The ``CI*`` elements use a **nodal** (Newton-Cotes / Lobatto) rule
(``NodalIntegration``): the integration points coincide with the nodes. With the high penalty stiffness of cohesive
laws, Gauss integration causes spurious oscillations of the traction along the interface, which the nodal rule
avoids. The ``CI*G`` elements use full Gauss integration (``FullIntegration``), e.g. for comparison. Reduced Gauss
integration is not registered, since a single integration point per face leaves zero-energy modes of the jump.

Explicit time integration
-------------------------

``computeKernelsExplicit`` performs the update of ``computeKernels`` without the tangent, through
``MaterialCohesive::computeTractionExplicit``. Since an explicit solver does not iterate, no Impl-Ex integration of
the cohesive law is needed; the implicit update of the law at the given jump is cheap and robust.

A zero-thickness interface has no volume and, by default, no mass: its nodes are carried by the adjacent bulk
elements. Its stable time increment, however, is set by the penalty springs between the two faces, which only the
element knows. As a bulk element takes its density from its material (``MaterialHypoElastic::getDensity``), the
interface takes its mass per unit area :math:`\mu` from the cohesive material (``MaterialCohesive::getMassPerArea``;
for ``COHESIVEFRICTIONALCRACK`` an optional material parameter). As for the density, only the inertia asks for it,
so implicit static analyses run without it, while a material without one throws when the inertia is needed. With
:math:`\mu > 0` the interface
carries an artificial mass of its own, :math:`\mu / 2` per unit area on each face (lumped:
:math:`\tfrac{\mu}{2} \int_\Gamma N_i\, dA` per node), and
``computeCriticalTimeStepForExplicitDynamics`` bounds the largest eigenvalue of
:math:`\mathbf{M}^{-1} \mathbf{K}` of the element, with
:math:`\mathbf{K} = \int_\Gamma \mathbf{N}_{\mathrm{jump}}^\mathsf{T} K_{max} \mathbf{N}_{\mathrm{jump}}\, dA` and
the maximum stiffness :math:`K_{max}` of the material (``MaterialCohesive::getMaximumStiffness``), by its Gershgorin
row sums. By default, :math:`K_{max}` is the largest diagonal entry of the tangent at the current state of each
quadrature point. A solver may evaluate the critical time step only once, at the start; a material whose stiffness
can grow later, e.g. through contact or dilatancy, overrides it with a bound over all states, as
``COHESIVEFRICTIONALCRACK`` does. For the nodal rule each node pair is an isolated spring, the bound is exact, and

.. math::

   \Delta t_{crit} = \sqrt{ \mu / K_{max} },

independent of the element size. Like a bulk element, the interface uses its own mass and stiffness only. By the
element eigenvalue theorem the minimum over all elements is then a lower bound of the critical time step of the
assembled model, also where several interfaces and bulk elements share a node; a time step based on the mass of the
adjacent bulk elements instead would count that mass twice and could overestimate the stable increment.

The mass is artificial. A target increment :math:`\Delta t` requires :math:`\mu = K_{max} \Delta t^2`, which should
be small compared to the mass per area :math:`\rho h` of the adjacent bulk elements, or acceptable as mass scaling
in a quasi-static analysis. Both are governed by the penalty number :math:`\alpha = K_c h / E` of the cohesive law:
even with the true nodal masses, the stable increment of the interface is roughly :math:`\Delta t_{bulk} /
\sqrt{\alpha}`. A penalty stiffness chosen for implicit analyses, e.g. :math:`K_c = 1000 E` per mm, is far stiffer
than needed for explicit ones; :math:`\alpha \approx 50` adds about 2 % to the compliance of the bulk elements.

Elements
--------

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Name
     - Nodes
     - Description
   * - ``CILINE2``
     - 4
     - Two-dimensional cohesive interface between two linear lines (2 bottom and 2 top nodes), for meshes of linear
       plane elements, nodal integration; written to result files as ``quad4``
   * - ``CITRI3``
     - 6
     - Three-dimensional cohesive interface between two linear triangles (3 bottom and 3 top nodes), for meshes of
       linear tetrahedra, nodal integration; written to result files as the EnSight wedge ``penta6``
   * - ``CIQUAD4``
     - 8
     - Three-dimensional cohesive interface between two bilinear quadrilaterals (4 bottom and 4 top nodes), for meshes
       of linear hexahedra, nodal integration; written to result files as ``hexa8``
   * - ``CILINE2G``, ``CITRI3G``, ``CIQUAD4G``
     - 4, 6, 8
     - As ``CILINE2``, ``CITRI3`` and ``CIQUAD4``, but with full Gauss integration (2, 3 and 2×2 points)

The class is a template in the spatial dimension and the number of nodes; all sizes follow from these two
parameters, so that a further interpolation of ``GeometryInterfaceElement`` only needs to be registered.

Properties
----------

The material of the section is created by the cohesive material factory. For elements embedded in 2D, the first
(optional) element property is the out-of-plane thickness; it defaults to 1.

The element has no named properties. The mass per unit area :math:`\mu`, required for explicit time integration
(see above), is a property of the cohesive material, as the density is one of a bulk material.

State variables
---------------

Per quadrature point the element stores, in this order: the local traction :math:`\mathbf{t}'` (``traction``,
:math:`n_{dim}` components), the local jump :math:`[\![\mathbf{u}]\!]'` (``jump``, :math:`n_{dim}` components),
the recoverable and the dissipated energy per unit area (``elastic energy density`` and ``dissipation density``, one
component each), followed by the state variables of the cohesive material.

Implementation
--------------

.. doxygenclass:: Marmot::Discretizations::Elements::CohesiveInterfaceElement
   :allow-dot-graphs:
