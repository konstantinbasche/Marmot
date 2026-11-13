.. _cohesivefrictionalcrack:

Cohesive-frictional crack model
===============================

Theory
------

``COHESIVEFRICTIONALCRACK`` is the traction-separation law of the discrete crack model for concrete of Gudžulić and
Meschke [GudzulicMeschke2025]_. It extends the cohesive law of Snozzi and Molinari by frictional sliding of the crack
faces and by roughness-induced dilatancy. It is used with the :ref:`cohesive interface element
<cohesiveinterfaceelement>`: cracks can only form along the interfaces between bulk elements, which remain linear
elastic.

All quantities are expressed in the local frame of the interface, with the normal jump :math:`u_x` and the
tangential jump :math:`\mathbf{u}_T = (u_y, u_z)`.

Traction
^^^^^^^^

.. math::

   \mathbf{t}' = (1 - d_c)\, K_c \begin{bmatrix} u_x - u_{dil} \\ \frac{\beta^2}{\kappa} u_y \\
   \frac{\beta^2}{\kappa} u_z \end{bmatrix}
   + d_c\, K_c \begin{bmatrix} u_x - u_{dil} \\ u_y - u_{\mu,y} \\ u_z - u_{\mu,z} \end{bmatrix} H,

where :math:`K_c` is the penalty stiffness that ties the intact interface together, :math:`d_c` the damage,
:math:`u_{dil}` the dilatant normal opening, :math:`\mathbf{u}_\mu` the frictional sliding and :math:`H` the contact
indicator. The first term is the cohesive response of the intact fraction of the interface, the second term the
frictional contact of the cracked fraction.

Damage
^^^^^^

Damage is driven by the history maximum :math:`\alpha_{max}` of the equivalent opening

.. math::

   \alpha_c = \sqrt{ \langle u_x \rangle^2 + \frac{\beta^2}{\kappa^2} \| \mathbf{u}_T \|^2 },

where :math:`\beta = f_s / f_t` is the ratio of shear to tensile strength and :math:`\kappa = G_{F,II} / G_F` the
ratio of the mode II to the mode I fracture energy. The bilinear softening law

.. math::

   q_c(\alpha_{max}) = \begin{cases}
     f_t \dfrac{\alpha_1 - \alpha_{max}}{\alpha_1 - \alpha_0}, & \alpha_0 \le \alpha_{max} \le \alpha_{int} \\[2ex]
     f_{t,int} \dfrac{\alpha_{ult} - \alpha_{max}}{\alpha_{ult} - \alpha_{int}},
       & \alpha_{int} < \alpha_{max} \le \alpha_{ult} \\[2ex]
     0, & \alpha_{max} > \alpha_{ult}
   \end{cases}

with :math:`\alpha_0 = f_t / K_c`, :math:`\alpha_1 = 2 G_f / f_t + \alpha_0`,
:math:`\alpha_{int} = \alpha_1 - (\alpha_1 - \alpha_0) f_{t,int} / f_t` and
:math:`\alpha_{ult} = 2 (G_F - G_f) / f_{t,int} + \alpha_{int}` defines the damage

.. math::

   d_c = \min\left( 1 - \frac{q_c}{K_c\, \alpha_{max}},\, 1 \right), \qquad d_c = 0 \text{ for }
   \alpha_{max} \le \alpha_0.

Unloading follows the secant to the origin. The parameters of the softening law follow [GudzulicMeschke2025]_, Eq.
19. With these definitions the area under the softening curve (the energy dissipated in pure mode I) is
:math:`G_F - \tfrac{1}{2} f_{t,int} (\alpha_1 - \alpha_{int})`, slightly less than :math:`G_F`.

Dilatancy and contact
^^^^^^^^^^^^^^^^^^^^^

Sliding over the asperities of a rough crack opens it in normal direction by

.. math::

   u_{dil} = \bar{u}_{dil} \tanh\left( \frac{\eta\, g(t_x)}{\bar{u}_{dil}} \| \mathbf{u}_T \| \right),

with the dilatancy coefficient :math:`\eta = \tan \psi_{dil}` and the limit opening :math:`\bar{u}_{dil}`. The
optional pressure function :math:`g` reduces the dilatancy linearly from 1 at :math:`t_x = 0` to 0 at the normal
traction :math:`t_{x,crit} < 0`. It is evaluated with the current normal traction: in contact,
:math:`t_x = K_c (u_x - u_{dil}(g(t_x)))` is a scalar implicit equation with a monotonic residual, solved by bisection
in the bracket :math:`[K_c (u_x - u_{dil}(1)), \min(0, K_c u_x)]` to the round-off of :math:`K_c (|u_x| + u_{dil}(1))`.
(Evaluating :math:`g` with the traction of the previous increment is unstable, since the gain of that explicit
coupling is of the order :math:`K_c \eta \| \mathbf{u}_T \| / |t_{x,crit}| \gg 1`.) The crack faces
are in contact, :math:`H = 1`, if :math:`u_x - u_{dil} \le 0`.

Friction
^^^^^^^^

The frictional traction :math:`\mathbf{t}_\mu = K_c (\mathbf{u}_T - \mathbf{u}_\mu)` is bounded by the Coulomb
condition

.. math::

   f_\mu = \| \mathbf{t}_\mu \| - \mu \langle -t_x \rangle \le 0,
   \qquad \dot{\mathbf{u}}_\mu = \dot\lambda_\mu \frac{\mathbf{t}_\mu}{\| \mathbf{t}_\mu \|},

with :math:`t_x = K_c (u_x - u_{dil}) H`. The flow rule is not dilatant, since the dilatancy is described separately.
The return mapping has a closed form. Out of contact the frictional resistance vanishes, so the sliding follows the
tangential jump, and the traction remains continuous when the crack closes again.

Energies
^^^^^^^^

The recoverable energy is that of the springs of Eq. 16 for the current damage and sliding,

.. math::

   \psi = \tfrac{1}{2} t_x (u_x - u_{dil}) + \tfrac{1}{2} (1 - d_c) K_c \frac{\beta^2}{\kappa} \| \mathbf{u}_T \|^2
   + d_c H \frac{\| \mathbf{t}_\mu \|^2}{2 K_c}.

The dissipation is the part of the work that is not recoverable,

.. math::

   D_{n+1} = D_n + \tfrac{1}{2} (\mathbf{t}_n + \mathbf{t}_{n+1}) \cdot \Delta [\![\mathbf{u}]\!]
   - (\psi_{n+1} - \psi_n),

with the work integrated by the trapezoidal rule, which is exact on the linear branches (elastic loading, secant
unloading, sticking, sliding at constant pressure). After complete separation in mode I, :math:`D` is the area under
the softening law.

The dilatancy is not derived from a potential: :math:`t_x` depends on :math:`\mathbf{u}_T` through :math:`u_{dil}`,
but :math:`\mathbf{t}_T` does not depend on :math:`u_x`. For :math:`\eta > 0` the work therefore depends on the
path, already in the elastic range. A closed elastic cycle
:math:`(u_x, 0) \to (u_x, \mathbf{u}_T) \to (0, \mathbf{u}_T) \to (0, 0)` does the net work
:math:`K_c u_{dil} u_x`, which appears in :math:`D`; the reversed cycle makes :math:`D` decrease. Shearing a
separated crack at :math:`u_x = 0` changes :math:`D` by
:math:`\tfrac{1}{2} K_c \eta (\mu - \eta - \mu^2 \eta) \| \mathbf{u}_T \|^2` while the dilatancy is linear,
which is negative for :math:`\eta > \mu / (1 + \mu^2)`.

Time integration
----------------

After each evaluation, the internal variables are integrated implicitly for the current jump, from the history
:math:`\alpha_{max}^{old}`, :math:`\mathbf{u}_\mu^{old}` at the beginning of the increment: the history maximum of
:math:`\alpha_c`, the frictional return mapping, and the dilatancy with the implicit pressure function. These values are
the state variables. The material provides two updates of the traction, which store the same state variables.

Implicit solvers: Impl-Ex
^^^^^^^^^^^^^^^^^^^^^^^^^

For implicit solvers (``computeTraction``) the law is integrated by the implicit-explicit scheme of Oliver et al.
[OliverEtAl2008]_, as used by Gudžulić and Meschke (Online Resource 1, Appendix C). A fully implicit integration does
not converge when many interfaces soften at once. The traction of an increment is evaluated with internal variables
extrapolated linearly from the two previous converged increments, with :math:`r = \Delta t_{n+1} / \Delta t_n`:

.. math::

   \tilde\alpha_{max} = \alpha_{max,n} + r\,(\alpha_{max,n} - \alpha_{max,n-1}), \qquad
   \Delta\tilde\lambda_\mu = r\,\Delta\lambda_{\mu,n}, \qquad
   \tilde u_T = u_{T,n} + r\,(u_{T,n} - u_{T,n-1}),

i.e. the damage :math:`\tilde d_c = d_c(\tilde\alpha_{max})`, the frictional sliding, applied in the direction of the
trial friction traction (the frictional traction vanishes if the extrapolated sliding exceeds the trial sliding), and
the tangential jump that drives the dilatancy.

The contact indicator :math:`H` is frozen at the contact state of the last converged increment. With :math:`H` from the
current jump, completely damaged interfaces switch between the normal stiffness 0 and :math:`K_c` during the
iterations, which makes Newton's method cycle. The frozen indicator lags one increment instead: closing faces
penetrate without resistance, and opening faces transmit the tension :math:`K_c (u_x - u_{dil})`, for one increment.

Within an increment, damage, friction and contact cannot switch state. The normal traction depends on :math:`u_x` only,
the tangential traction on :math:`\mathbf{u}_T` only, and the tangent

.. math::

   \frac{\partial t_x}{\partial u_x} = (1 - \tilde d_c + \tilde d_c H) K_c \left( 1 - \frac{\mathrm{d} u_{dil}}{\mathrm{d}
   u_x} \right), \qquad
   \frac{\partial \mathbf{t}_T}{\partial \mathbf{u}_T} = (1 - \tilde d_c) K_c \frac{\beta^2}{\kappa} \mathbf{I}
   + \tilde d_c H \frac{\partial \tilde{\mathbf{t}}_\mu}{\partial \mathbf{u}_T}

is symmetric and positive semi-definite, also in softening. In contact with :math:`t_{x,crit} < t_x < 0`, the normal
traction is the implicit function defined by :math:`r = t_x - K_c (u_x - u_{dil}(g(t_x), \tilde u_T)) = 0`, with
:math:`D = \partial r / \partial t_x = 1 + K_c \eta \tilde u_T (1 - \tanh^2 a) \, g'(t_x)` and
:math:`a = \eta g \tilde u_T / \bar u_{dil}`, so that :math:`\mathrm{d} u_{dil} / \mathrm{d} u_x = (D - 1) / D`; otherwise
it vanishes. With :math:`\mathbf{t}_\mu^{tr} = K_c (\mathbf{u}_T - \mathbf{u}_\mu^{old})`,
:math:`\mathbf{m} = \mathbf{t}_\mu^{tr} / \|\mathbf{t}_\mu^{tr}\|` and :math:`\Delta\tilde\lambda_\mu > 0`,

.. math::

   \frac{\partial \tilde{\mathbf{t}}_\mu}{\partial \mathbf{u}_T} = K_c \left[ \mathbf{I} - \frac{K_c
   \Delta\tilde\lambda_\mu}{\|\mathbf{t}_\mu^{tr}\|} (\mathbf{I} - \mathbf{m} \otimes \mathbf{m}) \right]

if :math:`\|\mathbf{t}_\mu^{tr}\| > K_c \Delta\tilde\lambda_\mu`, and 0 else.

The extrapolation error is of the order of the time increment, so small increments are required, and the time must
be a measure of the load (e.g. proportional to the prescribed displacement). Notable consequences:

- the first increment is not extrapolated (:math:`r = 0`), nor are increments with :math:`\Delta t = 0`; the first
  increment should therefore stay in the elastic range;
- an increment with :math:`\Delta t < 10^{-6} \Delta t_n`, such as the remainder that closes a step after the rounding
  of the summed time increments, does not shift the history of the extrapolation, since the next extrapolation would
  amplify its round-off by :math:`\Delta t_{n+1} / \Delta t`;
- the first increment after a reversal still extrapolates the previous one: e.g. the first unloading increment uses
  the damage of a further loading increment;
- the friction limit :math:`\mu \langle -t_x \rangle` is not enforced within an increment; frictional traction lags
  a change of the normal traction by one increment.

The energies are evaluated with the extrapolated damage and sliding of the traction; where the extrapolation is not
exact, the bookkeeping inherits its error, and :math:`D` can decrease when the damage of the next increment falls back
to the implicit one.

Explicit time integration
^^^^^^^^^^^^^^^^^^^^^^^^^

For explicit time integration (``computeTractionExplicit``) the traction is that of the implicitly integrated internal
variables for the current jump: nothing is extrapolated, the contact indicator follows the current jump, and no
tangent is computed. The time increment does not enter the traction.

Both updates evaluate the traction with the same function, for different internal variables :math:`\alpha_{max}`,
:math:`\Delta\lambda_\mu`, :math:`u_T` and :math:`H`: extrapolated for Impl-Ex, integrated implicitly here. For the
implicit slip multiplier :math:`\Delta\lambda_\mu = (\|\mathbf{t}_\mu^{tr}\| - \mu \langle -t_x \rangle) / K_c`,
the frictional traction :math:`\mathbf{t}_\mu^{tr} - K_c \Delta\lambda_\mu \mathbf{m}` is the one of the radial return.

Material parameters
-------------------

.. list-table::
   :header-rows: 1
   :widths: 10 20 70

   * - Index
     - Symbol
     - Description
   * - 0
     - :math:`K_c`
     - penalty stiffness [stress / length]
   * - 1
     - :math:`f_t`
     - tensile strength
   * - 2
     - :math:`G_F`
     - total mode I fracture energy
   * - 3
     - :math:`G_f`
     - initial mode I fracture energy, :math:`G_f < G_F`
   * - 4
     - :math:`f_{t,int}`
     - strength at the kink of the softening law, :math:`0 < f_{t,int} \le f_t`
   * - 5
     - :math:`\beta`
     - ratio of shear to tensile strength
   * - 6
     - :math:`\kappa`
     - ratio of mode II to mode I fracture energy
   * - 7
     - :math:`\mu`
     - Coulomb friction coefficient
   * - 8
     - :math:`\eta`
     - dilatancy coefficient :math:`\tan \psi_{dil}`; 0 switches dilatancy off
   * - 9
     - :math:`\bar{u}_{dil}`
     - limit of the dilatant normal opening
   * - 10
     - :math:`t_{x,crit}`
     - optional: normal traction (< 0) at which dilatancy ceases; omitted or :math:`\ge 0`: :math:`g = 1`
   * - 11
     - :math:`\mu_A`
     - optional: mass per unit area of the interface [mass / area], required only for the inertia, e.g. in explicit
       time integration; requires property 10, which may be set to 0. The interface element takes its mass from the material
       (``getMassPerArea``, which throws if it was not given), as a bulk element takes its density.

For concrete C35/45, [GudzulicMeschke2025]_ use (N, mm): :math:`K_c = 3.4 \cdot 10^7`, :math:`f_t = 2.2`,
:math:`G_F = 0.14`, :math:`G_f = 0.056`, :math:`f_{t,int} = 0.55`, :math:`\beta = 2`, :math:`\kappa = 10`,
:math:`\mu = 0.7`, :math:`\eta = 0.2` and :math:`\bar{u}_{dil} = 16`.

State variables
---------------

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Name
     - Size
     - Description
   * - ``alphaMax``
     - 1
     - history maximum :math:`\alpha_{max}` of the equivalent opening (implicit)
   * - ``alphaMaxOld``
     - 1
     - :math:`\alpha_{max}` of the increment before, for the extrapolation
   * - ``slip``
     - 2
     - frictional sliding :math:`\mathbf{u}_\mu` (implicit)
   * - ``slipMultiplier``
     - 1
     - increment of the frictional slip multiplier :math:`\Delta\lambda_{\mu,n}` of the last increment
   * - ``tangentialJump``
     - 1
     - norm of the tangential jump :math:`u_T`
   * - ``tangentialJumpOld``
     - 1
     - :math:`u_T` of the increment before
   * - ``dTOld``
     - 1
     - time increment of the last increment that shifted the history of the extrapolation
   * - ``contact``
     - 1
     - contact indicator :math:`H` of the implicit update, frozen in the next Impl-Ex increment
   * - ``damage``
     - 1
     - damage :math:`d_c` (implicit)
   * - ``alphaC``
     - 1
     - current equivalent opening :math:`\alpha_c` (Eq. 18), for output, e.g. as crack indicator

References
----------

.. [OliverEtAl2008] J. Oliver, A. E. Huespe, J. C. Cante. An implicit/explicit integration scheme to increase
   computability of non-linear material and contact/friction problems. Computer Methods in Applied Mechanics and
   Engineering 197:1865-1889, 2008. https://doi.org/10.1016/j.cma.2007.11.027

.. [GudzulicMeschke2025] V. Gudžulić, G. Meschke. Discrete modeling of cracking in reinforced concrete structures:
   formulation, size effect, and parameter sensitivity. Materials and Structures 58:295, 2025.
   https://doi.org/10.1617/s11527-025-02815-6

Implementation
--------------

.. doxygenclass:: Marmot::Materials::CohesiveFrictionalCrack
   :allow-dot-graphs:

.. doxygenclass:: Marmot::Materials::MaterialCohesive
   :allow-dot-graphs:
