/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Konstantin Viktor Basche konstantin.basche@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */

#pragma once

#include "Marmot/MarmotMaterialCohesive.h"
#include <Eigen/Core>

namespace Marmot::Materials {

  /**
   * @brief Cohesive-frictional traction-separation law with crack dilatancy for concrete, integrated by the
   * implicit-explicit (Impl-Ex) scheme for implicit solvers and directly for explicit time integration.
   *
   * Discrete crack model of Gudžulić & Meschke, *Discrete modeling of cracking in reinforced concrete structures:
   * formulation, size effect, and parameter sensitivity*, Materials and Structures 58:295 (2025), Eq. 16-26. It
   * extends the cohesive law of Snozzi & Molinari (2013) by frictional sliding and roughness-induced dilatancy.
   *
   * In the local frame of the interface (\f$ x' \f$ normal, \f$ y', z' \f$ tangential) the traction reads
   *
   * \f[ \boldsymbol{t}' = (1-d_c) K_c \begin{bmatrix} u_x - u_{dil} \\ \frac{\beta^2}{\kappa} u_y \\
   *     \frac{\beta^2}{\kappa} u_z \end{bmatrix}
   *     + d_c K_c \begin{bmatrix} u_x - u_{dil} \\ u_y - u_{\mu,y} \\ u_z - u_{\mu,z} \end{bmatrix} H \f]
   *
   * with
   * - the dilatant normal opening \f$ u_{dil} = \bar u_{dil} \tanh( \eta g / \bar u_{dil} \, \| \boldsymbol{u}_T \| )
   *   \f$ (Eq. 21, 26),
   * - the equivalent opening \f$ \alpha_c = \sqrt{ \langle u_x \rangle^2 + \beta^2/\kappa^2 \| \boldsymbol{u}_T \|^2 }
   *   \f$ and its history maximum \f$ \alpha_{max} \f$ (Eq. 18),
   * - the bilinear softening law \f$ q_c(\alpha_{max}) \f$ (Eq. 19) and the damage
   *   \f$ d_c = \min( 1 - q_c / (K_c \alpha_{max}), 1 ) \f$ (Eq. 20),
   * - the contact indicator \f$ H = [ u_x - u_{dil} \leq 0 ] \f$ (Eq. 17),
   * - the plastic sliding \f$ \boldsymbol{u}_\mu \f$ of the penalty-regularised Coulomb friction model
   *   \f$ \| \boldsymbol{t}_\mu \| - \mu \langle - t_x \rangle \leq 0 \f$ with a non-dilatant flow rule (Eq. 22, 23).
   *
   * @par Internal variables and traction
   * #evaluateTraction evaluates Eq. 16 for given #InternalVariables: \f$ \alpha_{max} \f$ (damage), the slip
   * multiplier \f$ \Delta\lambda_\mu \f$ (applied in the direction of the trial friction traction), the norm
   * \f$ u_T \f$ of the tangential jump (dilatancy) and the contact indicator \f$ H \f$. The two updates differ only in
   * where these come from: extrapolated in #computeTraction (Impl-Ex), integrated implicitly for the current jump in
   * #computeTractionExplicit. For the implicit values, the friction traction is that of the radial return.
   *
   * Both updates store the implicit integration for the current jump (#integrateInternalVariables) as the state:
   * friction by a closed-form radial return, and the pressure function \f$ g(t_x) \f$ with the current normal
   * traction, which in contact is the root of a scalar equation (#solveContactTraction). Out of contact the friction
   * limit is zero, so the sliding follows the tangential jump.
   *
   * @par Implicit solvers: Impl-Ex (computeTraction)
   * The internal variables are extrapolated linearly from the two previous converged increments (Oliver, Huespe &
   * Cante, CMAME 197:1865, 2008; Gudžulić & Meschke 2025, Online Resource 1, Appendix C), with
   * \f$ r = \Delta t_{n+1} / \Delta t_n \f$:
   * - the history maximum of the equivalent opening, \f$ \tilde\alpha_{max} = \alpha_{max,n} + r ( \alpha_{max,n}
   *   - \alpha_{max,n-1} ) \f$, which fixes the damage \f$ \tilde d_c = d_c( \tilde\alpha_{max} ) \f$;
   * - the frictional slip multiplier, \f$ \Delta\tilde\lambda_\mu = r \Delta\lambda_{\mu,n} \f$, applied in the
   *   direction of the trial friction traction, \f$ \tilde{\boldsymbol{t}}_\mu = K_c ( \boldsymbol{u}_T -
   *   \boldsymbol{u}_{\mu,n} ) - K_c \Delta\tilde\lambda_\mu \boldsymbol{m} \f$; if the extrapolated sliding exceeds
   *   the trial sliding, \f$ \tilde{\boldsymbol{t}}_\mu = 0 \f$;
   * - the norm of the tangential jump for the dilatancy only, \f$ \tilde u_T = \max( u_{T,n} + r ( u_{T,n} -
   *   u_{T,n-1} ), 0 ) \f$, so that \f$ u_{dil} \f$ does not depend on \f$ \boldsymbol{u}_T \f$ within the increment.
   *
   * The contact indicator \f$ H \f$ is frozen at the last converged increment (@c contact): otherwise a fully
   * damaged interface switches between the normal stiffness 0 and \f$ K_c \f$ during the iterations, and Newton's
   * method cycles. The price is a lag of one increment: closing faces penetrate without resistance, and opening faces
   * transmit the tension \f$ K_c ( u_x - u_{dil} ) \f$. With the internal variables fixed, \f$ t_x \f$ depends on
   * \f$ u_x \f$ only and \f$ \boldsymbol{t}_T \f$ on \f$ \boldsymbol{u}_T \f$ only, so the tangent is symmetric and
   * positive semi-definite, also in softening. The friction limit is not enforced within the increment; the
   * extrapolation error is controlled by the size of the time increment.
   *
   * Nothing is extrapolated (\f$ r = 0 \f$) in the first increment and for \f$ \Delta t = 0 \f$. An increment with
   * \f$ \Delta t = 0 \f$ or \f$ \Delta t < \f$ #minRelativeTimeIncrement \f$ \Delta t_n \f$, e.g. the remainder that
   * closes a step after rounding of the summed time increments, does not shift the history of the extrapolation:
   * its changes are round-off, which the next extrapolation would amplify by \f$ \Delta t_{n+1} / \Delta t \f$.
   *
   * @par Explicit time integration (computeTractionExplicit)
   * The traction of the implicitly integrated internal variables for the current jump, with \f$ H \f$ of the current
   * jump and without a tangent. The state variables are updated as in #computeTraction.
   *
   * @par Material parameters
   * - @b #Kc — penalty stiffness [stress / length] (@c materialProperties[0])
   * - @b #ft — tensile strength (@c materialProperties[1])
   * - @b #GF — total mode I fracture energy (@c materialProperties[2])
   * - @b #Gf — initial mode I fracture energy, \f$ G_f < G_F \f$ (@c materialProperties[3])
   * - @b #ftInt — strength at the kink of the softening law, \f$ 0 < f_{t,int} \leq f_t \f$
   *   (@c materialProperties[4])
   * - @b #beta — ratio of shear to tensile strength (@c materialProperties[5])
   * - @b #kappa — ratio of mode II to mode I fracture energy (@c materialProperties[6])
   * - @b #mu — Coulomb friction coefficient (@c materialProperties[7])
   * - @b #eta — dilatancy coefficient \f$ \tan \psi_{dil} \f$; 0 switches dilatancy off (@c materialProperties[8])
   * - @b #uDilMax — limit \f$ \bar u_{dil} \f$ of the dilatant normal opening (@c materialProperties[9])
   * - @b #txCrit — (optional) negative normal traction at which dilatancy ceases; if omitted or \f$ \geq 0 \f$,
   *   \f$ g = 1 \f$ (@c materialProperties[10])
   * - @b #massPerArea — (optional) mass per unit area of the interface [mass / area], required only for the inertia,
   *   e.g. in explicit time integration (@c materialProperties[11]; requires #txCrit, which may be set to 0)
   *
   * @par Energies
   * \f$ \psi \f$ is the energy of the springs of Eq. 16 for the internal variables of the traction (extrapolated in
   * #computeTraction, implicit in #computeTractionExplicit). The dissipation is the rest of the work, integrated by the
   * trapezoidal rule: \f$ D_{n+1} = D_n + \frac{1}{2} ( \boldsymbol{t}_n + \boldsymbol{t}_{n+1} ) \cdot \Delta
   * [\![ \boldsymbol{u} ]\!] - ( \psi_{n+1} - \psi_n ) \f$. It is exact for linear branches without dilatancy. With
   * Impl-Ex, \f$ D \f$ contains the extrapolation error and can decrease, e.g. after a load reversal.
   *
   * The dilatancy of Eq. 21 is not derived from a potential (\f$ t_x \f$ depends on \f$ \boldsymbol{u}_T \f$, but
   * \f$ \boldsymbol{t}_T \f$ not on \f$ u_x \f$), so for \f$ \eta > 0 \f$ the work depends on the path, already in
   * the elastic range, and \f$ D \f$ can decrease. Shearing a separated crack at \f$ u_x = 0 \f$, for instance, changes
   * \f$ D \f$ by \f$ \frac{1}{2} K_c \eta ( \mu - \eta - \mu^2 \eta ) \| \boldsymbol{u}_T \|^2 \f$ while the
   * dilatancy is linear, which is negative for \f$ \eta > \mu / (1 + \mu^2) \f$.
   *
   * @par State variables
   * - @b alphaMax — history maximum \f$ \alpha_{max} \f$ of the equivalent opening (implicit)
   * - @b alphaMaxOld — \f$ \alpha_{max} \f$ of the increment before
   * - @b slip — plastic sliding \f$ \boldsymbol{u}_\mu \f$ (2 components, implicit)
   * - @b slipMultiplier — slip multiplier \f$ \Delta\lambda_\mu = \| \Delta \boldsymbol{u}_\mu \| \f$ of the last
   *   increment that shifted the history
   * - @b tangentialJump — norm of the tangential jump \f$ u_T \f$ (implicit)
   * - @b tangentialJumpOld — \f$ u_T \f$ of the increment before
   * - @b dTOld — time increment of the last increment that shifted the history of the extrapolation
   * - @b contact — contact indicator \f$ H \f$ of the implicit update, frozen in the next Impl-Ex increment
   * - @b damage — damage \f$ d_c( \alpha_{max} ) \f$ (implicit, for output)
   * - @b alphaC — current equivalent opening \f$ \alpha_c \f$, Eq. 18 (for output, e.g. as crack indicator)
   */
  class CohesiveFrictionalCrack : public Marmot::Materials::MaterialCohesive {

  public:
    // cohesive law
    /** Penalty stiffness (read from @c materialProperties[0]) */
    const double Kc;
    /** Tensile strength (read from @c materialProperties[1]) */
    const double ft;
    /** Total mode I fracture energy (read from @c materialProperties[2]) */
    const double GF;
    /** Initial mode I fracture energy (read from @c materialProperties[3]) */
    const double Gf;
    /** Strength at the kink of the softening law (read from @c materialProperties[4]) */
    const double ftInt;
    /** Ratio of shear to tensile strength (read from @c materialProperties[5]) */
    const double beta;
    /** Ratio of mode II to mode I fracture energy (read from @c materialProperties[6]) */
    const double kappa;

    // friction and dilatancy
    /** Coulomb friction coefficient (read from @c materialProperties[7]) */
    const double mu;
    /** Dilatancy coefficient (read from @c materialProperties[8]) */
    const double eta;
    /** Limit of the dilatant normal opening (read from @c materialProperties[9]) */
    const double uDilMax;
    /** Normal traction at which dilatancy ceases (read from @c materialProperties[10], if provided) */
    const double txCrit;

    // numerics
    /** Time increment, relative to the last one, below which an increment does not shift the extrapolation history */
    static constexpr double minRelativeTimeIncrement = 1e-6;

    // inertia
    /** Mass per unit area of the interface (read from @c materialProperties[11], if provided) */
    const double massPerArea;

    // derived parameters of the bilinear softening law, Eq. 19
    /** Equivalent opening at damage initiation, \f$ f_t / K_c \f$ */
    const double alpha0;
    /** Zero of the initial softening branch */
    const double alpha1;
    /** Equivalent opening at the kink of the softening law */
    const double alphaInt;
    /** Equivalent opening at complete separation */
    const double alphaUlt;

    /**
     * @brief Construct the cohesive-frictional crack model.
     * @param materialProperties Array with parameters: #Kc, #ft, #GF, #Gf, #ftInt, #beta, #kappa, #mu, #eta,
     * #uDilMax, #txCrit (optional), #massPerArea (optional).
     * @param nMaterialProperties Length of @c materialProperties.
     * @param materialNumber Material number.
     * @throws std::invalid_argument if parameters are missing or inadmissible.
     */
    CohesiveFrictionalCrack( const double* materialProperties, int nMaterialProperties, int materialNumber );

    /// Impl-Ex traction with the frozen contact indicator and its algorithmic tangent, for implicit solvers.
    void computeTraction( Response& response, Tangents& tangents, const Increment& increment ) const override;

    /// Traction of the implicitly integrated internal variables for the current jump, for explicit time integration.
    void computeTractionExplicit( Response& response, const Increment& increment ) const override;

    /**
     * @brief Bound on the largest eigenvalue of the tangent over all states, for the critical time step of explicit
     * time integration.
     *
     * The diagonal stiffness is at most \f$ K_c \max( 1, \beta^2 / \kappa ) \f$; softening only lowers it. In the
     * explicit update (#computeTractionExplicit), the dilatancy, \f$ | \partial t_x / \partial \boldsymbol{u}_T | \leq
     * \eta K_c \f$, couples with the friction limit while sliding, \f$ | \partial \boldsymbol{t}_T / \partial u_x |
     * \leq \mu K_c \f$, which raises the eigenvalue by at most \f$ \sqrt{ \mu \eta } K_c \f$. The Impl-Ex tangent of
     * #computeTraction has no such coupling (extrapolated \f$ u_T \f$ and slip multiplier): it is symmetric and bounded
     * by \f$ K_c \max( 1, \beta^2 / \kappa ) \f$.
     *
     * The bound does not depend on the state, so a critical time step evaluated once stays valid when contact,
     * sliding or dilatancy start later.
     *
     * @return \f$ K_c ( \max( 1, \beta^2 / \kappa ) + \sqrt{ \mu \eta } ) \f$
     */
    double getMaximumStiffness( const Response& response, const Eigen::Vector3d& jump ) const override;

    /**
     * @brief Mass per unit area of the interface.
     * @return #massPerArea
     * @throws std::runtime_error if it was not provided, as the density of a bulk material.
     */
    double getMassPerArea( const double* stateVars ) const override;

    void initializeStateLayout()
    {
      stateLayout.add( "alphaMax", 1 );
      stateLayout.add( "alphaMaxOld", 1 );
      stateLayout.add( "slip", 2 );
      stateLayout.add( "slipMultiplier", 1 );
      stateLayout.add( "tangentialJump", 1 );
      stateLayout.add( "tangentialJumpOld", 1 );
      stateLayout.add( "dTOld", 1 );
      stateLayout.add( "contact", 1 );
      stateLayout.add( "damage", 1 );
      stateLayout.add( "alphaC", 1 );
      stateLayout.finalize();
    };

  protected:
    /// Value of the pressure function and its derivative.
    struct PressureFunctionResponse {
      double g;
      double dG_dTx;
    };

    /// Dilatant normal opening and its derivative with respect to the normal jump.
    struct DilatancyResponse {
      double uDil;
      double dUDil_dUx;
    };

    /// Increment of the frictional slip multiplier and plastic sliding of the radial return.
    struct SlipResponse {
      double          slipMultiplier;
      Eigen::Vector2d slip;
    };

    /// Internal variables that determine the traction of an increment, for a history fixed within the increment.
    struct InternalVariables {
      /// history maximum of the equivalent opening, which fixes the damage
      double alphaMax;
      /// increment of the frictional slip multiplier, applied in the direction of the trial friction traction
      double slipMultiplier;
      /// norm of the tangential jump, for the dilatancy
      double tangentialJump;
      /// contact indicator
      bool contact;
    };

    /// Internal variables of the implicit integration for the current jump, and the plastic sliding they imply.
    struct ImplicitUpdate {
      InternalVariables internalVariables;
      Eigen::Vector2d   slip;
    };

    /// Traction, its derivative with respect to the jump, and the recoverable energy.
    struct TractionResponse {
      Eigen::Vector3d traction;
      Eigen::Matrix3d dTraction_dJump;
      double          elasticEnergyDensity;
    };

    /**
     * @brief Extrapolate the internal variables linearly from the two previous increments (Impl-Ex).
     *
     * With \f$ r = \Delta t_{n+1} / \Delta t_n \f$ (0 in the first increment and for \f$ \Delta t = 0 \f$): the
     * history maximum and the tangential jump (not below 0) linearly, the slip multiplier as
     * \f$ r \Delta\lambda_{\mu,n} \f$, and the contact indicator frozen at the last converged increment.
     *
     * @param stateVars state variables at the beginning of the increment
     * @param dT time increment
     */
    InternalVariables extrapolateInternalVariables( const double* stateVars, double dT ) const;

    /**
     * @brief Integrate the internal variables, Eq. 17-23, implicitly for a given jump and the history at the
     * beginning of the increment.
     *
     * Combines #computeDilatancy, #computeAlphaMax and #computeSlip.
     *
     * @param jump total displacement jump in the local frame
     * @param alphaMaxOld history maximum of the equivalent opening at the beginning of the increment
     * @param slipOld plastic sliding at the beginning of the increment
     */
    ImplicitUpdate integrateInternalVariables( const Eigen::Vector3d& jump,
                                               double                 alphaMaxOld,
                                               const Eigen::Vector2d& slipOld ) const;

    /**
     * @brief Evaluate the traction-separation law, Eq. 16, for given internal variables.
     *
     * Only the pressure function follows the current jump. The recoverable energy is that of the springs of Eq. 16,
     * \f$ \psi = \frac{1}{2} t_x ( u_x - u_{dil} ) + \frac{1}{2} (1 - d_c) K_c \frac{\beta^2}{\kappa} \|
     * \boldsymbol{u}_T \|^2 + d_c H \frac{ \| \boldsymbol{t}_\mu \|^2 }{ 2 K_c } \f$. The derivative, with the internal
     * variables fixed, is the Impl-Ex tangent.
     *
     * @param jump total displacement jump in the local frame
     * @param internalVariables internal variables of the increment
     * @param slipOld plastic sliding at the beginning of the increment
     */
    TractionResponse evaluateTraction( const Eigen::Vector3d&   jump,
                                       const InternalVariables& internalVariables,
                                       const Eigen::Vector2d&   slipOld ) const;

    /**
     * @brief Store the implicit update as the state, and shift the extrapolation history unless the increment is
     * negligible (#minRelativeTimeIncrement).
     *
     * @param stateVars state variables of the material point
     * @param implicit implicit update for the jump of the increment
     * @param increment displacement jump and time of the increment
     */
    void storeState( double* stateVars, const ImplicitUpdate& implicit, const Increment& increment ) const;

    /**
     * @brief Update traction and energies of the response; the work of the increment is integrated by the
     * trapezoidal rule, and what is not stored is dissipated.
     */
    static void updateResponse( Response& response, const TractionResponse& traction, const Eigen::Vector3d& dJump );

    /// Damage \f$ d_c = \min( 1 - q_c / (K_c \alpha_{max}), 1 ) \f$ for a given history maximum, Eq. 20.
    double damageFromAlphaMax( double alphaMax ) const;

    /**
     * @brief Dilatant normal opening, Eq. 21, and its derivative with respect to the normal jump; the pressure
     * function of Eq. 26 is evaluated with the normal traction it produces (#solveContactTraction).
     *
     * @param ux normal jump
     * @param uT norm of the tangential jump
     */
    DilatancyResponse computeDilatancy( double ux, double uT ) const;

    /**
     * @brief History maximum of the equivalent opening, Eq. 18, for a given jump.
     *
     * @param jump total displacement jump in the local frame
     * @param alphaMaxOld history maximum of the equivalent opening at the beginning of the increment
     */
    double computeAlphaMax( const Eigen::Vector3d& jump, double alphaMaxOld ) const;

    /**
     * @brief Plastic sliding of the penalty-regularised Coulomb friction, Eq. 22-23, by a closed-form radial return
     * onto the friction limit \f$ \mu K_c \langle -u_{x,c} \rangle \f$, which vanishes out of contact.
     *
     * @param uT tangential jump
     * @param uxContact normal jump less the dilatancy, Eq. 17; negative in contact
     * @param slipOld plastic sliding at the beginning of the increment
     * @return the increment of the slip multiplier \f$ \Delta\lambda_\mu = ( \| \boldsymbol{t}_\mu^{tr} \| - \mu
     * \langle -t_x \rangle ) / K_c \f$ (0 if the faces stick) and the plastic sliding \f$ \boldsymbol{u}_{\mu,n} +
     * \Delta\lambda_\mu \boldsymbol{t}_\mu^{tr} / \| \boldsymbol{t}_\mu^{tr} \| \f$
     */
    SlipResponse computeSlip( const Eigen::Vector2d& uT, double uxContact, const Eigen::Vector2d& slipOld ) const;

    /**
     * @brief Normal traction in contact with the pressure function of Eq. 26 (\f$ t_{x,crit} < 0 \f$), i.e. the root of
     * \f$ t_x - K_c ( u_x - u_{dil}( g( t_x ), \| \boldsymbol{u}_T \| ) ) = 0 \f$, found by bisection in the bracket
     * \f$ [ K_c ( u_x - u_{dil}( 1 ) ), \min( 0, K_c u_x ) ] \f$ to the round-off of
     * \f$ K_c ( |u_x| + u_{dil}( 1 ) ) \f$.
     *
     * @param ux normal jump
     * @param uT norm of the tangential jump
     * @throws std::invalid_argument for a non-finite jump
     * @throws std::runtime_error if the bracket does not shrink to the tolerance (not expected for finite input)
     */
    double solveContactTraction( double ux, double uT ) const;

    /// Pressure function \f$ g(t_x) \f$ of the dilatancy, Eq. 26, and its derivative.
    PressureFunctionResponse pressureFunction( double tx ) const;

    /// Bilinear softening law \f$ q_c(\alpha) \f$, Eq. 19 (elastic branch \f$ K_c \alpha \f$ below \f$ \alpha_0 \f$).
    double softeningFunction( double alpha ) const;
  };

} // namespace Marmot::Materials
