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

#include "Marmot/MarmotStateHelpers.h"
#include <Eigen/Core>
#include <string>

namespace Marmot::Materials {

  /**
   * @brief Abstract base class for cohesive traction-separation laws of zero-thickness interfaces.
   *
   * A cohesive material relates the traction vector \f$ \boldsymbol{t}' \f$ acting on the middle surface of an
   * interface to the displacement jump \f$ [\![ \boldsymbol{u} ]\!]' = \boldsymbol{u}^+ - \boldsymbol{u}^- \f$ across
   * it. Both are expressed in the local frame of the interface, which always has three components:
   *
   * - \f$ x' \f$: the unit normal of the interface, pointing from the bottom to the top side,
   * - \f$ y', z' \f$: two orthonormal directions in the plane of the interface.
   *
   * A positive normal jump opens the interface. Elements embedded in 2D pass \f$ z' = 0 \f$ and use only the first two
   * components of the response.
   *
   * The constitutive relation is formulated in total form: the material receives the total jump at the end of the
   * increment together with its increment, and history enters only through the state variables.
   */
  class MaterialCohesive {

  protected:
    const double* materialProperties;  ///< Pointer to the array of material properties
    const int     nMaterialProperties; ///< Number of material properties

  public:
    const int materialNumber; ///< Integer identifier for this material instance

    /// Layout of the state variables
    Marmot::Utilities::StateLayoutDynamic stateLayout;

    /**
     * @brief Constructs the material with a given set of material properties and an identifier.
     * @param[in] matProperties_       Pointer to the array of material properties.
     * @param[in] nMaterialProperties_ Number of entries in @p matProperties_.
     * @param[in] materialNumber_      Integer identifying this material instance.
     */
    MaterialCohesive( const double* matProperties_, int nMaterialProperties_, int materialNumber_ )
      : materialProperties( matProperties_ ),
        nMaterialProperties( nMaterialProperties_ ),
        materialNumber( materialNumber_ )
    {
    }

    /// Default destructor
    virtual ~MaterialCohesive() = default;

    /// Input of a material point update: the displacement jump and the time of the increment.
    struct Increment {
      Eigen::Vector3d jump;    ///< total displacement jump at the end of the increment, in the local frame
      Eigen::Vector3d dJump;   ///< increment of the displacement jump, in the local frame
      double          timeOld; ///< time at the beginning of the increment
      double          dT;      ///< time increment
    };

    /// Response of the material point: holds the values of the previous increment on entry, the updated ones on exit.
    struct Response {
      Eigen::Vector3d traction;             ///< traction in the local frame
      double*         stateVars;            ///< pointer to the state variables
      double          elasticEnergyDensity; ///< recoverable energy per unit area of the interface
      double          dissipation;          ///< energy dissipated per unit area of the interface, accumulated
    };

    /// Algorithmic tangents of the material point.
    struct Tangents {
      /// \f$ \partial \boldsymbol{t}' / \partial [\![ \boldsymbol{u} ]\!]' \f$ in the local frame
      Eigen::Matrix3d dTraction_dJump = Eigen::Matrix3d::Zero();
    };

    /**
     * @brief Compute the traction and the algorithmic tangent for a given displacement jump.
     *
     * The work done on the interface, \f$ \int \boldsymbol{t}' \cdot d [\![ \boldsymbol{u} ]\!]' \f$ per unit
     * area, is split into the recoverable part @c elasticEnergyDensity and the remainder @c dissipation.
     *
     * @param[in,out] response  traction, state variables and energies, old on entry and updated on exit
     * @param[out]    tangents  algorithmic tangents in the local frame
     * @param[in]     increment displacement jump and time of the increment
     */
    virtual void computeTraction( Response& response, Tangents& tangents, const Increment& increment ) const = 0;

    /**
     * Explicit version of @ref computeTraction for use in explicit time integration schemes.
     * The algorithmic tangent is not needed in explicit schemes and will therefore not be computed.
     *
     * @param[in,out] response  traction, state variables and energies, old on entry and updated on exit
     * @param[in]     increment displacement jump and time of the increment
     *
     * @note The default implementation calls @ref computeTraction and ignores the algorithmic tangent.
     * @note Derived classes may override this method for efficiency reasons.
     */
    virtual void computeTractionExplicit( Response& response, const Increment& increment ) const
    {
      Tangents tangents;
      computeTraction( response, tangents, increment );
    }

    /**
     * @brief Get a view to the state variables.
     * @param stateName Name of the state variable
     * @param stateVars Pointer to the state variable array
     * @return Marmot::Utilities::StateView to access the state variable
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, double* stateVars ) const
    {
      return stateLayout.getStateView( stateVars, stateName );
    }

    /**
     * @brief Get the total number of required state variables.
     * @return Total number of required state variables
     */
    int getNumberOfRequiredStateVars() const { return stateLayout.totalSize(); }

    /**
     * @brief Initialize the state variables at a material point.
     * @param stateVars Pointer to the state variable array
     * @param nStateVars Number of state variables
     *
     * @note The default implementation initializes all state variables to zero.
     */
    virtual void initializeYourself( double* stateVars, int nStateVars )
    {
      for ( int i = 0; i < nStateVars; ++i ) {
        stateVars[i] = 0.0;
      }
    }

    /**
     * @brief Get the maximum stiffness of the material, e.g. for the critical time step of explicit time integration.
     * @param[in] response Current traction, state variables and energies of the material point
     * @param[in] jump     Current displacement jump in the local frame
     * @return Maximum stiffness [stress / length]
     * @details The default implementation computes the algorithmic tangent for a zero jump increment from a copy of
     *          the state and returns the largest diagonal entry of
     *          \f$ \partial \boldsymbol{t}' / \partial [\![ \boldsymbol{u} ]\!]' \f$. A solver that evaluates the
     *          critical time step only once needs a bound over all states the material can reach: a material whose
     *          stiffness may grow in the course of the analysis should override this method.
     */
    virtual double getMaximumStiffness( const Response& response, const Eigen::Vector3d& jump ) const;

    /**
     * @brief Get the mass per unit area of the interface, e.g. for explicit time integration.
     * @param stateVars Pointer to the state variable array
     * @return Mass per unit area of the interface [mass / area]
     * @details The counterpart of the density of a bulk material (MaterialHypoElastic::getDensity) for a
     *          zero-thickness interface, which has no volume. The mass is usually artificial: it gives the penalty
     *          springs between the two faces a stable time increment of their own. As for the density, only the
     *          inertia of the element asks for it, so a material may throw if no mass was given: implicit static
     *          analyses never call it.
     */
    virtual double getMassPerArea( const double* stateVars ) const = 0;
  };

} // namespace Marmot::Materials
