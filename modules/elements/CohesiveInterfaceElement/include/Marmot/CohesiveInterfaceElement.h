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

/**
 * @file CohesiveInterfaceElement.h
 * @brief Zero-thickness cohesive interface element for traction-separation laws.
 */
#pragma once

#include "Marmot/MarmotElement.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryInterfaceElement.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialCohesive.h"
#include "Marmot/MarmotMaterialCohesiveFactory.h"
#include "Marmot/MarmotStateVarVectorManager.h"

#include <Eigen/Dense>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Marmot::Discretizations::Elements {

  /**
   * @class CohesiveInterfaceElement
   * @tparam nDim Spatial embedding dimension.
   * @tparam nNodes Number of element nodes, i.e. of both interface sides.
   * @brief Zero-thickness interface element driven by a cohesive traction-separation law.
   *
   * The element is inserted between two faces of the (linear elastic) bulk mesh. Nodes `0 ... nNodes/2 - 1` form
   * the bottom face, nodes `nNodes/2 ... nNodes - 1` the top face, in the same order. The displacement jump
   *
   * \f[ [\![ \boldsymbol{u} ]\!] = \boldsymbol{u}^+ - \boldsymbol{u}^- = \boldsymbol{N}_{[\![u]\!]} \boldsymbol{q} \f]
   *
   * is rotated into the local frame \f$ [\![ \boldsymbol{u} ]\!]' = \boldsymbol{R} [\![ \boldsymbol{u} ]\!] \f$ of
   * the interface, whose first axis is the normal of the mid-surface \f$ (\boldsymbol{x}^- + \boldsymbol{x}^+)/2 \f$
   * (right-hand rule of the node order of a face, so that it points from the bottom towards the top face for a
   * correctly oriented element). The mid-surface also gives the integration measure; for coincident faces it is
   * either face, and for faces that start apart it treats both alike. The traction \f$ \boldsymbol{t}' \f$ of a
   * @ref Marmot::Materials::MaterialCohesive gives the internal force and the tangent
   *
   * \f[ \boldsymbol{P} = \int_\Gamma \boldsymbol{N}_{[\![u]\!]}^T \boldsymbol{R}^T \boldsymbol{t}' \, dA, \qquad
   *     \boldsymbol{K} = \int_\Gamma \boldsymbol{N}_{[\![u]\!]}^T \boldsymbol{R}^T
   *     \frac{\partial \boldsymbol{t}'}{\partial [\![ \boldsymbol{u} ]\!]'} \boldsymbol{R} \boldsymbol{N}_{[\![u]\!]}
   *     \, dA. \f]
   *
   * The kinematics are geometrically linear and the local frame is fixed in the reference configuration. The
   * integration rule of the surface is chosen at registration. The `CI*` elements use the nodal (Newton-Cotes /
   * Lobatto) rule, which avoids the spurious traction oscillations that Gauss integration produces with a high
   * penalty stiffness; the `CI*G` elements use full Gauss integration, e.g. for comparison.
   *
   * The element works with `nDim` local components: the normal and `nDim - 1` in-plane directions. The cohesive
   * material always has three; an element embedded in 2D passes a zero out-of-plane jump, as for plane strain, and
   * keeps the in-plane part of the response. Its out-of-plane thickness is the first element property (default 1).
   *
   * Every quadrature point stores the recoverable and the dissipated energy per unit area of the material; their
   * integral over the interface is the internal energy of the element, i.e. the work done on it.
   */
  template < int nDim, int nNodes >
  class CohesiveInterfaceElement : public Marmot::Discretizations::Elements::Element,
                                   public Marmot::Discretizations::Elements::GeometryInterfaceElement< nDim, nNodes > {

  public:
    /**
     * @brief Section type selector for the interface element.
     */
    enum SectionType {
      Interface,
    };

    static constexpr int sizeLoadVector = nNodes * nDim;

    /// Number of components of the local traction and jump of the cohesive material, for every nDim.
    static constexpr int nMaterialComponents = 3;

    using ParentGeometryElement = Marmot::Discretizations::Elements::GeometryInterfaceElement< nDim, nNodes >;

    using XiSized          = typename ParentGeometryElement::XiSized;
    using VectorDim        = typename ParentGeometryElement::VectorDim;
    using TensorDim        = typename ParentGeometryElement::TensorDim;
    using NJumpMatrixSized = typename ParentGeometryElement::NJumpMatrixSized;

    using RhsSized      = Eigen::Matrix< double, sizeLoadVector, 1 >;
    using KeSizedMatrix = Eigen::Matrix< double, sizeLoadVector, sizeLoadVector >;

    Eigen::Map< const Eigen::VectorXd > elementProperties;

    const int         elLabel;
    const SectionType sectionType;

    /**
     * @brief Container for all per-quadrature-point data.
     */
    struct QuadraturePoint {
      const XiSized xi;
      const double  weight;

      double J0xW;

      /// Rotation into the local frame; its rows are the local axes (normal first) in global coordinates.
      TensorDim R;

      /// Interpolates the displacement jump in the local frame from the element dofs: R * NJumpMatrix( N ), no
      /// derivatives; the rotation keeps the shape of the jump interpolation matrix.
      NJumpMatrixSized NJumpLocal;

      /**
       * @brief Named state-variable manager for cohesive quadrature points.
       *
       * The traction and the jump are stored in the local frame, with nDim components; the energies per unit area
       * of the interface.
       */
      class QPStateVarManager : public Marmot::Utilities::StateVarVectorManager {

        inline const static auto layout = makeLayout( {
          { .name = "traction", .length = nDim },
          { .name = "jump", .length = nDim },
          { .name = "elastic energy density", .length = 1 },
          { .name = "dissipation density", .length = 1 },
          { .name = "begin of material state", .length = 0 },
        } );

      public:
        Eigen::Map< VectorDim > traction;
        Eigen::Map< VectorDim > jump;
        double&                 elasticEnergyDensity;
        double&                 dissipationDensity;

        Eigen::Map< Eigen::VectorXd > materialStateVars;

        /**
         * @brief Number of state variables of the quadrature point itself, without the material state variables.
         */
        static int getNumberOfRequiredStateVarsQuadraturePointOnly() { return layout.nRequiredStateVars; }

        /**
         * @brief Map the state variables of one quadrature point.
         *
         * @param[in] theStateVarVector first state variable of the quadrature point
         * @param[in] nStateVars number of state variables of the quadrature point, including the material state
         */
        QPStateVarManager( double* theStateVarVector, int nStateVars )
          : Marmot::Utilities::StateVarVectorManager( theStateVarVector, layout ),
            traction( &find( "traction" ) ),
            jump( &find( "jump" ) ),
            elasticEnergyDensity( find( "elastic energy density" ) ),
            dissipationDensity( find( "dissipation density" ) ),
            materialStateVars( &find( "begin of material state" ),
                               nStateVars - getNumberOfRequiredStateVarsQuadraturePointOnly() )
        {
        }
      };

      std::unique_ptr< QPStateVarManager > managedStateVars;

      /**
       * @brief Assign state-variable memory to this quadrature point.
       */
      void assignStateVars( double* stateVars, int nStateVars )
      {
        managedStateVars = std::make_unique< QPStateVarManager >( stateVars, nStateVars );
      }

      /**
       * @brief Construct a quadrature-point data container.
       */
      QuadraturePoint( XiSized xi, double weight )
        : xi( xi ), weight( weight ), J0xW( 0.0 ), R( TensorDim::Zero() ), NJumpLocal( NJumpMatrixSized::Zero() )
      {
      }
    };

    std::vector< QuadraturePoint > qps;

    /**
     * @brief The cohesive material, shared by all quadrature points of the element.
     *
     * The material carries no per-quadrature-point data: the state variables are passed in by pointer.
     */
    std::unique_ptr< Marmot::Materials::MaterialCohesive > material;

    /**
     * @brief Construct a cohesive interface element.
     *
     * @param[in] elementID Element label.
     * @param[in] integrationType Integration rule of the surface: nodal, or full or reduced Gauss.
     * @param[in] sectionType Interface section type.
     */
    CohesiveInterfaceElement( int                                                             elementID,
                              Marmot::Discretizations::Elements::Quadrature::IntegrationTypes integrationType,
                              SectionType sectionType = SectionType::Interface );

    int getNumberOfRequiredStateVars();

    std::vector< std::vector< std::string > > getNodeFields();

    std::vector< int > getDofIndicesPermutationPattern();

    /** @brief Number of nodes of the element, i.e. of both interface sides. */
    int getNNodes() { return nNodes; }

    /** @brief Number of spatial dimensions of the element. */
    int getNSpatialDimensions() { return nDim; }

    /** @brief Number of degrees of freedom of the element. */
    int getNDofPerElement() { return sizeLoadVector; }

    /** @brief Result-file geometry keyword, delegated to the geometry element. */
    std::string getElementShape() { return ParentGeometryElement::getElementShape(); }

    void assignStateVars( double* stateVars, int nStateVars );

    /**
     * @brief Assign element-level properties: the out-of-plane thickness for elements embedded in 2D.
     */
    void assignProperty( const Marmot::Discretizations::ElementProperties& marmotElementProperty );

    /**
     * @brief Create the cohesive material of the section through the @ref Marmot::Registry::MaterialCohesiveFactory.
     */
    void assignProperty( const Marmot::Discretizations::MaterialSection& marmotElementProperty );

    void assignNodeCoordinates( const double* coordinates );

    /**
     * @brief Compute the local frames, the jump operators and the integration weights.
     */
    void initializeYourself();

    void setInitialConditions( StateTypes state, const double* values );

    /** @brief Distributed loads are not available for interface elements; throws. */
    void computeDistributedLoad( Marmot::Discretizations::Elements::Element::DistributedLoadTypes loadType,
                                 double*                                                          P,
                                 double*                                                          K,
                                 const int                                                        elementFace,
                                 const double*                                                    load,
                                 const double*                                                    QTotal,
                                 double                                                           time,
                                 double                                                           dT );

    /** @brief Body forces are not available for interface elements; throws. */
    void computeBodyForce( double* P, double* K, const double* load, const double* QTotal, double time, double dT );

    /**
     * @brief Assemble internal force and tangent for one increment.
     *
     * @param QTotal Total displacement vector at the end of the increment.
     * @param dQ Incremental displacement.
     * @param Pe Internal force vector (accumulated), Pe = +Pint.
     * @param Ke Tangent stiffness matrix (accumulated), Ke = +dPe/dQ.
     * @param time Time at the beginning of the increment, forwarded to the material.
     * @param dT Time increment.
     */
    void computeKernels( const double* QTotal, const double* dQ, double* Pe, double* Ke, double time, double dT );

    /**
     * @brief Zero consistent inertia: a zero-thickness interface has no volume, the mass of its nodes is carried by
     * the adjacent bulk elements.
     */
    void computeConsistentInertia( double* M );

    /**
     * @brief Zero lumped inertia: a zero-thickness interface has no volume, the mass of its nodes is carried by the
     * adjacent bulk elements.
     */
    void computeLumpedInertia( double* M );

    /**
     * @brief Internal energy of the element: the recoverable and the dissipated energy of all quadrature points,
     * i.e. the work done on the interface.
     */
    void computeInternalEnergy( double& internalEnergy );

    /**
     * @brief Access a named state view at a quadrature point.
     *
     * The quadrature point provides "traction" and "jump" (local frame, nDim components), "elastic energy density"
     * and "dissipation density" (per unit area); all other names are forwarded to the material.
     */
    Marmot::Utilities::StateView getStateView( const std::string& stateName, int qpNumber )
    {
      const auto& qp = qps[qpNumber];

      if ( qp.managedStateVars->contains( stateName ) ) {
        return qp.managedStateVars->getStateView( stateName );
      }

      return material->getStateView( stateName, qp.managedStateVars->materialStateVars.data() );
    }

    /**
     * @brief Coordinates of the element center (on the mid-surface).
     */
    std::vector< double > getCoordinatesAtCenter();

    /**
     * @brief Coordinates of all quadrature points (on the mid-surface).
     */
    std::vector< std::vector< double > > getCoordinatesAtQuadraturePoints();

    int getNumberOfQuadraturePoints();
  };

  /**
   * @name Template implementation (header-only)
   */
  ///@{

  template < int nDim, int nNodes >
  CohesiveInterfaceElement< nDim, nNodes >::CohesiveInterfaceElement(
    int                                                             elementID,
    Marmot::Discretizations::Elements::Quadrature::IntegrationTypes integrationType,
    SectionType                                                     sectionType )
    : ParentGeometryElement(),
      elementProperties( Eigen::Map< const Eigen::VectorXd >( nullptr, 0 ) ),
      elLabel( elementID ),
      sectionType( sectionType )
  {
    for ( const auto& qpInfo :
          Marmot::Discretizations::Elements::Quadrature::getQuadraturePointInfo( this->shape, integrationType ) )
      qps.emplace_back( XiSized( qpInfo.xi ), qpInfo.weight );
  }

  template < int nDim, int nNodes >
  int CohesiveInterfaceElement< nDim, nNodes >::getNumberOfRequiredStateVars()
  {
    const int nMaterialStateVars = material ? material->getNumberOfRequiredStateVars() : 0;

    return ( QuadraturePoint::QPStateVarManager::getNumberOfRequiredStateVarsQuadraturePointOnly() +
             nMaterialStateVars ) *
           qps.size();
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< std::string > > CohesiveInterfaceElement< nDim, nNodes >::getNodeFields()
  {
    return std::vector< std::vector< std::string > >( nNodes, { "displacement" } );
  }

  template < int nDim, int nNodes >
  std::vector< int > CohesiveInterfaceElement< nDim, nNodes >::getDofIndicesPermutationPattern()
  {
    std::vector< int > permutationPattern( nNodes * nDim );
    for ( int i = 0; i < nNodes * nDim; i++ )
      permutationPattern[i] = i;

    return permutationPattern;
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::assignStateVars( double* stateVars, int nStateVars )
  {
    const int nQpStateVars = nStateVars / qps.size();

    for ( size_t i = 0; i < qps.size(); i++ )
      qps[i].assignStateVars( stateVars + ( i * nQpStateVars ), nQpStateVars );
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::assignProperty(
    const Marmot::Discretizations::ElementProperties& elementPropertiesInfo )
  {
    new ( &elementProperties ) Eigen::Map< const Eigen::VectorXd >( elementPropertiesInfo.elementProperties,
                                                                    elementPropertiesInfo.nElementProperties );
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::assignProperty(
    const Marmot::Discretizations::MaterialSection& section )
  {
    material = std::unique_ptr< Marmot::Materials::MaterialCohesive >(
      Marmot::Registry::MaterialCohesiveFactory::createMaterial( section.materialName,
                                                                 section.materialProperties,
                                                                 section.nMaterialProperties,
                                                                 elLabel ) );
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::assignNodeCoordinates( const double* coordinates )
  {
    ParentGeometryElement::assignNodeCoordinates( coordinates );
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::initializeYourself()
  {
    const double thickness = ( nDim == 2 && elementProperties.size() > 0 ) ? elementProperties[0] : 1.0;

    for ( QuadraturePoint& qp : qps ) {
      const auto N     = this->N( qp.xi );
      const auto dNdXi = this->dNdXi( qp.xi );
      const auto J     = this->surfaceJacobian( dNdXi, ParentGeometryElement::InterfaceSurface::Mid );

      qp.R          = this->localFrame( J );
      qp.NJumpLocal = qp.R * this->NJumpMatrix( N );
      qp.J0xW       = qp.weight * this->sqrtDetMetric( J ) * thickness;
    }
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::computeKernels( const double* QTotal_,
                                                                 const double* dQ_,
                                                                 double*       Pe_,
                                                                 double*       Ke_,
                                                                 double        time,
                                                                 double        dT )
  {
    Eigen::Map< const RhsSized > QTotal( QTotal_ );
    Eigen::Map< const RhsSized > dQ( dQ_ );
    Eigen::Map< RhsSized >       Pe( Pe_ );
    Eigen::Map< KeSizedMatrix >  Ke( Ke_ );

    using MaterialVector = Eigen::Matrix< double, nMaterialComponents, 1 >;

    // the components of the material beyond nDim, i.e. the out-of-plane jump in 2D, are zero
    const auto toMaterial = []( const VectorDim& v ) {
      MaterialVector v3          = MaterialVector::Zero();
      v3.template head< nDim >() = v;
      return v3;
    };

    for ( QuadraturePoint& qp : qps ) {
      auto& qpState = *qp.managedStateVars;

      const VectorDim displacementJumpLocal          = qp.NJumpLocal * QTotal;
      const VectorDim displacementJumpIncrementLocal = qp.NJumpLocal * dQ;

      Marmot::Materials::MaterialCohesive::Response        response{ toMaterial( qpState.traction ),
                                                              qpState.materialStateVars.data(),
                                                              qpState.elasticEnergyDensity,
                                                              qpState.dissipationDensity };
      Marmot::Materials::MaterialCohesive::Tangents        tangents;
      const Marmot::Materials::MaterialCohesive::Increment increment{ toMaterial( displacementJumpLocal ),
                                                                      toMaterial( displacementJumpIncrementLocal ),
                                                                      time,
                                                                      dT };

      material->computeTraction( response, tangents, increment );

      const VectorDim traction = response.traction.template head< nDim >();
      const TensorDim C        = tangents.dTraction_dJump.template topLeftCorner< nDim, nDim >();

      qpState.traction             = traction;
      qpState.jump                 = displacementJumpLocal;
      qpState.elasticEnergyDensity = response.elasticEnergyDensity;
      qpState.dissipationDensity   = response.dissipation;

      Pe += qp.NJumpLocal.transpose() * traction * qp.J0xW;
      Ke += qp.NJumpLocal.transpose() * C * qp.NJumpLocal * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::setInitialConditions( StateTypes state, const double* values )
  {
    switch ( state ) {
    case Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization: {
      for ( QuadraturePoint& qp : qps ) {
        material->initializeYourself( qp.managedStateVars->materialStateVars.data(),
                                      qp.managedStateVars->materialStateVars.size() );
      }
      break;
    }

    default:
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << ": invalid initial condition for CohesiveInterfaceElement" );
    }
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::computeDistributedLoad(
    Marmot::Discretizations::Elements::Element::DistributedLoadTypes,
    double*,
    double*,
    const int,
    const double*,
    const double*,
    double,
    double )
  {
    throw std::invalid_argument( Marmot::Utilities::MakeString() << __PRETTY_FUNCTION__
                                                                 << ": distributed loads are not implemented for "
                                                                    "CohesiveInterfaceElement." );
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::computeBodyForce( double*,
                                                                   double*,
                                                                   const double*,
                                                                   const double*,
                                                                   double,
                                                                   double )
  {
    throw std::invalid_argument( Marmot::Utilities::MakeString()
                                 << __PRETTY_FUNCTION__
                                 << ": body forces are not implemented for CohesiveInterfaceElement." );
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::computeConsistentInertia( double* M )
  {
    Eigen::Map< KeSizedMatrix >( M ).setZero();
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::computeLumpedInertia( double* M )
  {
    Eigen::Map< RhsSized >( M ).setZero();
  }

  template < int nDim, int nNodes >
  void CohesiveInterfaceElement< nDim, nNodes >::computeInternalEnergy( double& internalEnergy )
  {
    internalEnergy = 0.0;
    for ( const QuadraturePoint& qp : qps ) {
      const auto& qpState = *qp.managedStateVars;
      internalEnergy += ( qpState.elasticEnergyDensity + qpState.dissipationDensity ) * qp.J0xW;
    }
  }

  template < int nDim, int nNodes >
  std::vector< double > CohesiveInterfaceElement< nDim, nNodes >::getCoordinatesAtCenter()
  {
    std::vector< double > coords( nDim );

    // the parametric center: the centroid of a triangle, the origin of the bi-unit line and quadrilateral
    XiSized centerXi = XiSized::Zero();
    if constexpr ( ParentGeometryElement::getInterfaceElementShape() ==
                   Marmot::Discretizations::Elements::ElementShapes::Tri3 )
      centerXi.setConstant( 1. / 3 );

    Eigen::Map< VectorDim >( coords.data() ) = this->NMatrix( this->N( centerXi ) ) *
                                               this->getSideCoordinates( ParentGeometryElement::InterfaceSurface::Mid );

    return coords;
  }

  template < int nDim, int nNodes >
  std::vector< std::vector< double > > CohesiveInterfaceElement< nDim, nNodes >::getCoordinatesAtQuadraturePoints()
  {
    std::vector< std::vector< double > > listedCoords;

    const auto xSide = this->getSideCoordinates( ParentGeometryElement::InterfaceSurface::Mid );

    for ( const auto& qp : qps ) {
      std::vector< double > coords( nDim );
      Eigen::Map< VectorDim >( coords.data() ) = this->NMatrix( this->N( qp.xi ) ) * xSide;
      listedCoords.push_back( coords );
    }

    return listedCoords;
  }

  template < int nDim, int nNodes >
  int CohesiveInterfaceElement< nDim, nNodes >::getNumberOfQuadraturePoints()
  {
    return static_cast< int >( qps.size() );
  }

  ///@}

} // namespace Marmot::Discretizations::Elements
