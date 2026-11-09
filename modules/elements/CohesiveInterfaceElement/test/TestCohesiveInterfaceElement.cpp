#include "Marmot/CohesiveInterfaceElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotElementProperty.h"
#include "Marmot/MarmotJournal.h"
#include "Marmot/MarmotMaterialCohesiveFactory.h"
#include "Marmot/MarmotTesting.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::Discretizations::Elements;

namespace {

  /// Linear penalty law t' = diag(kn, ks, ks) [[u]]', to check the element kinematics exactly; the optional third
  /// parameter is the mass per area, which, as the density of a bulk material, throws if asked for but not given.
  class LinearPenaltyCohesive : public Marmot::Materials::MaterialCohesive {
  public:
    LinearPenaltyCohesive( const double* props, int nProps, int number )
      : Marmot::Materials::MaterialCohesive( props, nProps, number )
    {
      stateLayout.add( "calls", 1 );
      stateLayout.finalize();
    }

    void computeTraction( Response& response, Tangents& tangents, const Increment& increment ) const override
    {
      tangents.dTraction_dJump = Eigen::Vector3d( materialProperties[0], materialProperties[1], materialProperties[1] )
                                   .asDiagonal();
      response.traction             = tangents.dTraction_dJump * increment.jump;
      response.elasticEnergyDensity = 0.5 * response.traction.dot( increment.jump );
      stateLayout.getAs< double& >( response.stateVars, "calls" ) += 1;
    }

    double getMassPerArea( const double* ) const override
    {
      if ( nMaterialProperties <= 2 )
        throw std::runtime_error( "TESTLINEARPENALTYCOHESIVE: mass per area not specified." );
      return materialProperties[2];
    }
  };

  const static bool
    linearPenaltyIsRegistered = Marmot::Registry::MaterialCohesiveFactory::registerMaterial< LinearPenaltyCohesive >(
      "TESTLINEARPENALTYCOHESIVE" );

  /// The cohesive-frictional law with the parameters of Gudzulic & Meschke (2025), Table 2, but a softer penalty.
  const std::vector< double > cohesiveFrictionalProperties =
    { 3.4e4, 2.2, 0.14, 0.056, 0.55, 2.0, 10.0, 0.7, 0.2, 16.0 };

  /**
   * The cohesive-frictional law is a separate module: the tests that need a nonlinear, dissipative law run only if
   * it is built and registered.
   */
  bool cohesiveFrictionalCrackIsAvailable()
  {
    try {
      std::unique_ptr< Marmot::Materials::MaterialCohesive > material(
        Marmot::Registry::MaterialCohesiveFactory::createMaterial( "COHESIVEFRICTIONALCRACK",
                                                                   cohesiveFrictionalProperties.data(),
                                                                   cohesiveFrictionalProperties.size(),
                                                                   1 ) );
      return material != nullptr;
    }
    catch ( const std::exception& ) {
      return false;
    }
  }

  template < class A, class B >
  bool isClose( const A& a, const B& b, double tol )
  {
    return ( a - b ).template lpNorm< Eigen::Infinity >() <= tol;
  }

  /**
   * An element geometry with coincident bottom and top faces, unless top offsets are given, together with the
   * expected integration measure and local frame of the coincident faces.
   */
  template < int nDim_, int nNodes_ >
  struct Geometry {
    static constexpr int nDim       = nDim_;
    static constexpr int nNodes     = nNodes_;
    static constexpr int nFaceNodes = nNodes / 2;
    static constexpr int nDof       = nDim * nNodes;

    using Element   = CohesiveInterfaceElement< nDim, nNodes >;
    using Vector    = Eigen::Matrix< double, nDim, 1 >;
    using DofVector = Eigen::Matrix< double, nDof, 1 >;
    using DofMatrix = Eigen::Matrix< double, nDof, nDof >;

    std::string                      name;
    std::string                      shape;
    std::array< Vector, nFaceNodes > faceNodes;
    std::vector< double >            elementProperties; ///< the out-of-plane thickness in 2D
    double                           measure;           ///< area, or length times thickness in 2D
    Vector                           normal;            ///< the first local axis
    Vector                           firstEdge;         ///< the second local axis at the first node
    std::vector< Vector >            topOffsets;        ///< top face node minus bottom face node; empty: coincident

    std::vector< double > coordinates() const
    {
      std::vector< double > x( nDof );
      for ( int side = 0; side < 2; side++ )
        for ( int A = 0; A < nFaceNodes; A++ )
          for ( int i = 0; i < nDim; i++ )
            x[( side * nFaceNodes + A ) * nDim + i] = faceNodes[A]( i ) +
                                                      ( side == 1 && !topOffsets.empty() ? topOffsets[A]( i ) : 0.0 );
      return x;
    }

    /// Nodal displacements with a uniform jump: the top face moves by `jump` relative to the bottom face.
    DofVector uniformJump( const Vector& jump, const Vector& translation ) const
    {
      DofVector Q;
      for ( int A = 0; A < nFaceNodes; A++ ) {
        Q.template segment< nDim >( nDim * A )                  = translation;
        Q.template segment< nDim >( nDim * ( nFaceNodes + A ) ) = translation + jump;
      }
      return Q;
    }
  };

  /// A line in the plane with an out-of-plane thickness.
  Geometry< 2, 4 > inclinedLine()
  {
    Geometry< 2, 4 > g;
    g.name                 = "CILINE2";
    g.shape                = "quad4";
    g.faceNodes            = { Eigen::Vector2d( 0.1, -0.2 ), Eigen::Vector2d( 2.1, 0.7 ) };
    const double thickness = 2.5;
    g.elementProperties    = { thickness };
    g.measure              = ( g.faceNodes[1] - g.faceNodes[0] ).norm() * thickness;
    g.firstEdge            = ( g.faceNodes[1] - g.faceNodes[0] ).normalized();
    g.normal               = Eigen::Vector2d( -g.firstEdge( 1 ), g.firstEdge( 0 ) );
    return g;
  }

  /// A triangle in a general plane.
  Geometry< 3, 6 > inclinedTriangle()
  {
    Geometry< 3, 6 > g;
    g.name                   = "CITRI3";
    g.shape                  = "penta6";
    g.faceNodes              = { Eigen::Vector3d( 0.1, -0.2, 0.3 ),
                                 Eigen::Vector3d( 2.1, 0.4, -0.5 ),
                                 Eigen::Vector3d( 0.4, 1.7, 0.9 ) };
    const Eigen::Vector3d e1 = g.faceNodes[1] - g.faceNodes[0], e2 = g.faceNodes[2] - g.faceNodes[0];
    g.measure   = 0.5 * e1.cross( e2 ).norm();
    g.normal    = e1.cross( e2 ).normalized();
    g.firstEdge = e1.normalized();
    return g;
  }

  /// A plane, distorted (non-parallelogram) quadrilateral in a general plane.
  Geometry< 3, 8 > inclinedQuadrilateral()
  {
    const Geometry< 3, 6 > triangle = inclinedTriangle();
    const Eigen::Vector3d  e1       = triangle.firstEdge;
    const Eigen::Vector3d  e2       = triangle.normal.cross( e1 );

    // counter-clockwise in the plane coordinates (e1, e2)
    const std::array< Eigen::Vector2d, 4 > p = { Eigen::Vector2d( 0.0, 0.0 ),
                                                 Eigen::Vector2d( 2.0, 0.0 ),
                                                 Eigen::Vector2d( 2.3, 1.5 ),
                                                 Eigen::Vector2d( -0.2, 1.1 ) };

    Geometry< 3, 8 > g;
    g.name    = "CIQUAD4";
    g.shape   = "hexa8";
    g.measure = 0;
    for ( int A = 0; A < 4; A++ ) {
      g.faceNodes[A] = triangle.faceNodes[0] + p[A]( 0 ) * e1 + p[A]( 1 ) * e2;
      g.measure += 0.5 * ( p[A]( 0 ) * p[( A + 1 ) % 4]( 1 ) - p[( A + 1 ) % 4]( 0 ) * p[A]( 1 ) );
    }
    g.normal    = triangle.normal;
    g.firstEdge = e1;
    return g;
  }

  /// An initialized element together with the memory it works on.
  template < class G >
  struct ElementFixture {
    using Element = typename G::Element;

    std::vector< double >      coordinates;
    std::vector< double >      elementProperties;
    std::vector< double >      materialProperties;
    std::unique_ptr< Element > element;
    std::vector< double >      stateVars;

    ElementFixture( const G&                                                      geometry,
                    const std::string&                                            materialName,
                    std::vector< double >                                         props,
                    Marmot::Discretizations::Elements::Quadrature::QuadratureRule quadratureRule =
                      Marmot::Discretizations::Elements::Quadrature::newtonCotesExactToDegree( 1 ) )
      : coordinates( geometry.coordinates() ),
        elementProperties( geometry.elementProperties ),
        materialProperties( std::move( props ) ),
        element( std::make_unique< Element >( 7, quadratureRule, Element::SectionType::Interface ) )
    {
      element->assignNodeCoordinates( coordinates.data() );
      if ( !elementProperties.empty() )
        element->assignProperty(
          Marmot::Discretizations::ElementProperties( elementProperties.data(), elementProperties.size() ) );
      element->assignProperty( Marmot::Discretizations::MaterialSection( materialName,
                                                                         materialProperties.data(),
                                                                         materialProperties.size() ) );
      stateVars.resize( element->getNumberOfRequiredStateVars() );
      element->assignStateVars( stateVars.data(), stateVars.size() );
      element->initializeYourself();
      element->setInitialConditions( Marmot::Discretizations::Elements::Element::MarmotMaterialInitialization,
                                     nullptr );
    }

    /// Internal force and tangent for a total displacement, starting from the current history.
    std::pair< typename G::DofVector, typename G::DofMatrix > compute( const typename G::DofVector& Q,
                                                                       const typename G::DofVector& dQ )
    {
      typename G::DofVector P = G::DofVector::Zero();
      typename G::DofMatrix K = G::DofMatrix::Zero();
      element->computeKernels( Q.data(), dQ.data(), P.data(), K.data(), 0.0, 1.0 );
      return { P, K };
    }
  };

  template < class G >
  void testFactory( const G& g )
  {
    using namespace Marmot::Discretizations::Elements::Quadrature;

    // the CI* elements are integrated with the nodal rule, the CI*G elements with full Gauss integration
    const auto shape = G::Element::ParentGeometryElement::getInterfaceElementShape();
    for ( const auto& [name, nQuadraturePoints] :
          { std::pair{ g.name, G::nFaceNodes },
            std::pair{ g.name + "G", int( getQuadraturePoints( shape, gaussExactToDegree( 2 ) ).size() ) } } ) {

      std::unique_ptr< Marmot::Discretizations::Elements::Element > element(
        Marmot::Registry::ElementFactory::createElement( name, 1 ) );
      throwExceptionOnFailure( element != nullptr, name + " is not registered with the element factory." );
      throwExceptionOnFailure( element->getElementShape() == g.shape,
                               name + " must be written as the EnSight cell " + g.shape + "." );
      throwExceptionOnFailure( element->getNNodes() == G::nNodes && element->getNDofPerElement() == G::nDof &&
                                 element->getNSpatialDimensions() == G::nDim,
                               name + " descriptors are wrong." );
      throwExceptionOnFailure( int( element->getNodeFields().size() ) == G::nNodes &&
                                 element->getNodeFields()[G::nNodes - 1][0] == "displacement",
                               name + " node fields are wrong." );
      throwExceptionOnFailure( dynamic_cast< typename G::Element& >( *element ).getNumberOfQuadraturePoints() ==
                                 nQuadraturePoints,
                               Marmot::Utilities::MakeString()
                                 << name << " must have " << nQuadraturePoints << " integration points." );
    }
  }

  template < class G >
  void testNodalIntegrationAndLocalFrame( const G& g )
  {
    using Matrix = Eigen::Matrix< double, G::nDim, G::nDim >;

    ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0 } );

    throwExceptionOnFailure( f.element->getNumberOfQuadraturePoints() == G::nFaceNodes,
                             g.name + " must have one integration point per face node." );

    double weightSum = 0;
    for ( const auto& qp : f.element->qps ) {
      weightSum += qp.J0xW;
      throwExceptionOnFailure( isClose( qp.R * qp.R.transpose(), Matrix::Identity(), 1e-14 ),
                               g.name + ": the local frame must be orthonormal." );
      throwExceptionOnFailure( isClose( qp.R.row( 0 ).transpose(), g.normal, 1e-14 ),
                               g.name + ": the first local axis must be the normal of the face." );
    }
    throwExceptionOnFailure( isClose( f.element->qps[0].R.row( 1 ).transpose(), g.firstEdge, 1e-14 ),
                             g.name + ": the second local axis must follow the first edge." );
    throwExceptionOnFailure( checkIfEqual( weightSum, g.measure, 1e-14 ),
                             Marmot::Utilities::MakeString()
                               << g.name << ": the weights sum to " << weightSum << " instead of the measure "
                               << g.measure << " of the face." );

    const auto qpCoords = f.element->getCoordinatesAtQuadraturePoints();
    for ( int A = 0; A < G::nFaceNodes; A++ )
      throwExceptionOnFailure( isClose( Eigen::Map< const typename G::Vector >( qpCoords[A].data() ),
                                        g.faceNodes[A],
                                        1e-14 ),
                               g.name + ": the integration points must coincide with the nodes, in node order." );

    typename G::Vector centroid = G::Vector::Zero();
    for ( const auto& x : g.faceNodes )
      centroid += x / G::nFaceNodes;
    const auto center = f.element->getCoordinatesAtCenter();
    throwExceptionOnFailure( isClose( Eigen::Map< const typename G::Vector >( center.data() ), centroid, 1e-14 ),
                             g.name + ": the element center must be the mean of the face nodes." );
  }

  /**
   * With faces that start apart, the geometry is that of the mid-surface: the element must have the same local
   * frames, weights and coordinates as an element whose coincident faces lie on the mid-surface.
   */
  template < class G >
  void testFacesApartUseMidSurface( const G& g )
  {
    G gApart = g, gMid = g;
    for ( int A = 0; A < G::nFaceNodes; A++ ) {
      // a gap that varies over the face and a tangential misalignment, so that the faces are not parallel
      const typename G::Vector offset = 0.1 * ( 1 + A ) * g.normal + 0.03 * ( A - 1 ) * g.firstEdge;
      gApart.topOffsets.push_back( offset );
      gMid.faceNodes[A] += 0.5 * offset;
    }

    ElementFixture fApart( gApart, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0 } );
    ElementFixture fMid( gMid, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0 } );

    for ( int q = 0; q < G::nFaceNodes; q++ ) {
      const auto& qpApart = fApart.element->qps[q];
      const auto& qpMid   = fMid.element->qps[q];
      throwExceptionOnFailure( isClose( qpApart.R, qpMid.R, 1e-14 ),
                               g.name + ": with faces apart, the local frame must be that of the mid-surface." );
      throwExceptionOnFailure( checkIfEqual( qpApart.J0xW, qpMid.J0xW, 1e-14 ),
                               g.name + ": with faces apart, the weights must be those of the mid-surface." );
      throwExceptionOnFailure( isClose( Eigen::Map< const typename G::Vector >(
                                          fApart.element->getCoordinatesAtQuadraturePoints()[q].data() ),
                                        gMid.faceNodes[q],
                                        1e-14 ),
                               g.name + ": with faces apart, the integration points must lie on the mid-surface." );
    }

    const auto centerApart = fApart.element->getCoordinatesAtCenter();
    const auto centerMid   = fMid.element->getCoordinatesAtCenter();
    throwExceptionOnFailure( isClose( Eigen::Map< const typename G::Vector >( centerApart.data() ),
                                      Eigen::Map< const typename G::Vector >( centerMid.data() ),
                                      1e-14 ),
                             g.name + ": with faces apart, the element center must lie on the mid-surface." );
  }

  template < class G >
  void testRigidTranslationGivesZeroForce( const G& g )
  {
    ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { 1e6, 3e5 } );

    const auto Q      = g.uniformJump( G::Vector::Zero(), G::Vector::LinSpaced( 0.3, -1.2 ) );
    const auto [P, K] = f.compute( Q, Q );
    throwExceptionOnFailure( P.norm() < 1e-9, g.name + ": a rigid translation must not produce internal forces." );
    throwExceptionOnFailure( ( K * Q ).norm() < 1e-9,
                             g.name + ": rigid translations must lie in the kernel of the stiffness." );
  }

  template < class G >
  void testUniformJumpGivesPenaltyForce( const G& g )
  {
    const double   kn = 1e6, ks = 3e5, delta = 2e-4;
    ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { kn, ks } );

    // normal opening plus sliding along the first edge, superposed on a translation
    const typename G::Vector jump = delta * g.normal + 0.5 * delta * g.firstEdge;
    const auto               Q    = g.uniformJump( jump, G::Vector::LinSpaced( 1.0, 3.0 ) );
    const auto [P, K]             = f.compute( Q, Q );

    // the law is isotropic in the plane of the interface, so the global traction does not depend on the in-plane axes
    const typename G::Vector traction = kn * delta * g.normal + ks * 0.5 * delta * g.firstEdge;
    double                   forceSum = 0;
    for ( int A = 0; A < G::nFaceNodes; A++ ) {
      const typename G::Vector nodalForce = traction * f.element->qps[A].J0xW;
      throwExceptionOnFailure( isClose( P.template segment< G::nDim >( G::nDim * ( G::nFaceNodes + A ) ),
                                        nodalForce,
                                        1e-9 ),
                               Marmot::Utilities::MakeString()
                                 << g.name << ": wrong force at top node " << A << ": "
                                 << P.template segment< G::nDim >( G::nDim * ( G::nFaceNodes + A ) ).transpose() );
      throwExceptionOnFailure( isClose( P.template segment< G::nDim >( G::nDim * A ), -nodalForce, 1e-9 ),
                               g.name + ": the bottom nodes must carry the opposite force." );
      forceSum += P.template segment< G::nDim >( G::nDim * ( G::nFaceNodes + A ) ).dot( g.normal );
    }
    throwExceptionOnFailure( checkIfEqual( forceSum, kn * delta * g.measure, 1e-9 ),
                             g.name + ": the normal force must be the normal traction times the measure of the face." );

    // the stored traction and jump are local, with nDim components
    const auto jumpView     = f.element->getStateView( "jump", 0 );
    const auto tractionView = f.element->getStateView( "traction", 0 );
    throwExceptionOnFailure( jumpView.stateSize == G::nDim && tractionView.stateSize == G::nDim,
                             g.name + ": the stored jump and traction must have nDim components." );
    throwExceptionOnFailure( checkIfEqual( jumpView.stateLocation[0], delta, 1e-15 ) &&
                               checkIfEqual( jumpView.stateLocation[1], 0.5 * delta, 1e-15 ),
                             g.name + ": the stored jump must be the local jump." );
    throwExceptionOnFailure( checkIfEqual( tractionView.stateLocation[0], kn * delta, 1e-9 ) &&
                               checkIfEqual( tractionView.stateLocation[1], ks * 0.5 * delta, 1e-9 ),
                             g.name + ": the stored traction must be the local traction." );

    // the material state is reachable through the element and follows the jump directly
    throwExceptionOnFailure( checkIfEqual( f.element->getStateView( "calls", G::nFaceNodes - 1 ).stateLocation[0],
                                           1.0,
                                           0.0 ),
                             g.name + ": material state views must be forwarded to the material." );
    throwExceptionOnFailure( G::Element::QuadraturePoint::QPStateVarManager::
                                 getNumberOfRequiredStateVarsQuadraturePointOnly() == 2 * G::nDim + 2,
                             g.name + ": the material state must follow the traction, the jump and the energies." );
  }

  template < class G >
  void testStiffnessMatchesNumericalDifferentiation( const G& g )
  {
    using namespace Marmot::Discretizations::Elements::Quadrature;

    for ( auto quadratureRule : { newtonCotesExactToDegree( 1 ), gaussExactToDegree( 2 ) } ) {
      ElementFixture f( g, "COHESIVEFRICTIONALCRACK", cohesiveFrictionalProperties, quadratureRule );

      // a mixed-mode jump that varies over the face
      typename G::DofVector Q = G::DofVector::Zero();
      for ( int A = 0; A < G::nFaceNodes; A++ )
        for ( int i = 0; i < G::nDim; i++ ) {
          const int k = G::nDim * ( G::nFaceNodes + A ) + i;
          Q( k )      = 1e-2 * ( 0.6 + 0.2 * A ) * g.normal( i ) + 2e-3 * std::sin( 1.3 * k + 0.4 );
        }
      f.compute( Q, Q );

      const std::vector< double > stateVarsOld = f.stateVars;
      typename G::DofVector       QNext        = 1.05 * Q;
      QNext( 0 ) += 2e-4;

      const auto [P, K] = f.compute( QNext, QNext - Q );

      typename G::DofMatrix KNum;
      const double          h = 1e-8;
      for ( int j = 0; j < G::nDof; j++ ) {
        typename G::DofVector dQ = G::DofVector::Zero();
        dQ( j )                  = h;

        f.stateVars   = stateVarsOld;
        const auto Pp = f.compute( QNext + dQ, QNext + dQ - Q ).first;

        f.stateVars   = stateVarsOld;
        const auto Pm = f.compute( QNext - dQ, QNext - dQ - Q ).first;

        KNum.col( j ) = ( Pp - Pm ) / ( 2 * h );
      }

      const double err = ( K - KNum ).template lpNorm< Eigen::Infinity >() / K.template lpNorm< Eigen::Infinity >();
      throwExceptionOnFailure( err < 1e-6,
                               Marmot::Utilities::MakeString()
                                 << g.name << " ("
                                 << ( quadratureRule == newtonCotesExactToDegree( 1 ) ? "nodal" : "Gauss" ) << " rule"
                                 << "): the stiffness differs from numerical differentiation: " << err );
    }
  }

  /**
   * The internal energy of the element, recoverable plus dissipated, is the work of the element forces along the load
   * path, integrated by the trapezoidal rule like the work of the tractions in the material. The path loads
   * elastically, softens, unloads, and closes and slides the interface under compression.
   */
  template < class G >
  void testInternalEnergyIsTheWorkOfTheElementForces( const G& g )
  {
    ElementFixture f( g, "COHESIVEFRICTIONALCRACK", cohesiveFrictionalProperties );

    // the normal and tangential jump amplitudes of the waypoints, and the number of increments to reach them
    const double alpha0 = cohesiveFrictionalProperties[1] / cohesiveFrictionalProperties[0];
    const std::vector< std::tuple< double, double, int > > waypoints = {
      { 0.5 * alpha0, 0.2 * alpha0, 1 }, // elastic
      { 2e-2, 5e-3, 20 },                // softening
      { 1e-2, 2.5e-3, 5 },               // unloading
      { -1e-4, 2e-2, 10 },               // closing and sliding
      { -1e-4, 5e-2, 10 },               // sliding
    };

    // a jump that varies over the face: the top nodes move, the bottom nodes are fixed
    const auto displacement = [&]( double normal, double tangential ) {
      typename G::DofVector Q = G::DofVector::Zero();
      for ( int A = 0; A < G::nFaceNodes; A++ )
        Q.template segment< G::nDim >( G::nDim * ( G::nFaceNodes + A ) ) = normal * ( 1 + 0.2 * A ) * g.normal +
                                                                           tangential * ( 1 + 0.1 * A ) * g.firstEdge;
      return Q;
    };

    double                internalEnergy, work = 0, normalOld = 0, tangentialOld = 0;
    typename G::DofVector QOld = G::DofVector::Zero(), POld = G::DofVector::Zero();
    for ( size_t w = 0; w < waypoints.size(); w++ ) {
      const auto [normal, tangential, nIncrements] = waypoints[w];
      for ( int i = 1; i <= nIncrements; i++ ) {
        const double s = double( i ) / nIncrements;
        const auto   Q = displacement( normalOld + s * ( normal - normalOld ),
                                     tangentialOld + s * ( tangential - tangentialOld ) );
        const auto   P = f.compute( Q, Q - QOld ).first;
        work += 0.5 * ( POld + P ).dot( Q - QOld );
        QOld = Q;
        POld = P;
      }
      normalOld     = normal;
      tangentialOld = tangential;

      f.element->computeInternalEnergy( internalEnergy );
      throwExceptionOnFailure( checkIfEqual( internalEnergy, work, 1e-10 * work ),
                               Marmot::Utilities::MakeString()
                                 << g.name << ", waypoint " << w << ": the internal energy " << internalEnergy
                                 << " differs from the work " << work << "." );
    }

    for ( int q = 0; q < G::nFaceNodes; q++ ) {
      const auto dissipation = f.element->getStateView( "dissipation density", q );
      throwExceptionOnFailure( dissipation.stateSize == 1 && dissipation.stateLocation[0] > 0,
                               g.name + ": cracking and sliding must dissipate at every integration point." );
    }
  }

  /**
   * Full Gauss integration integrates the same face as the nodal rule, but couples the node pairs: with the nodal
   * rule each integration point sees the jump of one node pair only, so the stiffness has no blocks between
   * different node pairs.
   */
  template < class G >
  void testGaussIntegrationCouplesNodePairs( const G& g )
  {
    using namespace Marmot::Discretizations::Elements::Quadrature;

    const double   kn = 1e6, ks = 3e5, delta = 2e-4;
    ElementFixture fNodal( g, "TESTLINEARPENALTYCOHESIVE", { kn, ks }, newtonCotesExactToDegree( 1 ) );
    ElementFixture fGauss( g, "TESTLINEARPENALTYCOHESIVE", { kn, ks }, gaussExactToDegree( 2 ) );

    double weightSum = 0;
    for ( const auto& qp : fGauss.element->qps )
      weightSum += qp.J0xW;
    throwExceptionOnFailure( checkIfEqual( weightSum, g.measure, 1e-14 ),
                             g.name + ": the Gauss weights must sum to the measure of the face." );

    // a uniform jump: both rules give the normal traction times the measure as the resultant on the top face
    const typename G::Vector jump = delta * g.normal;
    const auto               Q    = g.uniformJump( jump, G::Vector::Zero() );
    const auto [PNodal, KNodal]   = fNodal.compute( Q, Q );
    const auto [PGauss, KGauss]   = fGauss.compute( Q, Q );

    for ( const auto& P : { PNodal, PGauss } ) {
      typename G::Vector resultant = G::Vector::Zero();
      for ( int A = 0; A < G::nFaceNodes; A++ )
        resultant += P.template segment< G::nDim >( G::nDim * ( G::nFaceNodes + A ) );
      throwExceptionOnFailure( isClose( resultant, kn * delta * g.measure * g.normal, 1e-9 ),
                               g.name + ": the resultant force must be the traction times the measure of the face." );
    }

    // the block between the top nodes 0 and 1
    const auto block = []( const typename G::DofMatrix& K ) {
      return K.template block< G::nDim, G::nDim >( G::nDim * G::nFaceNodes, G::nDim * ( G::nFaceNodes + 1 ) );
    };
    throwExceptionOnFailure( block( KNodal ).norm() == 0.0,
                             g.name + ": with nodal integration, different node pairs must not be coupled." );
    throwExceptionOnFailure( block( KGauss ).norm() > 1e-3 * KGauss.norm(),
                             g.name + ": with Gauss integration, different node pairs must be coupled." );
  }

  template < class G >
  void testUnsupportedOperationsThrow( const G& g )
  {
    ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0 } );

    auto expectThrow = [&]( const std::function< void() >& call, const std::string& what ) {
      bool thrown = false;
      try {
        call();
      }
      catch ( const std::exception& ) {
        thrown = true;
      }
      throwExceptionOnFailure( thrown, g.name + ": " + what + " must throw." );
    };

    std::vector< double > dummy( G::nDof * G::nDof );
    expectThrow( [&] { f.element->computeBodyForce( dummy.data(), dummy.data(), dummy.data(), dummy.data(), 0, 1 ); },
                 "computeBodyForce" );
    expectThrow(
      [&] {
        f.element->setInitialConditions( Marmot::Discretizations::Elements::Element::GeostaticStress, dummy.data() );
      },
      "an unsupported initial condition" );
    expectThrow(
      [&] {
        const std::vector< double > props = { 1.0 };
        f.element->assignProperty(
          Marmot::Discretizations::MaterialSection( "NOTACOHESIVEMATERIAL", props.data(), 1 ) );
      },
      "an unknown cohesive material" );
  }

  /**
   * The inertia asks the material for its mass per area, as a bulk element asks for the density: without one it
   * throws, so only analyses with inertia need it. A massless material gives zero inertia, overwriting what the
   * output held before.
   */
  template < class G >
  void testInertiaNeedsMassPerArea( const G& g )
  {
    ElementFixture fWithoutMass( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0 } );
    for ( const bool lumped : { false, true } ) {
      bool thrown = false;
      try {
        typename G::DofMatrix M;
        if ( lumped )
          fWithoutMass.element->computeLumpedInertia( M.data() );
        else
          fWithoutMass.element->computeConsistentInertia( M.data() );
      }
      catch ( const std::runtime_error& ) {
        thrown = true;
      }
      throwExceptionOnFailure( thrown, g.name + ": the inertia without a mass per area of the material must throw." );
    }

    ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0, 0.0 } );

    typename G::DofMatrix M = G::DofMatrix::Constant( 1.0 );
    f.element->computeConsistentInertia( M.data() );
    throwExceptionOnFailure( M.isZero( 0.0 ), g.name + ": the consistent inertia must be zero." );

    typename G::DofVector MLumped = G::DofVector::Constant( 1.0 );
    f.element->computeLumpedInertia( MLumped.data() );
    throwExceptionOnFailure( MLumped.isZero( 0.0 ), g.name + ": the lumped inertia must be zero." );
  }

  /**
   * With a material that has a mass per area, half of the mass sits on each face: the mass of each face and direction
   * is half the mass per area times the measure of the face. For the Gauss rule the lumped inertia is the row sum of
   * the consistent one; the nodal rule lumps with its own points, which on a distorted face differ from the row sums
   * of the exact consistent inertia.
   */
  template < class G >
  void testArtificialMass( const G& g )
  {
    using namespace Marmot::Discretizations::Elements::Quadrature;

    const double massPerArea = 2.5e-3;
    for ( const auto quadratureRule : { newtonCotesExactToDegree( 1 ), gaussExactToDegree( 2 ) } ) {
      ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0, massPerArea }, quadratureRule );

      typename G::DofMatrix M = G::DofMatrix::Constant( 1.0 );
      f.element->computeConsistentInertia( M.data() );
      typename G::DofVector MLumped = G::DofVector::Constant( 1.0 );
      f.element->computeLumpedInertia( MLumped.data() );

      if ( quadratureRule == gaussExactToDegree( 2 ) )
        throwExceptionOnFailure( isClose( MLumped, M.rowwise().sum(), 1e-15 ),
                                 g.name + ": the lumped inertia must be the row sums of the consistent inertia." );
      throwExceptionOnFailure( isClose( M, M.transpose(), 0.0 ), g.name + ": the consistent inertia is symmetric." );
      for ( int side = 0; side < 2; side++ )
        for ( int i = 0; i < G::nDim; i++ ) {
          double faceMass = 0;
          for ( int A = 0; A < G::nFaceNodes; A++ )
            faceMass += MLumped( ( side * G::nFaceNodes + A ) * G::nDim + i );
          throwExceptionOnFailure( checkIfEqual( faceMass, 0.5 * massPerArea * g.measure, 1e-15 ),
                                   g.name + ": each face must carry half of the mass." );
        }
    }
  }

  /**
   * The consistent inertia is integrated with the Gauss rule of degree 2p, whatever the rule of the element: the
   * nodal and the Gauss element have the same, exact consistent mass, which couples the nodes of a face. For linear
   * faces it is known in closed form: mu / 2 * L / 6 * ( 1 + delta_AB ) for a line of length L, and
   * mu / 2 * A / 12 * ( 1 + delta_AB ) for a triangle of area A, in each direction.
   */
  template < class G >
  void testConsistentMassIsExact( const G& g )
  {
    using namespace Marmot::Discretizations::Elements::Quadrature;

    const double massPerArea = 2.5e-3;

    ElementFixture fNodal( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0, massPerArea }, newtonCotesExactToDegree( 1 ) );
    ElementFixture fGauss( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0, massPerArea }, gaussExactToDegree( 2 ) );
    typename G::DofMatrix MNodal, MGauss;
    fNodal.element->computeConsistentInertia( MNodal.data() );
    fGauss.element->computeConsistentInertia( MGauss.data() );

    throwExceptionOnFailure( isClose( MNodal, MGauss, 0.0 ),
                             g.name + ": the nodal and the Gauss element must have the same consistent inertia." );

    constexpr int nFaceDofs = G::nFaceNodes * G::nDim;
    throwExceptionOnFailure( MNodal.template topRightCorner< nFaceDofs, nFaceDofs >().isZero( 0.0 ) &&
                               MNodal.template bottomLeftCorner< nFaceDofs, nFaceDofs >().isZero( 0.0 ),
                             g.name + ": the consistent inertia must not couple the two faces." );

    // the exact mass of a linear line or triangle; the distorted quadrilateral is covered by the Gauss element
    if constexpr ( G::nFaceNodes == 2 || G::nFaceNodes == 3 ) {
      const double unit = G::nFaceNodes == 2 ? g.measure / 6 : g.measure / 12;

      Eigen::Matrix< double, nFaceDofs, nFaceDofs > MFace = Eigen::Matrix< double, nFaceDofs, nFaceDofs >::Zero();
      for ( int A = 0; A < G::nFaceNodes; A++ )
        for ( int B = 0; B < G::nFaceNodes; B++ )
          for ( int i = 0; i < G::nDim; i++ )
            MFace( A * G::nDim + i, B * G::nDim + i ) = 0.5 * massPerArea * unit * ( A == B ? 2.0 : 1.0 );

      throwExceptionOnFailure( isClose( MNodal.template topLeftCorner< nFaceDofs, nFaceDofs >(), MFace, 1e-17 ) &&
                                 isClose( MNodal.template bottomRightCorner< nFaceDofs, nFaceDofs >(), MFace, 1e-17 ),
                               g.name + ": the consistent inertia of a linear face must be the exact one." );
    }
  }

  /**
   * The critical time step bounds the stable increment of the element from its own mass: it is at most 2 / omega_max
   * of the generalized eigenproblem K v = omega^2 M v with the lumped mass, and equal to it for the nodal rule, where
   * it is sqrt( mass per area / maximum stiffness ).
   */
  template < class G >
  void testCriticalTimeStep( const G& g )
  {
    using namespace Marmot::Discretizations::Elements::Quadrature;

    const double kn = 1e6, ks = 3e5, massPerArea = 2.5e-3;
    for ( const auto quadratureRule : { newtonCotesExactToDegree( 1 ), gaussExactToDegree( 2 ) } ) {
      ElementFixture fMassless( g, "TESTLINEARPENALTYCOHESIVE", { kn, ks, 0.0 }, quadratureRule );

      bool thrown = false;
      try {
        double dt;
        fMassless.element->computeCriticalTimeStepForExplicitDynamics( dt, nullptr );
      }
      catch ( const std::invalid_argument& ) {
        thrown = true;
      }
      throwExceptionOnFailure( thrown, g.name + ": the critical time step of an element without mass must throw." );

      ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { kn, ks, massPerArea }, quadratureRule );
      double         dt;
      f.element->computeCriticalTimeStepForExplicitDynamics( dt, nullptr );

      // the default MaterialCohesive::getMaximumStiffness() evaluates the material on a copy of the state
      throwExceptionOnFailure( f.element->getStateView( "calls", 0 ).stateLocation[0] == 0.0,
                               g.name + ": the critical time step must not change the state of the material." );

      const typename G::DofVector Q = G::DofVector::Zero();
      const auto                  K = f.compute( Q, Q ).second;
      typename G::DofVector       M;
      f.element->computeLumpedInertia( M.data() );
      const typename G::DofVector MInvSqrt = M.cwiseSqrt().cwiseInverse();
      const double                omegaMax = std::sqrt(
        Eigen::SelfAdjointEigenSolver< typename G::DofMatrix >( MInvSqrt.asDiagonal() * K * MInvSqrt.asDiagonal() )
          .eigenvalues()
          .maxCoeff() );

      throwExceptionOnFailure( dt <= 2.0 / omegaMax * ( 1 + 1e-12 ),
                               Marmot::Utilities::MakeString()
                                 << g.name << " ("
                                 << ( quadratureRule == newtonCotesExactToDegree( 1 ) ? "nodal" : "Gauss" ) << " rule"
                                 << "): the critical time step " << dt << " exceeds the stable increment "
                                 << 2.0 / omegaMax << "." );
      if ( quadratureRule == newtonCotesExactToDegree( 1 ) )
        throwExceptionOnFailure( checkIfEqual( dt, std::sqrt( massPerArea / kn ), 1e-12 * dt ) &&
                                   checkIfEqual( dt, 2.0 / omegaMax, 1e-12 * dt ),
                                 g.name + ": with the nodal rule the critical time step is sqrt( mu / K_max )." );
    }
  }

  /// The mass per area belongs to the material, as the density of a bulk element: the element has no named
  /// properties, and the base element rejects any name, the former element property "mass per area" included.
  template < class G >
  void testNamedPropertiesAreRejected( const G& g )
  {
    ElementFixture f( g, "TESTLINEARPENALTYCOHESIVE", { 1.0, 1.0 } );

    throwExceptionOnFailure( f.element->getPropertyNames().empty(), g.name + ": the element has no named properties." );

    // through the base, as a solver assigns it: the element's own overloads of assignProperty hide the named one
    Marmot::Discretizations::Elements::Element& element = *f.element;
    for ( const std::string name : { "mass per area", "not a property" } ) {
      const double value  = 1.0;
      bool         thrown = false;
      try {
        element.assignProperty( name, &value, 1 );
      }
      catch ( const std::invalid_argument& ) {
        thrown = true;
      }
      throwExceptionOnFailure( thrown, g.name + ": the named property '" + name + "' must throw." );
    }
  }

  /**
   * The explicit kernels are the implicit ones without the tangent: along a load path that opens, softens, closes
   * and slides the interface, forces and states are identical.
   */
  template < class G >
  void testExplicitKernelsMatchImplicit( const G& g )
  {
    const double alpha0 = cohesiveFrictionalProperties[1] / cohesiveFrictionalProperties[0];

    std::vector< std::pair< std::string, std::vector< double > > > materials = {
      { "TESTLINEARPENALTYCOHESIVE", { 1e6, 3e5 } } };
    if ( cohesiveFrictionalCrackIsAvailable() )
      materials.emplace_back( "COHESIVEFRICTIONALCRACK", cohesiveFrictionalProperties );

    for ( const auto& [materialName, props] : materials ) {
      ElementFixture fImplicit( g, materialName, props );
      ElementFixture fExplicit( g, materialName, props );

      typename G::DofVector QOld = G::DofVector::Zero();
      for ( const auto& [normal, tangential] : std::vector< std::pair< double, double > >{
              { 0.5 * alpha0, 0.2 * alpha0 }, // elastic
              { 2e-2, 5e-3 },                 // softening
              { 1e-2, 2.5e-3 },               // unloading
              { -1e-4, 2e-2 },                // closing and sliding
            } ) {
        typename G::DofVector Q = G::DofVector::Zero();
        for ( int A = 0; A < G::nFaceNodes; A++ )
          Q.template segment< G::nDim >( G::nDim * ( G::nFaceNodes + A ) ) = normal * ( 1 + 0.2 * A ) * g.normal +
                                                                             tangential * ( 1 + 0.1 * A ) * g.firstEdge;
        const typename G::DofVector dQ = Q - QOld;

        const auto            PImplicit = fImplicit.compute( Q, dQ ).first;
        typename G::DofVector PExplicit = G::DofVector::Zero();
        fExplicit.element->computeKernelsExplicit( Q.data(), dQ.data(), PExplicit.data(), 0.0, 1.0 );
        QOld = Q;

        throwExceptionOnFailure( PExplicit == PImplicit && fExplicit.stateVars == fImplicit.stateVars,
                                 g.name + " with " + materialName +
                                   ": the explicit kernels must give the forces and states of the implicit ones." );
      }
    }
  }

  template < class G >
  std::function< void() > forElement( void ( *test )( const G& ), G ( *geometry )() )
  {
    return [=] { test( geometry() ); };
  }

  /// All tests of one element type.
  template < class G >
  std::vector< std::function< void() > > testsFor( G ( *geometry )() )
  {
    std::vector< std::function< void() > > tests = { forElement( testFactory< G >, geometry ),
                                                     forElement( testNodalIntegrationAndLocalFrame< G >, geometry ),
                                                     forElement( testFacesApartUseMidSurface< G >, geometry ),
                                                     forElement( testRigidTranslationGivesZeroForce< G >, geometry ),
                                                     forElement( testUniformJumpGivesPenaltyForce< G >, geometry ),
                                                     forElement( testGaussIntegrationCouplesNodePairs< G >, geometry ),
                                                     forElement( testUnsupportedOperationsThrow< G >, geometry ),
                                                     forElement( testInertiaNeedsMassPerArea< G >, geometry ),
                                                     forElement( testArtificialMass< G >, geometry ),
                                                     forElement( testConsistentMassIsExact< G >, geometry ),
                                                     forElement( testCriticalTimeStep< G >, geometry ),
                                                     forElement( testNamedPropertiesAreRejected< G >, geometry ),
                                                     forElement( testExplicitKernelsMatchImplicit< G >, geometry ) };
    if ( cohesiveFrictionalCrackIsAvailable() ) {
      tests.push_back( forElement( testStiffnessMatchesNumericalDifferentiation< G >, geometry ) );
      tests.push_back( forElement( testInternalEnergyIsTheWorkOfTheElementForces< G >, geometry ) );
    }
    return tests;
  }

} // namespace

/**
 * A CILINE2 with thickness t must respond exactly like a CIQUAD4 that extrudes the line by t out of plane, with all
 * out-of-plane displacements zero: the 2D element passes a zero out-of-plane jump to the cohesive material. The
 * load path opens the interface beyond the strength, closes it again and slides it under compression.
 */
void TestLineMatchesExtrudedQuadrilateral()
{
  const Geometry< 2, 4 > line      = inclinedLine();
  const double           thickness = line.elementProperties[0];

  // the quadrilateral a, b, b - t e_z, a - t e_z has the same normal and first edge as the line
  Geometry< 3, 8 > strip = inclinedQuadrilateral();
  const auto       lift  = []( const Eigen::Vector2d& x, double z ) { return Eigen::Vector3d( x( 0 ), x( 1 ), z ); };
  strip.faceNodes        = { lift( line.faceNodes[0], 0 ),
                             lift( line.faceNodes[1], 0 ),
                             lift( line.faceNodes[1], -thickness ),
                             lift( line.faceNodes[0], -thickness ) };

  // Q3D = T Q2D: the nodes of the strip follow the line node they are extruded from
  const std::array< int, 8 >     lineNode = { 0, 1, 1, 0, 2, 3, 3, 2 };
  Eigen::Matrix< double, 24, 8 > T        = Eigen::Matrix< double, 24, 8 >::Zero();
  for ( int m = 0; m < 8; m++ )
    for ( int i = 0; i < 2; i++ )
      T( 3 * m + i, 2 * lineNode[m] + i ) = 1;

  ElementFixture fLine( line, "COHESIVEFRICTIONALCRACK", cohesiveFrictionalProperties );
  ElementFixture fStrip( strip, "COHESIVEFRICTIONALCRACK", cohesiveFrictionalProperties );

  // the local jumps (normal, tangential) at the two ends of the line, which differ to vary the jump over the face
  const std::vector< std::array< Eigen::Vector2d, 2 > > path = {
    { Eigen::Vector2d( 2e-5, 0 ), Eigen::Vector2d( 3e-5, 1e-5 ) },      // elastic
    { Eigen::Vector2d( 2e-2, 5e-3 ), Eigen::Vector2d( 3e-2, 1e-2 ) },   // softening, open
    { Eigen::Vector2d( -1e-4, 2e-2 ), Eigen::Vector2d( -2e-4, 3e-2 ) }, // closed, sliding
    { Eigen::Vector2d( -1e-4, 1e-1 ), Eigen::Vector2d( -2e-4, 1e-1 ) }, // closed, sliding further
  };

  Eigen::Matrix< double, 8, 1 > QOld = Eigen::Matrix< double, 8, 1 >::Zero();
  for ( size_t step = 0; step < path.size(); step++ ) {
    Eigen::Matrix< double, 8, 1 > Q = Eigen::Matrix< double, 8, 1 >::Zero();
    for ( int A = 0; A < 2; A++ )
      Q.segment< 2 >( 2 * ( 2 + A ) ) = path[step][A]( 0 ) * line.normal + path[step][A]( 1 ) * line.firstEdge;

    const auto [PLine, KLine]   = fLine.compute( Q, Q - QOld );
    const auto [PStrip, KStrip] = fStrip.compute( T * Q, T * ( Q - QOld ) );
    QOld                        = Q;

    const double scale = std::max( PLine.lpNorm< Eigen::Infinity >(), 1e-30 );
    throwExceptionOnFailure( isClose( PLine, T.transpose() * PStrip, 1e-12 * scale ),
                             Marmot::Utilities::MakeString()
                               << "Step " << step << ": the CILINE2 force differs from the extruded CIQUAD4:\n"
                               << PLine.transpose() << "\n"
                               << ( T.transpose() * PStrip ).transpose() );
    throwExceptionOnFailure( isClose( KLine,
                                      Eigen::Matrix< double, 8, 8 >( T.transpose() * KStrip * T ),
                                      1e-12 * KLine.lpNorm< Eigen::Infinity >() ),
                             Marmot::Utilities::MakeString()
                               << "Step " << step << ": the CILINE2 stiffness differs from the extruded CIQUAD4." );
    for ( int m = 0; m < 8; m++ )
      throwExceptionOnFailure( std::abs( PStrip( 3 * m + 2 ) ) <= 1e-12 * scale,
                               Marmot::Utilities::MakeString()
                                 << "Step " << step << ": the extruded CIQUAD4 has an out-of-plane force." );
  }
}

int main()
{
  std::vector< std::function< void() > > tests;
  for ( auto testsOfElement :
        { testsFor( inclinedLine ), testsFor( inclinedTriangle ), testsFor( inclinedQuadrilateral ) } )
    tests.insert( tests.end(), testsOfElement.begin(), testsOfElement.end() );
  if ( cohesiveFrictionalCrackIsAvailable() )
    tests.push_back( TestLineMatchesExtrudedQuadrilateral );
  else
    std::cout << "COHESIVEFRICTIONALCRACK is not built: the tests with a dissipative cohesive law are skipped.\n";

  executeTestsAndCollectExceptions( tests );

  return 0;
}
