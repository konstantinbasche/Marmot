#include "Marmot/MarmotFiniteElement.h"
#include "Marmot/MarmotGeometryElement.h"
#include "Marmot/MarmotTesting.h"
#include <Eigen/Dense>

using namespace Marmot;
using namespace Marmot::Testing;
using namespace Marmot::FiniteElement;

// ---------------------------------------------------------------------------------------------
// getElementShapeByMetric()
// ---------------------------------------------------------------------------------------------

void testGetElementShapeByMetricAllValidCombinations()
{
  throwExceptionOnFailure( getElementShapeByMetric( 1, 2 ) == Bar2, "nDim=1, nNodes=2 must be Bar2." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 3 ) == Tri3, "nDim=2, nNodes=3 must be Tri3." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 4 ) == Quad4, "nDim=2, nNodes=4 must be Quad4." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 8 ) == Quad8, "nDim=2, nNodes=8 must be Quad8." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 9 ) == Quad9, "nDim=2, nNodes=9 must be Quad9." );
  throwExceptionOnFailure( getElementShapeByMetric( 2, 16 ) == Quad16, "nDim=2, nNodes=16 must be Quad16." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 4 ) == Tetra4, "nDim=3, nNodes=4 must be Tetra4." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 10 ) == Tetra10, "nDim=3, nNodes=10 must be Tetra10." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 8 ) == Hexa8, "nDim=3, nNodes=8 must be Hexa8." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 20 ) == Hexa20, "nDim=3, nNodes=20 must be Hexa20." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 27 ) == Hexa27, "nDim=3, nNodes=27 must be Hexa27." );
  throwExceptionOnFailure( getElementShapeByMetric( 3, 64 ) == Hexa64, "nDim=3, nNodes=64 must be Hexa64." );
}

void testGetElementShapeByMetricThrowsForInvalidCombinations()
{
  const std::vector< std::pair< int, int > > invalidCombinations = {
    { 1, 3 }, // no 3-node 1D shape
    { 2, 5 }, // no 5-node 2D shape
    { 3, 5 }, // no 5-node 3D shape
    { 4, 8 }, // no 4-D shapes at all
  };

  for ( const auto& [nDim, nNodes] : invalidCombinations ) {
    bool threw = false;
    try {
      getElementShapeByMetric( nDim, nNodes );
    }
    catch ( const std::invalid_argument& ) {
      threw = true;
    }
    throwExceptionOnFailure( threw,
                             "getElementShapeByMetric() must throw for nDim=" + std::to_string( nDim ) +
                               ", nNodes=" + std::to_string( nNodes ) + "." );
  }
}

// ---------------------------------------------------------------------------------------------
// NB() / Jacobian(): dynamic (runtime-sized) versions have no callers anywhere in the codebase --
// every element uses the templated (compile-time-sized) overloads declared alongside them in
// MarmotFiniteElement.h instead. Cross-validated against those already-trusted templated versions
// for a known Quad4 configuration, rather than re-deriving expected values independently.
// ---------------------------------------------------------------------------------------------

void testDynamicNBMatchesTemplatedNBForQuad4()
{
  MarmotGeometryElement< 2, 4 > geo;
  const auto                    N = geo.N( Eigen::Vector2d( 0.2, 0.3 ) ); // NSized = Matrix<double,1,4>

  const auto templatedResult = Marmot::FiniteElement::NB< 2, 4 >( N );
  const auto dynamicResult   = Marmot::FiniteElement::NB( Eigen::VectorXd( N.transpose() ), 2 );

  throwExceptionOnFailure( checkIfEqual( Eigen::MatrixXd( dynamicResult ), Eigen::MatrixXd( templatedResult ), 1e-12 ),
                           "The dynamic NB() must match the templated NB<nDim,nNodes>() for the same input." );
}

void testDynamicJacobianMatchesTemplatedJacobianForQuad4()
{
  const std::vector< double > nodeCoordsVec = { 0.0, 0.0, 2.0, 0.0, 3.0, 2.0, 0.0, 1.0 }; // general quadrilateral

  MarmotGeometryElement< 2, 4 > geo;
  geo.assignNodeCoordinates( nodeCoordsVec.data() );

  const Eigen::Vector2d xi( 0.2, 0.3 );
  const auto            dNdXi = geo.dNdXi( xi );

  const auto templatedJ = geo.Jacobian( dNdXi ); // uses the templated Jacobian<2,4> internally

  const Eigen::MatrixXd dynamicDNdXi = dNdXi;
  const Eigen::VectorXd dynamicCoords( Eigen::Map< const Eigen::VectorXd >( nodeCoordsVec.data(), 8 ) );
  const Eigen::MatrixXd dynamicJ = Marmot::FiniteElement::Jacobian( dynamicDNdXi, dynamicCoords );

  throwExceptionOnFailure( checkIfEqual( dynamicJ, Eigen::MatrixXd( templatedJ ), 1e-12 ),
                           "The dynamic Jacobian() must match the templated Jacobian<nDim,nNodes>() for the same "
                           "input." );
}

// ---------------------------------------------------------------------------------------------
// Quadrature::getGaussPointInfo() / getNumGaussPoints(): the branches actually reached by every
// existing element and boundary test all use FullIntegration for Bar3/Tetra4/Tetra10 (or, for
// Bar2/Quad4, never use ReducedIntegration) and never construct a Tetra4/Tetra10 element (there is
// no registered element type using either shape) or pass an unsupported shape at all.
// ---------------------------------------------------------------------------------------------

void testGetGaussPointInfoReducedIntegrationBranches()
{
  using namespace Marmot::FiniteElement::Quadrature;

  throwExceptionOnFailure( getGaussPointInfo( Bar2, ReducedIntegration ).size() == 1,
                           "Bar2 reduced integration must use the 1-point rule." );
  throwExceptionOnFailure( getGaussPointInfo( Bar3, ReducedIntegration ).size() == 2,
                           "Bar3 reduced integration must use the 2-point rule." );
  throwExceptionOnFailure( getGaussPointInfo( Quad4, ReducedIntegration ).size() == 1,
                           "Quad4 reduced integration must use the 1x1-point rule." );
}

void testGetGaussPointInfoTetrahedralShapes()
{
  using namespace Marmot::FiniteElement::Quadrature;

  // Tetra4/Tetra10 quadrature rules don't distinguish integration type; both calls must succeed
  // and return a non-empty rule.
  throwExceptionOnFailure( getGaussPointInfo( Tetra4, FullIntegration ).size() > 0,
                           "Tetra4 must return a non-empty quadrature rule." );
  throwExceptionOnFailure( getGaussPointInfo( Tetra10, ReducedIntegration ).size() > 0,
                           "Tetra10 must return a non-empty quadrature rule." );
}

void testGetGaussPointInfoThrowsForUnsupportedShape()
{
  using namespace Marmot::FiniteElement::Quadrature;

  bool threw = false;
  try {
    getGaussPointInfo( Quad9, FullIntegration );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw, "getGaussPointInfo() must throw for an unsupported shape (Quad9)." );
}

void testTriangleAndNodalRulesIntegratePolynomialsExactly()
{
  using namespace Marmot::FiniteElement::Quadrature;

  // integral over the reference triangle of 1, x and x*y: 1/2, 1/6 and 1/24
  auto integrate = []( const std::vector< QuadraturePointInfo >& rule, auto f ) {
    double result = 0;
    for ( const auto& qp : rule )
      result += qp.weight * f( qp.xi );
    return result;
  };

  for ( auto type : { FullIntegration, ReducedIntegration } ) {
    const auto& rule = getGaussPointInfo( Tri3, type );
    throwExceptionOnFailure( checkIfEqual( integrate( rule, []( const Eigen::VectorXd& ) { return 1.0; } ),
                                           0.5,
                                           1e-15 ),
                             "Tri3 Gauss rules must integrate the area of the reference triangle." );
    throwExceptionOnFailure( checkIfEqual( integrate( rule, []( const Eigen::VectorXd& x ) { return x( 0 ); } ),
                                           1. / 6,
                                           1e-15 ),
                             "Tri3 Gauss rules must integrate linear functions exactly." );
  }
  throwExceptionOnFailure( checkIfEqual( integrate( getGaussPointInfo( Tri3, FullIntegration ),
                                                    []( const Eigen::VectorXd& x ) { return x( 0 ) * x( 1 ); } ),
                                         1. / 24,
                                         1e-15 ),
                           "The 3-point Tri3 rule must integrate quadratic functions exactly." );

  // nodal rules: points at the nodes, in node order, and exact for linear functions
  const std::vector< std::pair< ElementShapes, double > > shapesAndMeasures = { { Bar2, 2.0 },
                                                                                { Quad4, 4.0 },
                                                                                { Tri3, 0.5 } };
  for ( const auto& [shape, measure] : shapesAndMeasures ) {
    const auto& rule = getQuadraturePointInfo( shape, NodalIntegration );
    throwExceptionOnFailure( checkIfEqual( integrate( rule, []( const Eigen::VectorXd& ) { return 1.0; } ),
                                           measure,
                                           1e-15 ),
                             "Nodal rules must integrate the measure of the reference element." );
  }

  for ( size_t i = 0; i < 3; i++ ) {
    const Eigen::Vector2d xi = getQuadraturePointInfo( Tri3, NodalIntegration )[i].xi;
    const auto            N  = FiniteElement::Spatial2D::Tri3::N( xi );
    throwExceptionOnFailure( checkIfEqual( N( i ), 1.0, 1e-15 ) && checkIfEqual( N.sum(), 1.0, 1e-15 ),
                             "The Tri3 nodal rule must list the nodes in node order." );
  }
  for ( size_t i = 0; i < 4; i++ ) {
    const Eigen::Vector2d xi = getQuadraturePointInfo( Quad4, NodalIntegration )[i].xi;
    throwExceptionOnFailure( checkIfEqual( FiniteElement::Spatial2D::Quad4::N( xi )( i ), 1.0, 1e-15 ),
                             "The Quad4 nodal rule must list the nodes in node order." );
  }

  bool threw = false;
  try {
    getQuadraturePointInfo( Quad8, NodalIntegration );
  }
  catch ( const std::invalid_argument& ) {
    threw = true;
  }
  throwExceptionOnFailure( threw,
                           "getQuadraturePointInfo() must throw for nodal integration of an unsupported shape "
                           "(Quad8)." );
}

void testTri3ShapeFunctionDerivatives()
{
  const Eigen::Vector2d xi( 0.2, 0.3 );
  const double          h = 1e-7;
  for ( int j = 0; j < 2; j++ ) {
    Eigen::Vector2d dXi            = Eigen::Vector2d::Zero();
    dXi( j )                       = h;
    const Eigen::RowVector3d dNNum = ( FiniteElement::Spatial2D::Tri3::N( xi + dXi ) -
                                       FiniteElement::Spatial2D::Tri3::N( xi - dXi ) ) /
                                     ( 2 * h );
    throwExceptionOnFailure( checkIfEqual< double >( Eigen::MatrixXd( dNNum ),
                                                     Eigen::MatrixXd(
                                                       FiniteElement::Spatial2D::Tri3::dNdXi( xi ).row( j ) ),
                                                     1e-8 ),
                             "Tri3 shape function derivatives must match central differences." );
  }
  throwExceptionOnFailure( checkIfEqual( FiniteElement::Spatial2D::Tri3::N( xi ).sum(), 1.0, 1e-15 ),
                           "Tri3 shape functions must form a partition of unity." );
}

void testGetQuadraturePointInfoSelectsTheRule()
{
  using namespace Marmot::FiniteElement::Quadrature;

  // Gauss integration is delegated to the very same tables as getGaussPointInfo()
  for ( auto shape : { Bar2, Quad4, Tri3, Hexa8 } )
    for ( auto type : { FullIntegration, ReducedIntegration } )
      throwExceptionOnFailure( &getQuadraturePointInfo( shape, type ) == &getGaussPointInfo( shape, type ),
                               "getQuadraturePointInfo() must return the Gauss rule for Gauss integration." );

  // getGaussPointInfo() must not fall back to a reduced Gauss rule for nodal integration
  auto throwsInvalidArgument = []( auto call ) {
    try {
      call();
    }
    catch ( const std::invalid_argument& ) {
      return true;
    }
    return false;
  };
  throwExceptionOnFailure( throwsInvalidArgument( [] { getGaussPointInfo( Hexa8, NodalIntegration ); } ),
                           "getGaussPointInfo() must throw for nodal integration." );
  throwExceptionOnFailure( throwsInvalidArgument( [] { getNumGaussPoints( Quad4, NodalIntegration ); } ),
                           "getNumGaussPoints() must throw for nodal integration." );
}

void testGetNumGaussPointsMatchesGaussPointInfoSize()
{
  using namespace Marmot::FiniteElement::Quadrature;

  throwExceptionOnFailure( getNumGaussPoints( Hexa8, FullIntegration ) ==
                             static_cast< int >( getGaussPointInfo( Hexa8, FullIntegration ).size() ),
                           "getNumGaussPoints() must match getGaussPointInfo().size()." );
}

int main()
{
  auto tests = std::vector< std::function< void() > >{
    testGetElementShapeByMetricAllValidCombinations,
    testGetElementShapeByMetricThrowsForInvalidCombinations,
    testDynamicNBMatchesTemplatedNBForQuad4,
    testDynamicJacobianMatchesTemplatedJacobianForQuad4,
    testGetGaussPointInfoReducedIntegrationBranches,
    testGetGaussPointInfoTetrahedralShapes,
    testGetGaussPointInfoThrowsForUnsupportedShape,
    testGetNumGaussPointsMatchesGaussPointInfoSize,
    testTriangleAndNodalRulesIntegratePolynomialsExactly,
    testGetQuadraturePointInfoSelectsTheRule,
    testTri3ShapeFunctionDerivatives,
  };

  executeTestsAndCollectExceptions( tests );

  return 0;
}
