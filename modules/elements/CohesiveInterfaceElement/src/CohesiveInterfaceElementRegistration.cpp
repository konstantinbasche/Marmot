#include "Marmot/CohesiveInterfaceElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace {

  using namespace Marmot::Discretizations::Elements;

  /**
   * @brief Create the factory function of a cohesive interface element type.
   *
   * @tparam T element type
   * @tparam quadratureRule quadrature rule of the surface of the created elements
   * @return a function that creates an element of type T from an element number
   */
  template < class T, Quadrature::QuadratureRule quadratureRule >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeCohesiveFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, quadratureRule, T::SectionType::Interface );
    };
  }

  // nodal integration, which avoids traction oscillations with a high penalty stiffness

  const static bool CILINE2_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CILINE2",
                     makeCohesiveFactoryFunction< CohesiveInterfaceElement< 2, 4 >,
                                                  Quadrature::newtonCotesExactToDegree( 1 ) >() );

  const static bool CITRI3_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CITRI3",
                     makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 6 >,
                                                  Quadrature::newtonCotesExactToDegree( 1 ) >() );

  const static bool CIQUAD4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CIQUAD4",
                     makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 8 >,
                                                  Quadrature::newtonCotesExactToDegree( 1 ) >() );

  // Gauss integration, e.g. for comparison; a single point per face (gaussExactToDegree( 1 )) is not registered, since
  // it leaves zero-energy modes of the jump

  const static bool CILINE2G_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CILINE2G",
                     makeCohesiveFactoryFunction< CohesiveInterfaceElement< 2, 4 >,
                                                  Quadrature::gaussExactToDegree( 2 ) >() );

  const static bool CITRI3G_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CITRI3G",
                     makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 6 >,
                                                  Quadrature::gaussExactToDegree( 2 ) >() );

  const static bool CIQUAD4G_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CIQUAD4G",
                     makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 8 >,
                                                  Quadrature::gaussExactToDegree( 2 ) >() );

} // namespace
