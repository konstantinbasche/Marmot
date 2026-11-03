#include "Marmot/CohesiveInterfaceElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace {

  using namespace Marmot::Discretizations::Elements;
  using namespace Marmot::Discretizations::Elements::Quadrature;

  /**
   * @brief Create the factory function of a cohesive interface element type.
   *
   * @tparam T element type
   * @tparam integrationType integration rule of the surface of the created elements
   * @return a function that creates an element of type T from an element number
   */
  template < class T, IntegrationTypes integrationType >
  Marmot::Registry::ElementFactory::elementFactoryFunction makeCohesiveFactoryFunction()
  {
    return []( int elementID ) -> Marmot::Discretizations::Elements::Element* {
      return new T( elementID, integrationType, T::SectionType::Interface );
    };
  }

  // nodal integration, which avoids traction oscillations with a high penalty stiffness

  const static bool CILINE2_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CILINE2", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 2, 4 >, NodalIntegration >() );

  const static bool CITRI3_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CITRI3", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 6 >, NodalIntegration >() );

  const static bool CIQUAD4_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CIQUAD4", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 8 >, NodalIntegration >() );

  // full Gauss integration, e.g. for comparison; reduced Gauss integration is not registered, since a single point
  // per face leaves zero-energy modes of the jump

  const static bool CILINE2G_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CILINE2G", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 2, 4 >, FullIntegration >() );

  const static bool CITRI3G_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CITRI3G", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 6 >, FullIntegration >() );

  const static bool CIQUAD4G_isRegistered = Marmot::Registry::ElementFactory::
    registerElement( "CIQUAD4G", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 8 >, FullIntegration >() );

} // namespace
