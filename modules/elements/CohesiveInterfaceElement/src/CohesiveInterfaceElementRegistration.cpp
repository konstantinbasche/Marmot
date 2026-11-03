#include "Marmot/CohesiveInterfaceElement.h"
#include "Marmot/MarmotElementFactory.h"
#include "Marmot/MarmotFiniteElement.h"

namespace Marmot::Elements::Registration {

  using namespace MarmotLibrary;
  using namespace Marmot::FiniteElement::Quadrature;

  /**
   * @brief Create the factory function of a cohesive interface element type.
   *
   * @tparam T element type
   * @tparam integrationType integration rule of the surface of the created elements
   * @return a function that creates an element of type T from an element number
   */
  template < class T, IntegrationTypes integrationType >
  MarmotElementFactory::elementFactoryFunction makeCohesiveFactoryFunction()
  {
    return
      []( int elementID ) -> MarmotElement* { return new T( elementID, integrationType, T::SectionType::Interface ); };
  }

  // nodal integration, which avoids traction oscillations with a high penalty stiffness

  const static bool CILINE2_isRegistered = MarmotElementFactory::
    registerElement( "CILINE2", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 2, 4 >, NodalIntegration >() );

  const static bool CITRI3_isRegistered = MarmotElementFactory::
    registerElement( "CITRI3", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 6 >, NodalIntegration >() );

  const static bool CIQUAD4_isRegistered = MarmotElementFactory::
    registerElement( "CIQUAD4", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 8 >, NodalIntegration >() );

  // full Gauss integration, e.g. for comparison; reduced Gauss integration is not registered, since a single point
  // per face leaves zero-energy modes of the jump

  const static bool CILINE2G_isRegistered = MarmotElementFactory::
    registerElement( "CILINE2G", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 2, 4 >, FullIntegration >() );

  const static bool CITRI3G_isRegistered = MarmotElementFactory::
    registerElement( "CITRI3G", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 6 >, FullIntegration >() );

  const static bool CIQUAD4G_isRegistered = MarmotElementFactory::
    registerElement( "CIQUAD4G", makeCohesiveFactoryFunction< CohesiveInterfaceElement< 3, 8 >, FullIntegration >() );

} // namespace Marmot::Elements::Registration
