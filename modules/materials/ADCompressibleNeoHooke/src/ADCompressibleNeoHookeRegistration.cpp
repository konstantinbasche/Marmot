#include "Marmot/ADCompressibleNeoHooke.h"
#include "Marmot/MarmotMaterialFiniteStrainFactory.h"

namespace Marmot::Materials {

  namespace Registration {

    using namespace MarmotLibrary;

    const static bool
      ADCompressibleNeoHookeRegistered = MarmotMaterialFiniteStrainFactory::registerMaterial< ADCompressibleNeoHooke >(
        "ADCOMPRESSIBLENEOHOOKE" );

  } // namespace Registration
} // namespace Marmot::Materials
