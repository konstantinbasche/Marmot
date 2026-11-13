#include "Marmot/CohesiveFrictionalCrack.h"
#include "Marmot/MarmotMaterialCohesiveFactory.h"

namespace {

  using namespace Marmot::Materials;

  const static bool CohesiveFrictionalCrackIsRegistered = Marmot::Registry::MaterialCohesiveFactory::registerMaterial<
    CohesiveFrictionalCrack >( "COHESIVEFRICTIONALCRACK" );

} // namespace
