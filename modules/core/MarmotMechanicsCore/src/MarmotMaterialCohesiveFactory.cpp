#include "Marmot/MarmotMaterialCohesiveFactory.h"
#include "Marmot/MarmotJournal.h"

using namespace MarmotLibrary;

MarmotMaterialCohesive* MarmotMaterialCohesiveFactory::createMaterial( const std::string& materialName,
                                                                       const double*      materialProperties,
                                                                       int                nMaterialProperties,
                                                                       int                materialNumber )
{
  auto& map = materialFactoryFunctionByName();
  auto  it  = map.find( materialName );
  if ( it == map.end() ) {
    std::string reg = "Registered cohesive materials are: ";
    for ( const auto& pair : map ) {
      reg += pair.first + ", ";
    }
    throw std::invalid_argument( MakeString() << __PRETTY_FUNCTION__
                                              << " Cohesive material " + materialName + " not registered! " + reg );
  }

  return it->second( materialProperties, nMaterialProperties, materialNumber );
}

MarmotMaterialCohesiveFactory::MaterialFactoryMap& MarmotMaterialCohesiveFactory::materialFactoryFunctionByName()
{
  static MaterialFactoryMap map;
  return map;
}
