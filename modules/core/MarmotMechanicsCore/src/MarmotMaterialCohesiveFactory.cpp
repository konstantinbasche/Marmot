#include "Marmot/MarmotMaterialCohesiveFactory.h"
#include "Marmot/MarmotJournal.h"

namespace Marmot::Registry {

  Marmot::Materials::MaterialCohesive* MaterialCohesiveFactory::createMaterial( const std::string& materialName,
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
      throw std::invalid_argument( Marmot::Utilities::MakeString()
                                   << __PRETTY_FUNCTION__
                                   << " Cohesive material " + materialName + " not registered! " + reg );
    }

    return it->second( materialProperties, nMaterialProperties, materialNumber );
  }

  MaterialCohesiveFactory::MaterialFactoryMap& MaterialCohesiveFactory::materialFactoryFunctionByName()
  {
    static MaterialFactoryMap map;
    return map;
  }

} // namespace Marmot::Registry
