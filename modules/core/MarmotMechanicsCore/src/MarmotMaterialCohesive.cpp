#include "Marmot/MarmotMaterialCohesive.h"
#include <algorithm>
#include <vector>

double Marmot::Materials::MaterialCohesive::getMaximumStiffness( const Response&        response,
                                                                 const Eigen::Vector3d& jump ) const
{
  const int nStateVars = getNumberOfRequiredStateVars();

  std::vector< double > stateVarsCopy( nStateVars, 0.0 );
  if ( response.stateVars != nullptr && nStateVars > 0 ) {
    std::copy_n( response.stateVars, nStateVars, stateVarsCopy.begin() );
  }

  Response responseCopy{ response.traction, stateVarsCopy.data(), response.elasticEnergyDensity, response.dissipation };
  Tangents tangents;
  const Increment increment{ jump, Eigen::Vector3d::Zero(), 0.0, 1.0 };

  computeTraction( responseCopy, tangents, increment );

  return std::max( 0.0, tangents.dTraction_dJump.diagonal().maxCoeff() );
}
