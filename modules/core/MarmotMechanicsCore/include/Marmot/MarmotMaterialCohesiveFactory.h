/* ---------------------------------------------------------------------
 *                                       _
 *  _ __ ___   __ _ _ __ _ __ ___   ___ | |_
 * | '_ ` _ \ / _` | '__| '_ ` _ \ / _ \| __|
 * | | | | | | (_| | |  | | | | | | (_) | |_
 * |_| |_| |_|\__,_|_|  |_| |_| |_|\___/ \__|
 *
 * Unit of Strength of Materials and Structural Analysis
 * University of Innsbruck,
 * 2020 - today
 *
 * festigkeitslehre@uibk.ac.at
 *
 * Konstantin Viktor Basche konstantin.basche@uibk.ac.at
 *
 * This file is part of the MAteRialMOdellingToolbox (marmot).
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of marmot.
 * ---------------------------------------------------------------------
 */
#pragma once
#include "Marmot/MarmotMaterialCohesive.h"
#include "Marmot/MarmotPortability.h"
#include <cassert>
#include <functional>
#include <string>
#include <unordered_map>

namespace MarmotLibrary {

  /**
   * @class MarmotMaterialCohesiveFactory
   * @brief Factory for cohesive traction-separation laws (@ref MarmotMaterialCohesive).
   *
   * Cohesive materials register themselves by name, and cohesive interface elements create them from the name of
   * their material section.
   */
  class MARMOT_API MarmotMaterialCohesiveFactory {
  public:
    /// Function signature for material factory functions registered in the map.
    using materialFactoryFunction = std::function<
      MarmotMaterialCohesive*( const double* materialProperties, int nMaterialProperties, int materialNumber ) >;

    MarmotMaterialCohesiveFactory() = delete;

    /**
     * @brief Create a material instance based on its name and properties.
     * @param[in] materialName        Registered name of the material.
     * @param[in] materialProperties  Array of material properties.
     * @param[in] nMaterialProperties Number of properties in the array.
     * @param[in] materialNumber      Unique identifier for the material instance.
     * @return Pointer to the created material instance.
     * @throws std::invalid_argument if no material is registered under @p materialName.
     */
    static MarmotMaterialCohesive* createMaterial( const std::string& materialName,
                                                   const double*      materialProperties,
                                                   int                nMaterialProperties,
                                                   int                materialNumber );

    /**
     * @brief Register a material with its name and an auto-generated factory function.
     * @param[in] materialName Unique name under which the material is registered.
     * @return True if registration was successful.
     */
    template < class T >
    static bool registerMaterial( const std::string& materialName )
    {
      auto& map = materialFactoryFunctionByName();

      assert( map.find( materialName ) == map.end() && "Material already registered!" );

      map[materialName] =
        []( const double* materialProperties, int nMaterialProperties, int materialNumber ) -> MarmotMaterialCohesive* {
        return new T( materialProperties, nMaterialProperties, materialNumber );
      };
      return true;
    }

  private:
    using MaterialFactoryMap = std::unordered_map< std::string, materialFactoryFunction >;
    static MaterialFactoryMap& materialFactoryFunctionByName();
  };
} // namespace MarmotLibrary
