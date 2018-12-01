// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Kirsten Weber   (2012 - 2013)
//   René Fritze     (2012 - 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_STARCD_HH
#define DUNE_XT_GRID_GRIDPROVIDER_STARCD_HH

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <type_traits>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/grid/gridprovider.hh>

namespace Dune {
namespace XT {
namespace Grid {


static inline std::string starcd_gridprovider_id()
{
  return "xt.grid.gridprovider.starcd";
}


static inline Common::Configuration starcd_gridprovider_default_config()
{
  Common::Configuration config;
  config["type"] = starcd_gridprovider_id();
  config["filename_prefix"] = "grid";
  return config;
}


/**
 * \brief   StarCD grid provider
 *
 *          Implemented for dimensions 1, 2, and 3.
 *
 */
template <class GridType>
class StarCDGridProviderFactory
{
  static const size_t dimDomain = GridType::dimension;

public:
  static const bool available = true;

  static const std::string static_id()
  {
    return starcd_gridprovider_id();
  }

  static Common::Configuration default_config()
  {
    return starcd_gridprovider_default_config();
  }

  static GridProvider<GridType> create(const std::string& filename_prefix)
  {
    using std::pow;
    std::ostream& out = Dune::XT::Common::Logger().devnull();

    // set up the grid factory
    GridFactory<GridType> factory;

    // read the vertices
    const std::string vertexFileName = filename_prefix + ".vrt";
    std::ifstream vertexFile(vertexFileName);
    if (!vertexFile)
      DUNE_THROW(Dune::IOError, "Could not open " << vertexFileName);
    std::string line;
    if (!std::getline(vertexFile, line))
      DUNE_THROW(Dune::IOError, "File " << vertexFileName << " is too short!");
    if (!(line == "PROSTAR_VERTEX"))
      DUNE_THROW(Dune::IOError,
                 "First line of File " << vertexFileName << " (" << line << "is not equal to 'PROSTAR_VERTEX' !");
    if (!std::getline(vertexFile, line))
      DUNE_THROW(Dune::IOError, "File " << vertexFileName << " is too short!");
    size_t numberOfVertices = 0;
    Dune::FieldVector<double, dimDomain> position;
    out << "Reading " << vertexFileName << " ...   " << std::flush;
    while (std::getline(vertexFile, line)) {
      numberOfVertices++;
      const std::vector<double> items = Dune::XT::Common::tokenize<double>(line, " ");
      if (items.size() != dimDomain + 1)
        DUNE_THROW(Dune::IOError, "Error: " << items.size() << " = items.size() != dim + 1 = " << dimDomain + 1 << "!");
      for (size_t ii = 0; ii < dimDomain; ++ii)
        position[ii] = items[ii + 1];
      factory.insertVertex(position);
    }
    out << "done: " << numberOfVertices << " vertices read." << std::endl;

    // read the elements
    std::string elementFileName = filename_prefix + ".cel";
    std::ifstream elementFile(elementFileName);
    if (!elementFile)
      DUNE_THROW(Dune::IOError, "Could not open " << elementFileName);
    if (!std::getline(elementFile, line))
      DUNE_THROW(Dune::IOError, "File " << elementFileName << " is too short!");
    if (!(line == "PROSTAR_CELL"))
      DUNE_THROW(Dune::IOError,
                 "First line of File " << elementFileName << " (" << line << "is not equal to 'PROSTAR_CELL' !");
    if (!std::getline(elementFile, line))
      DUNE_THROW(Dune::IOError, "File " << elementFileName << " is too short!");
    size_t numberOfElements = 0;
    size_t numberOfPrisms = 0;
    size_t numberOfCubes = 0;
    size_t numberOfVerticesCube = pow(2, dimDomain);
    size_t numberOfVerticesPrism = 6;
    std::vector<unsigned int> cubeVertices(numberOfVerticesCube); // unsigned int required by the grid factory
    std::vector<unsigned int> prismVertices(numberOfVerticesPrism);
    std::string firstLine;
    std::string secondLine;
    out << "Reading " << elementFileName << " ...   " << std::flush;
    while (std::getline(elementFile, firstLine)) {
      if (!std::getline(elementFile, secondLine))
        DUNE_THROW(Dune::IOError,
                   "No vertex data available in file " << elementFileName << " for element " << numberOfElements + 1
                                                       << "!");
      numberOfElements++;

      const auto items1 = Dune::XT::Common::tokenize<int>(firstLine, " ");
      auto items2 = Dune::XT::Common::tokenize<int>(secondLine, " ");

      // Erase zeros at the beginning of the line
      items2.erase(std::remove(items2.begin(), items2.end(), 0), items2.end());

      if (!(items1[0] == items2[0]) || !(numberOfElements == items2[0]))
        DUNE_THROW(Dune::IOError, "Elementindices do not correspond!");

      if (items2.size() == numberOfVerticesCube + 1) { // cube
        numberOfCubes++;
        for (size_t k = 0; k < numberOfVerticesCube; k++)
          cubeVertices[k] = items2[k + 1] - 1;
        if (dimDomain > 1)
          std::swap(cubeVertices[2], cubeVertices[3]);
        if (dimDomain == 3)
          std::swap(cubeVertices[6], cubeVertices[7]);
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, dimDomain), cubeVertices);
      } else if ((items2.size() == numberOfVerticesPrism + 1) && (dimDomain == 3)) { // prism
        numberOfPrisms++;
        for (size_t k = 0; k < numberOfVerticesPrism; k++)
          prismVertices[k] = items2[k + 1] - 1;
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism, dimDomain), prismVertices);
      } else // neither cube or prism
        DUNE_THROW(Dune::IOError, "Type of element " << numberOfElements << " is not cube or prism!");
    } // while loop

    if (!(numberOfElements == (numberOfCubes + numberOfPrisms)))
      DUNE_THROW(Dune::IOError,
                 "Number of Elements (" << numberOfElements << ") is not equal to number of cubes (" << numberOfCubes
                                        << ") and number of prisms (" << numberOfPrisms << ").");

    out << "done: " << numberOfElements << " elements read (" << numberOfPrisms << " prisms and " << numberOfCubes
        << " cubes)." << std::endl;

    // finish the construction of the grid object
    out << "Starting createGrid() ... " << std::endl;

    return GridProvider<GridType>(factory.createGrid());
  } // ... create(...)

  static GridProvider<GridType> create(const Common::Configuration cfg = default_config())
  {
    return create(cfg.get("filename_prefix", default_config().template get<std::string>("filename_prefix")));
  }
}; // class StarCDGridProviderFactory


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDPROVIDER_STARCD_HH
