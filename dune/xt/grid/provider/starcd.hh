// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2015)
//   Kirsten Weber   (2012 - 2013)
//   Rene Milk       (2012 - 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_PROVIDER_STARCD_HH
#define DUNE_XT_GRID_PROVIDER_STARCD_HH

#if HAVE_DUNE_GRID

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <type_traits>

#include <boost/assign/list_of.hpp>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/logging.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/grid/provider/interface.hh>

namespace Dune {
namespace XT {

/**
 * \brief   StarCD grid provider
 *
 *          Implemented for dimensions 1, 2, and 3.
 *
 */
template <class GridImp>
class GridProviderStarCD : public Grid::ProviderInterface<GridImp>
{
  typedef Grid::ProviderInterface<GridImp> BaseType;
  typedef GridProviderStarCD<GridImp> ThisType;

public:
  using typename BaseType::GridType;
  using BaseType::dimDomain;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".starcd";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config("filename_prefix", "sample");
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    return Common::make_unique<ThisType>(cfg.get("filename", default_cfg.get<std::string>("filename")));
  }

  GridProviderStarCD(const std::string& filename)
  {
    std::ostream& out = Dune::XT::Common::Logger().devnull();

    // set up the grid factory
    GridFactory<GridType> factory;

    // read the vertices
    const std::string vertexFileName = filename + ".vrt";
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
    std::string elementFileName = filename + ".cel";
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
    size_t numberOfElements      = 0;
    size_t numberOfPrisms        = 0;
    size_t numberOfCubes         = 0;
    size_t numberOfVerticesCube  = pow(2, dimDomain);
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
      auto items2       = Dune::XT::Common::tokenize<int>(secondLine, " ");

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
                                        << ") and number of prisms ("
                                        << numberOfPrisms
                                        << ").");

    out << "done: " << numberOfElements << " elements read (" << numberOfPrisms << " prisms and " << numberOfCubes
        << " cubes)." << std::endl;

    // finish the construction of the grid object
    out << "Starting createGrid() ... " << std::endl;

    grid_ = std::shared_ptr<GridType>(factory.createGrid());
  } // GridProviderStarCD(...)

  GridProviderStarCD(ThisType&& source) = default;
  GridProviderStarCD(const ThisType& other) = default;

  virtual ~GridProviderStarCD() = default;

  ThisType& operator=(ThisType&& source) = default;
  ThisType& operator=(const ThisType& other) = default;

  virtual const GridType& grid() const override final
  {
    return *grid_;
  }

  virtual GridType& grid() override final
  {
    return *grid_;
  }

  const std::shared_ptr<const GridType> grid_ptr() const
  {
    return grid_;
  }

  std::shared_ptr<GridType> grid_ptr()
  {
    return grid_;
  }

  virtual std::unique_ptr<Grid::ConstProviderInterface<GridType>> copy() const override final
  {
    DUNE_THROW(NotImplemented, "");
    return nullptr;
  }

  virtual std::unique_ptr<Grid::ProviderInterface<GridType>> copy() override final
  {
    DUNE_THROW(NotImplemented, "");
    return nullptr;
  }

private:
  std::shared_ptr<GridType> grid_;
}; // class GridProviderStarCD

} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_XT_GRID_PROVIDER_STARCD_HH
