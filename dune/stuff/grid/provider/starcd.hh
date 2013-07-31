#ifndef DUNE_STUFF_GRID_PROVIDER_STARCD_HH
#define DUNE_STUFF_GRID_PROVIDER_STARCD_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID
//#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
//#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined
// UGGRID

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <type_traits>

#include <boost/assign/list_of.hpp>

#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/sgrid.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/logging.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {


/**
 * \brief   StarCD grid provider
 *
 *          Implemented for dimensions 1, 2, and 3.
 *
 */
#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class GridProviderStarCD : public GridProviderInterface<GridImp>
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  typedef GridProviderInterface<GridType> BaseType;

  typedef GridProviderStarCD<GridType> ThisType;

  //! Dimension of the provided grid.
  static const unsigned int dim = BaseType::dim;

  //! Type of the grids coordinates.
  typedef typename BaseType::CoordinateType CoordinateType;

  //! Unique identifier: \c stuff.grid.provider.starcd
  static const std::string id()
  {
    return BaseType::id() + ".starcd";
  }

  GridProviderStarCD(const std::string filename, std::ostream& out = Dune::Stuff::Common::Logger().devnull())
  {
    // set up the grid factory
    GridFactory<GridType> factory;

    // set the name of the vertex file
    const std::string vertexFileName = filename + ".vrt";

    // set the vertex input stream
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

    // read the vertices
    int numberOfVertices = 0;
    Dune::FieldVector<double, dim> position;

    out << "Reading " << vertexFileName << " ...   " << std::flush;
    while (std::getline(vertexFile, line)) {
      numberOfVertices++;
      const std::vector<double> items = Dune::Stuff::Common::tokenize<double>(line, " ");
      if (items.size() != dim + 1)
        DUNE_THROW(Dune::IOError, "Error: " << items.size() << " = items.size() != dim + 1 = " << dim + 1 << "!");
      for (unsigned int ii = 0; ii < dim; ++ii)
        position[ii] = items[ii + 1];
      factory.insertVertex(position);
    }
    out << "done: " << numberOfVertices << " vertices read." << std::endl;

    // set the name of the element file
    std::string elementFileName = filename + ".cel";

    // set the element input stream
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

    // read the elements
    unsigned int numberOfElements      = 0;
    unsigned int numberOfPrisms        = 0;
    unsigned int numberOfCubes         = 0;
    unsigned int numberOfVerticesCube  = (unsigned int)pow(2, dim);
    unsigned int numberOfVerticesPrism = 6;

    std::vector<unsigned int> cubeVertices(numberOfVerticesCube);
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

      const std::vector<int> items1 = Dune::Stuff::Common::tokenize<int>(firstLine, " ");
      std::vector<int> items2       = Dune::Stuff::Common::tokenize<int>(secondLine, " ");

      // Erase zeros at the beginning of the line
      items2.erase(std::remove(items2.begin(), items2.end(), 0), items2.end());

      if ((!items1[0] == items2[0]) || (!numberOfElements == items2[0]))
        DUNE_THROW(Dune::IOError, "Elementindices do not correspond!");

      if (items2.size() == numberOfVerticesCube + 1) // cube
      {
        numberOfCubes++;
        for (unsigned int k = 0; k < numberOfVerticesCube; k++) {
          cubeVertices[k] = items2[k + 1] - 1;
        }
        if (dim > 1)
          std::swap(cubeVertices[2], cubeVertices[3]);
        if (dim == 3)
          std::swap(cubeVertices[6], cubeVertices[7]);
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, dim), cubeVertices);
      } else if ((items2.size() == numberOfVerticesPrism + 1) && (dim == 3)) // prism
      {
        numberOfPrisms++;
        for (unsigned int k = 0; k < numberOfVerticesPrism; k++) {
          prismVertices[k] = items2[k + 1] - 1;
        }
        factory.insertElement(Dune::GeometryType(Dune::GeometryType::prism, dim), prismVertices);
      } else // not cube or prism
      {
        DUNE_THROW(Dune::IOError, "Type of element " << numberOfElements << " is not cube or prism!");
      }
    } // while loop

    if (!numberOfElements == (numberOfCubes + numberOfPrisms))
      DUNE_THROW(Dune::IOError,
                 "Number of Elements (" << numberOfElements << ") is not equal to number of cubes (" << numberOfCubes
                                        << ") and number of prisms ("
                                        << numberOfPrisms
                                        << ").");

    out << "done: " << numberOfElements << " elements read (" << numberOfPrisms << " prisms and " << numberOfCubes
        << " cubes)." << std::endl;

    // finish off the construction of the grid object
    out << "Starting createGrid() ... " << std::endl;

    grid_ = std::shared_ptr<GridType>(factory.createGrid());
  } // constructor


  GridProviderStarCD(ThisType& other)
    : grid_(other.grid_)
  {
  }

  GridProviderStarCD(const ThisType& other)
    : grid_(other.grid_)
  {
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["filename"] = "path_to_starcd_filename_prefix";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... createDefaultSettings(...)

  static ThisType* create(const Dune::ParameterTree& _settings, const std::string subName = id())
  {
    // get correct _settings
    Dune::Stuff::Common::ExtendedParameterTree settings;
    if (_settings.hasSub(subName))
      settings = _settings.sub(subName);
    else
      settings = _settings;
    // create and return
    if (!settings.hasKey("filename"))
      DUNE_THROW(Dune::RangeError,
                 "\nMissing key 'filename' in the following Dune::ParameterTree:\n" << settings.reportString("  "));
    const std::string filename = settings.get("filename", "meaningless_default_value");
    return new ThisType(filename);
  } // ... create(...)

  ThisType& operator=(ThisType& other)
  {
    if (this != &other) {
      grid_ = other.grid();
    }
    return this;
  }

  //! access to shared ptr
  virtual std::shared_ptr<GridType> grid()
  {
    return grid_;
  }

  //! const access to shared ptr
  virtual const std::shared_ptr<const GridType> grid() const
  {
    return grid_;
  }

private:
  std::shared_ptr<GridType> grid_;
}; // class GridProviderStarCD

} // namespace Stuff
} // namespace Dune

//#endif // defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined
// UGGRID
//#endif // HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_STARCD_HH
