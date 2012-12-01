#ifndef DUNE_STUFF_GRID_PROVIDER_GMSH_HH
#define DUNE_STUFF_GRID_PROVIDER_GMSH_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID

#include <sstream>
#include <type_traits>

#include <boost/assign/list_of.hpp>

#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/stuff/common/parameter/tree.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Provider {

/**
 * \brief   Gmsh grid provider
 */
#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Gmsh : public Interface<GridImp>
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  typedef Interface<GridType> BaseType;

  typedef Gmsh<GridType> ThisType;

  //! Dimension of the provided grid.
  static const unsigned int dim = BaseType::dim;

  //! Type of the grids coordinates.
  typedef typename BaseType::CoordinateType CoordinateType;

  //! Unique identifier: \c stuff.grid.provider.gmsh
  static const std::string id()
  {
    return BaseType::id() + ".gmsh";
  }

  Gmsh(const std::string filename)
  {
    dune_static_assert(!(Dune::is_same<GridType, Dune::YaspGrid<dim>>::value),
                       "GmshReader does not work with YaspGrid!");
    dune_static_assert(!(Dune::is_same<GridType, Dune::SGrid<2, 2>>::value), "GmshReader does not work with SGrid!");
    grid_ = Dune::shared_ptr<GridType>(GmshReader<GridType>::read(filename));
  }

  Gmsh(ThisType& other)
    : grid_(other.grid_)
  {
  }

  Gmsh(const ThisType& other)
    : grid_(other.grid_)
  {
  }

  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Dune::Stuff::Common::ExtendedParameterTree extendedParamTree;
    if (paramTree.hasSub(subName))
      extendedParamTree = paramTree.sub(subName);
    else
      extendedParamTree = paramTree;
    // create and return
    if (!extendedParamTree.hasKey("filename"))
      DUNE_THROW(Dune::RangeError,
                 "\nMissing key 'filename' in the following Dune::ParameterTree:\n"
                     << extendedParamTree.reportString("  "));
    const std::string filename = extendedParamTree.get("filename", "meaningless_default_value");
    return Gmsh(filename);
  }

  ThisType& operator=(ThisType& other)
  {
    if (this != &other) {
      grid_ = other.grid();
    }
    return this;
  }

  //! access to shared ptr
  virtual Dune::shared_ptr<GridType> grid()
  {
    return grid_;
  }

  //! const access to shared ptr
  virtual const Dune::shared_ptr<const GridType> grid() const
  {
    return grid_;
  }

private:
  Dune::shared_ptr<GridType> grid_;
}; // class Gmsh

} // namespace Provider
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_GMSH_HH
