#ifndef DUNE_STUFF_GRID_PROVIDER_GMSH_HH
#define DUNE_STUFF_GRID_PROVIDER_GMSH_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID

#include <memory>
#include <sstream>
#include <type_traits>

#include <boost/assign/list_of.hpp>

#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/stuff/common/parameter/tree.hh>

#include "interface.hh"

namespace Dune {
namespace Stuff {


/**
 * \brief   Gmsh grid provider
 */
#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::GridSelector::GridType>
#else
template <class GridImp = Dune::SGrid<2, 2>>
#endif
class GridProviderGmsh : public GridProviderInterface<GridImp>
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  typedef GridProviderInterface<GridType> BaseType;

  typedef GridProviderGmsh<GridType> ThisType;

  //! Dimension of the provided grid.
  static const unsigned int dim = BaseType::dim;

  //! Type of the grids coordinates.
  typedef typename BaseType::CoordinateType CoordinateType;

  //! Unique identifier: \c stuff.grid.provider.gmsh
  static const std::string id()
  {
    return BaseType::id() + ".gmsh";
  }

  GridProviderGmsh(const std::string filename)
  {
    dune_static_assert(!(Dune::is_same<GridType, Dune::YaspGrid<dim>>::value),
                       "GmshReader does not work with YaspGrid!");
    dune_static_assert(!(Dune::is_same<GridType, Dune::SGrid<2, 2>>::value), "GmshReader does not work with SGrid!");
    grid_ = std::shared_ptr<GridType>(GmshReader<GridType>::read(filename));
  }

  GridProviderGmsh(ThisType& other)
    : grid_(other.grid_)
  {
  }

  GridProviderGmsh(const ThisType& other)
    : grid_(other.grid_)
  {
  }

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["filename"] = "path_to_g.msh";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... createSampleDescription(...)

  static ThisType* create(const Dune::ParameterTree& paramTree, const std::string subName = id())
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
    return new ThisType(filename);
  }

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
}; // class GridProviderGmsh


} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID
#endif // defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined
// UGGRID
#endif // HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG

#endif // DUNE_STUFF_GRID_PROVIDER_GMSH_HH
