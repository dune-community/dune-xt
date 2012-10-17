#ifndef DUNE_STUFF_GRID_PROVIDER_GMSH_HH
#define DUNE_STUFF_GRID_PROVIDER_GMSH_HH

// system
#include <sstream>
#include <type_traits>
#include <boost/assign/list_of.hpp>

// dune-common
#include <dune/common/parametertree.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

// dune-grid
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#ifdef HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
//#include <dune/grid/common/mcmgmapper.hh>
//#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/gmshreader.hh>

// local
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
#else
template <class GridImp>
#endif
class Gmsh : public Interface<GridImp>
{
public:
  //! Type of the provided grid.
  typedef GridImp GridType;

  //! Dimension of the provided grid.
  static const unsigned int dim = GridType::dimension;

  //! Type of the grids coordinates.
  typedef Dune::FieldVector<typename GridType::ctype, dim> CoordinateType;

  //! Unique identifier: \c stuff.grid.provider.gmsh
  static const std::string static_id;

  Gmsh(const Dune::ParameterTree paramTree)
  {
    const std::string filename = paramTree.get("mshfile", "sample.msh");
    // read gmshfile
    grid_ = Dune::shared_ptr<GridType>(GmshReader<GridType>::read(filename));
  }

  virtual std::string id() const
  {
    return static_id;
  }

  /**
    \brief  Provides access to the created grid.
    \return Reference to the grid.
    **/
  virtual GridType& grid()
  {
    return *grid_;
  }

  /**
   *  \brief  Provides const access to the created grid.
   *  \return Reference to the grid.
   **/
  virtual const GridType& grid() const
  {
    return *grid_;
  }

  //! access to shared ptr
  virtual Dune::shared_ptr<GridType> gridPtr()
  {
    return grid_;
  }

  //! const access to shared ptr
  virtual const Dune::shared_ptr<const GridType> gridPtr() const
  {
    return grid_;
  }

private:
  const std::string filename_;
  Dune::shared_ptr<GridType> grid_;
}; // class Gmsh

template <typename GridImp>
const std::string Gmsh<GridImp>::static_id = "stuff.grid.provider.gmsh";


} // namespace Provider
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PROVIDER_GMSH_HH
