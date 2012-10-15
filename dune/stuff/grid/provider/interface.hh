
#ifndef DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
#define DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH

// dune-common
#include <dune/common/fvector.hh>
#include <dune/common/shared_ptr.hh>

namespace Dune {

namespace Stuff {

namespace Grid {

namespace Provider {

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridImp>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
class Interface
{
public:
  typedef GridImp GridType;

  typedef Dune::FieldVector<typename GridType::ctype, GridType::dimension> CoordinateType;

  virtual std::string id() const = 0;

  virtual GridType& grid() = 0;

  virtual const GridType& grid() const = 0;

  virtual Dune::shared_ptr<GridType> gridPtr() = 0;

  virtual const Dune::shared_ptr<const GridType> gridPtr() const = 0;

  virtual void visualize(const std::string) const = 0;

}; // class Interface

} // namespace Provider

} // namespace Grid

} // namespace Stuff

} // namespace Dune

#endif // DUNE_STUFF_GRID_PROVIDER_INTERFACE_HH
