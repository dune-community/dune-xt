#ifndef DUNE_STUFF_GRID_PROVIDER_HH
#define DUNE_STUFF_GRID_PROVIDER_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID

#include <functional>
#include <dune/common/shared_ptr.hh>
#include <dune/stuff/aliases.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"
#include "provider/gmsh.hh"
#include "provider/starcd.hh"

namespace Dune {
namespace Stuff {
namespace Grid {
namespace Provider {


std::vector<std::string> types()
{
  std::vector<std::string> ret;
  ret.push_back("stuff.grid.provider.cube");
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
  ret.push_back("stuff.grid.provider.gmsh");
#endif // defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined
// UGGRID
#endif // HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
  ret.push_back("stuff.grid.provider.starcd");
  return ret;
} // std::vector< const std::string > types()


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
Dune::ParameterTree createSampleDescription(const std::string type)
{
  if (type == "stuff.grid.provider.cube") {
    return Dune::Stuff::Grid::Provider::GenericCube<GridType>::createSampleDescription();
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
  } else if (type == "stuff.grid.provider.gmsh") {
    typedef Dune::Stuff::Grid::Provider::Gmsh<GridType> ProviderType;
    return ProviderType::createSampleDescription();
#endif // defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined
// UGGRID
#endif // HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
  } else if (type == "stuff.grid.provider.starcd") {
    typedef Dune::Stuff::Grid::Provider::StarCD<GridType> ProviderType;
    return ProviderType::createSampleDescription();
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown grid provider '" << type
                    << "' requested!");
} // ... create(...)

template <class GridType>
struct FunctionObject
{
  typedef std::function<Interface<GridType>*(const std::string)> Type;
};

template <class ProviderType>
std::pair<std::string, typename FunctionObject<typename ProviderType::GridType>::Type>
make(const Dune::ParameterTree paramTree)
{
  return std::make_pair(ProviderType::id(),
                        std::bind(&ProviderType::createFromDescription, paramTree, std::placeholders::_1));
}

#define DSGP_MAKE(type) make<type<GridType>>(paramTree)

#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
Interface<GridType>* create(const std::string& type = "stuff.grid.provider.cube",
                            const Dune::ParameterTree paramTree = Dune::ParameterTree())
{
  using namespace DSG::Provider;
  typedef std::map<std::string, std::pair<std::string, typename FunctionObject<GridType>::Type>> MapType;

  MapType ptr_map = {{"stuff.grid.provider.cube", DSGP_MAKE(GenericCube)},
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
                     {"stuff.grid.provider.gmsh", DSGP_MAKE(Gmsh)},
#endif
                     {"stuff.grid.provider.starcd", DSGP_MAKE(StarCD)}};
  auto id_func_it = ptr_map.find(type);
  if (id_func_it == ptr_map.end())
    DUNE_THROW(Dune::RangeError, "\nERROR: unknown grid provider '" << type << "' requested!");
  auto id_func = id_func_it->second;
  return id_func.second(id_func.first);
} // Interface< GridImp >* create(const std::string& type, const Dune::ParameterTree paramTree = Dune::ParameterTree())

} // namespace Provider
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_HH
