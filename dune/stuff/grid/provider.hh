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
#include <dune/common/parametertree.hh>

#include <dune/stuff/common/color.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"
#include "provider/gmsh.hh"
#include "provider/starcd.hh"

namespace Dune {
namespace Stuff {


std::vector<std::string> availableGridProviders()
{
  return {"gridprovider.cube"
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
          ,
          "gridprovider.gmsh"
#endif
          ,
          "gridprovider.starcd"};
} // ... availableGridProviders()


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else
template <class GridType = Dune::SGrid<2, 2>>
#endif
Dune::ParameterTree createSampleGridProviderDescription(const std::string type)
{
  if (type == "gridprovider.cube") {
    return GridProviderCube<GridType>::createSampleDescription();
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
  } else if (type == "gridprovider.gmsh") {
    return GridProviderGmsh<GridType>::createSampleDescription();
#endif
  } else if (type == "gridprovider.starcd") {
    return GridProviderStarCD<GridType>::createSampleDescription();
  } else
    DUNE_THROW(Dune::RangeError,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown grid provider '" << type
                    << "' requested!");
} // ... createSampleGridProviderDescription(...)


template <class GridType>
struct DSGP_FunctionObject
{
  typedef std::function<GridProviderInterface<GridType>*(const std::string)> Type;
};

template <class ProviderType>
std::pair<std::string, typename DSGP_FunctionObject<typename ProviderType::GridType>::Type>
make(const Dune::ParameterTree paramTree)
{
  return std::make_pair(ProviderType::id(), std::bind(&ProviderType::create, paramTree, std::placeholders::_1));
}

#define DSGP_MAKE(type) make<type<GridType>>(description)


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else
template <class GridType = Dune::SGrid<2, 2>>
#endif
GridProviderInterface<GridType>* createGridProvider(const std::string& type = "grid.provider.cube",
                                                    const Dune::ParameterTree description = Dune::ParameterTree())
{
  typedef std::map<std::string, std::pair<std::string, typename DSGP_FunctionObject<GridType>::Type>> MapType;

  MapType ptr_map = {{"gridprovider.cube", DSGP_MAKE(GridProviderCube)}
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
                     ,
                     {"gridprovider.gmsh", DSGP_MAKE(GridProviderGmsh)}
#endif
                     ,
                     {"gridprovider.starcd", DSGP_MAKE(GridProviderStarCD)}};
  auto id_func_it = ptr_map.find(type);
  if (id_func_it == ptr_map.end())
    DUNE_THROW(Dune::RangeError,
               "\n" << Common::colorStringRed("ERROR:") << " unknown grid provider '" << type << "' requested!");
  auto id_func = id_func_it->second;
  return id_func.second(id_func.first);
} // ... createGridProvider(...)


} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_HH
