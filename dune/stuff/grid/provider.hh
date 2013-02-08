#ifndef DUNE_STUFF_GRID_PROVIDER_HH
#define DUNE_STUFF_GRID_PROVIDER_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID

#include <dune/common/shared_ptr.hh>

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
    typedef Dune::Stuff::Grid::Provider::Cube<GridType> ProviderType;
    return ProviderType::createSampleDescription();
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


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::SGrid<2, 2>>
#endif // defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
Interface<GridType>* create(const std::string& type = "stuff.grid.provider.cube",
                            const Dune::ParameterTree paramTree = Dune::ParameterTree())
{
  // choose provider
  if (type == "stuff.grid.provider.cube") {
    typedef Dune::Stuff::Grid::Provider::Cube<GridType> CubeProviderType;
    CubeProviderType* cubeProvider = new CubeProviderType(CubeProviderType::createFromParamTree(paramTree));
    return cubeProvider;
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
  } else if (type == "stuff.grid.provider.gmsh") {
    typedef Dune::Stuff::Grid::Provider::Gmsh<GridType> GmshProviderType;
    GmshProviderType* gmshProvider = new GmshProviderType(GmshProviderType::createFromParamTree(paramTree));
    return gmshProvider;
#endif // defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined
// UGGRID
#endif // HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
  } else if (type == "stuff.grid.provider.starcd") {
    typedef Dune::Stuff::Grid::Provider::StarCD<GridType> StarCDProviderType;
    StarCDProviderType* starCDProvider = new StarCDProviderType(StarCDProviderType::createFromParamTree(paramTree));
    return starCDProvider;
  } else
    DUNE_THROW(Dune::RangeError, "\nERROR: unknown grid provider '" << type << "' requested!");
} // Interface< GridImp >* create(const std::string& type, const Dune::ParameterTree paramTree = Dune::ParameterTree())

} // namespace Provider
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_HH
