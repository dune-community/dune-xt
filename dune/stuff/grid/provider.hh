#ifndef DUNE_STUFF_GRID_PROVIDER_HH
#define DUNE_STUFF_GRID_PROVIDER_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_GRID

#include <dune/common/parametertree.hh>

#include <dune/grid/sgrid.hh>

#include <dune/stuff/common/color.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
#include "provider/gmsh.hh"
#endif
#include "provider/starcd.hh"

namespace Dune {
namespace Stuff {


#if defined HAVE_CONFIG_H || defined HAVE_CMAKE_CONFIG
template <class GridType = Dune::GridSelector::GridType>
#else
template <class GridType = Dune::SGrid<2, 2>>
#endif
class GridProviders
{
public:
  static std::vector<std::string> available()
  {
    return {"gridprovider.cube"
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
            ,
            "gridprovider.gmsh"
#endif
            ,
            "gridprovider.starcd"};
  } // ... available()

  static Dune::ParameterTree createSampleDescription(const std::string type, const std::string subname = "")
  {
    if (type == "gridprovider.cube") {
      return GridProviderCube<GridType>::createSampleDescription(subname);
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
    } else if (type == "gridprovider.gmsh") {
      return GridProviderGmsh<GridType>::createSampleDescription(subname);
#endif
    } else if (type == "gridprovider.starcd") {
      return GridProviderStarCD<GridType>::createSampleDescription(subname);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... createSampleDescription(...)

  static GridProviderInterface<GridType>* create(const std::string& type = available()[0],
                                                 const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "gridprovider.cube")
      return GridProviderCube<GridType>::create(description);
#ifdef HAVE_UNSTRUCTURED_GRIDFACTORY
    else if (type == "gridprovider.gmsh")
      return GridProviderGmsh<GridType>::create(description);
#endif
    else if (type == "gridprovider.starcd")
      return GridProviderStarCD<GridType>::create(description);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... create(...)
}; // class GridProviders


} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_HH
