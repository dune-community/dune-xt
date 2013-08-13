#ifndef DUNE_STUFF_GRID_PROVIDER_HH
#define DUNE_STUFF_GRID_PROVIDER_HH

#if HAVE_DUNE_GRID

#include <dune/common/parametertree.hh>

#include <dune/grid/sgrid.hh>

#include <dune/stuff/common/color.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
#include "provider/gmsh.hh"
#endif
#endif
#include "provider/starcd.hh"

namespace Dune {
namespace Stuff {

template <class GridType = Dune::SGrid<2, 2>>
class GridProviders
{
public:
  static std::vector<std::string> available()
  {
    return
    {
      "gridprovider.cube"
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
          ,
          "gridprovider.gmsh"
#endif
#endif
          ,
          "gridprovider.starcd"
    };
  } // ... available()

  static Dune::ParameterTree defaultSettings(const std::string type, const std::string subname = "")
  {
    if (type == "gridprovider.cube") {
      return GridProviderCube<GridType>::defaultSettings(subname);
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
    } else if (type == "gridprovider.gmsh") {
      return GridProviderGmsh<GridType>::defaultSettings(subname);
#endif
#endif
    } else if (type == "gridprovider.starcd") {
      return GridProviderStarCD<GridType>::defaultSettings(subname);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... createDefaultSettings(...)

  static GridProviderInterface<GridType>* create(const std::string& type = available()[0],
                                                 const Dune::ParameterTree settings = Dune::ParameterTree())
  {
    if (type == "gridprovider.cube")
      return GridProviderCube<GridType>::create(settings);
#if HAVE_ALUGRID || HAVE_ALBERTA || HAVE_UG
#if defined ALUGRID_CONFORM || defined ALUGRID_CUBE || defined ALUGRID_SIMPLEX || defined ALBERTAGRID || defined UGGRID
    else if (type == "gridprovider.gmsh")
      return GridProviderGmsh<GridType>::create(settings);
#endif
#endif
    else if (type == "gridprovider.starcd")
      return GridProviderStarCD<GridType>::create(settings);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... create(...)
}; // class GridProviders


} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_HH
