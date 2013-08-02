#ifndef DUNE_STUFF_GRID_PROVIDER_HH
#define DUNE_STUFF_GRID_PROVIDER_HH

#if HAVE_DUNE_GRID

#include <dune/common/parametertree.hh>

#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/color.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"
#include "provider/gmsh.hh"
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
    return {"gridprovider.cube", "gridprovider.gmsh", "gridprovider.starcd"};
  } // ... available()

  static Dune::ParameterTree defaultSettings(const std::string type, const std::string subname = "")
  {
    if (type == "gridprovider.cube") {
      return GridProviderCube<GridType>::defaultSettings(subname);
    } else if (type == "gridprovider.gmsh") {
      return GridProviderGmsh<GridType>::defaultSettings(subname);
    } else if (type == "gridprovider.starcd") {
      return GridProviderStarCD<GridType>::defaultSettings(subname);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... defaultSettings(...)

  static GridProviderInterface<GridType>* create(const std::string& type = available()[0],
                                                 const Dune::ParameterTree settings = Dune::ParameterTree())
  {
    if (type == "gridprovider.cube")
      return GridProviderCube<GridType>::create(settings);
    else if (type == "gridprovider.gmsh")
      return GridProviderGmsh<GridType>::create(settings);
    else if (type == "gridprovider.starcd")
      return GridProviderStarCD<GridType>::create(settings);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... create(...)
}; // class GridProviders


template <int dim>
class GridProviders<Dune::SGrid<dim, dim>>
{
public:
  static std::vector<std::string> available()
  {
    return {"gridprovider.cube", "gridprovider.starcd"};
  } // ... available()

  static Dune::ParameterTree defaultSettings(const std::string type, const std::string subname = "")
  {
    if (type == "gridprovider.cube") {
      return GridProviderCube<Dune::SGrid<dim, dim>>::defaultSettings(subname);
    } else if (type == "gridprovider.starcd") {
      return GridProviderStarCD<Dune::SGrid<dim, dim>>::defaultSettings(subname);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... createDefaultSettings(...)

  static GridProviderInterface<Dune::SGrid<dim, dim>>*
  create(const std::string& type = available()[0], const Dune::ParameterTree settings = Dune::ParameterTree())
  {
    if (type == "gridprovider.cube")
      return GridProviderCube<Dune::SGrid<dim, dim>>::create(settings);
    else if (type == "gridprovider.starcd")
      return GridProviderStarCD<Dune::SGrid<dim, dim>>::create(settings);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... create(...)
}; // class GridProviders< ... SGrid ... >


template <int dim>
class GridProviders<Dune::YaspGrid<dim>>
{
public:
  static std::vector<std::string> available()
  {
    return {"gridprovider.cube", "gridprovider.starcd"};
  } // ... available()

  static Dune::ParameterTree defaultSettings(const std::string type, const std::string subname = "")
  {
    if (type == "gridprovider.cube") {
      return GridProviderCube<Dune::YaspGrid<dim>>::defaultSettings(subname);
    } else if (type == "gridprovider.starcd") {
      return GridProviderStarCD<Dune::YaspGrid<dim>>::defaultSettings(subname);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... createDefaultSettings(...)

  static GridProviderInterface<Dune::YaspGrid<dim>>* create(const std::string& type = available()[0],
                                                            const Dune::ParameterTree settings = Dune::ParameterTree())
  {
    if (type == "gridprovider.cube")
      return GridProviderCube<Dune::YaspGrid<dim>>::create(settings);
    else if (type == "gridprovider.starcd")
      return GridProviderStarCD<Dune::YaspGrid<dim>>::create(settings);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Common::colorStringRed("ERROR:") << " unknown gridprovider '" << type << "' requested!");
  } // ... create(...)
}; // class GridProviders< ... SGrid ... >


} // namespace Stuff
} // namespace Dune

#endif // HAVE_DUNE_GRID

#endif // DUNE_STUFF_GRID_PROVIDER_HH
