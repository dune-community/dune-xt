// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2015)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_STUFF_GRID_PROVIDER_HH
#define DUNE_STUFF_GRID_PROVIDER_HH

#include <memory>

#include <type_traits>

#include <dune/stuff/common/configuration.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"

namespace Dune {
namespace Stuff {

#if HAVE_DUNE_GRID

template <class GridType = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>
class GridProviders
{
public:
  typedef Stuff::Grid::ProviderInterface<GridType> InterfaceType;

protected:
  template <class GridProviderType>
  static std::unique_ptr<InterfaceType> call_create(const Common::Configuration& config)
  {
    if (config.empty())
      return GridProviderType::create();
    else
      return GridProviderType::create(config);
  } // ... call_create(...)

public:
  static std::vector<std::string> available()
  {
    namespace Providers = Stuff::Grid::Providers;
    return {Providers::Cube<GridType>::static_id()};
  } // ... available()

  static Common::Configuration default_config(const std::string type, const std::string subname = "")
  {
    namespace Providers = Stuff::Grid::Providers;
    if (type == Providers::Cube<GridType>::static_id())
      return Providers::Cube<GridType>::default_config(subname);
    else
      DUNE_THROW(Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const Common::Configuration& config)
  {
    return create(config.get<std::string>("type"), config);
  }

  static std::unique_ptr<InterfaceType> create(const std::string& type = available()[0],
                                               const Common::Configuration config = Common::Configuration())
  {
    namespace Providers = Stuff::Grid::Providers;
    if (type == Providers::Cube<GridType>::static_id())
      return call_create<Providers::Cube<GridType>>(config);
    else
      DUNE_THROW(Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
}; // class GridProviders

#else // HAVE_DUNE_GRID

template <class GridType>
class GridProviders
{
  static_assert(AlwaysFalse<GridType>::value, "You are missing dune-grid!");
};

#endif // HAVE_DUNE_GRID

} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_PROVIDER_HH
