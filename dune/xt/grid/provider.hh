// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/configuration.hh>

#include "provider/interface.hh"
#include "provider/cube.hh"

namespace Dune {
namespace XT {


template <class GridType = Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>>
class GridProviders
{
public:
  typedef XT::Grid::ProviderInterface<GridType> InterfaceType;

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
    namespace Providers = XT::Grid::Providers;
    return {Providers::Cube<GridType>::static_id()};
  } // ... available()

  static Common::Configuration default_config(const std::string type, const std::string subname = "")
  {
    namespace Providers = XT::Grid::Providers;
    if (type == Providers::Cube<GridType>::static_id())
      return Providers::Cube<GridType>::default_config(subname);
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const Common::Configuration& config)
  {
    return create(config.get<std::string>("type"), config);
  }

  static std::unique_ptr<InterfaceType> create(const std::string& type = available()[0],
                                               const Common::Configuration config = Common::Configuration())
  {
    namespace Providers = XT::Grid::Providers;
    if (type == Providers::Cube<GridType>::static_id())
      return call_create<Providers::Cube<GridType>>(config);
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
}; // class GridProviders


} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_HH
