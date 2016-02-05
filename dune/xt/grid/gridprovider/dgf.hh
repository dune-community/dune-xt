// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2016)
//   Kirsten Weber   (2012)
//   Rene Milk       (2012 - 2013, 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_PROVIDER_DGF_HH
#define DUNE_XT_GRID_PROVIDER_DGF_HH

#include <memory>

#if HAVE_ALUGRID
#include <dune/grid/io/file/dgfparser/dgfalu.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfgridfactory.hh> // How convenient that GridPtr requires DGFGridFactory but ...
#include <dune/grid/io/file/dgfparser/dgfoned.hh>
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#include <dune/grid/io/file/dgfparser/gridptr.hh> // ... does not include it!

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/dgf.hh>
#endif

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/string.hh>

#include "provider.hh"

namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


template <class GridType>
static GridProvider<GridType> create_dgf_grid(const std::string& filename)
{
  return GridProvider<GridType>(GridPtr<GridType>(filename).release());
}


} // namespace internal


static inline Common::Configuration dgf_gridprovider_default_config()
{
  Common::Configuration config;
  config["type"]     = "xt.grid.gridprovider.dgf";
  config["filename"] = "dgf_1d_interval.dgf";
  return config;
}


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType>>::type
make_dgf_grid(const std::string& filename)
{
  return internal::create_dgf_grid<GridType>(filename);
}


template <class GridType>
typename std::enable_if<is_grid<GridType>::value, GridProvider<GridType>>::type
make_dgf_grid(const Common::Configuration& cfg = dgf_gridprovider_default_config())
{
  auto filename = cfg.get("filename", dgf_gridprovider_default_config().get<std::string>("filename"));
  return internal::create_dgf_grid<GridType>(filename);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_DGF_HH
