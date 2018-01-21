// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2017)
//   Kirsten Weber   (2012)
//   Rene Milk       (2012 - 2013, 2015 - 2016, 2018)
//   Tobias Leibner  (2014 - 2016)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_GMSH_HH
#define DUNE_XT_GRID_GRIDPROVIDER_GMSH_HH

#include <memory>
#include <type_traits>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>

#include <dune/xt/grid/grids.hh>

#include "provider.hh"

namespace Dune {
namespace XT {
namespace Grid {


static inline std::string gmsh_gridprovider_id()
{
  return "xt.grid.gridprovider.gmsh";
}


static inline Common::Configuration gmsh_gridprovider_default_config()
{
  Common::Configuration config;
  config["type"] = gmsh_gridprovider_id();
  config["filename"] = "g.msh";
  return config;
}


template <class GridType>
class GmshGridProviderFactory;


template <int dim, class Coordinates>
class GmshGridProviderFactory<Dune::YaspGrid<dim, Coordinates>>
{
public:
  static const bool available = false;
};


template <class GridType>
class GmshGridProviderFactory
{
public:
  static const bool available = false;

  static std::string static_id()
  {
    return gmsh_gridprovider_id();
  }

  static Common::Configuration default_config()
  {
    auto cfg = gmsh_gridprovider_default_config();
#if HAVE_DUNE_ALUGRID
    if (std::is_same<ALUGrid<2, 2, simplex, conforming>, GridType>::value
        || std::is_same<ALUGrid<2, 2, simplex, nonconforming>, GridType>::value) {
      cfg["filename"] = "gmsh_2d_simplices.msh";
    }
#endif // HAVE_DUNE_ALUGRID
    return cfg;
  }

  static GridProvider<GridType> create(const std::string& filename)
  {
    return GridProvider<GridType>(GmshReader<GridType>::read(filename));
  }

  static GridProvider<GridType> create(const Common::Configuration& cfg = default_config())
  {
    return create(cfg.get("filename", default_config().template get<std::string>("filename")));
  }
}; // class GmshGridProviderFactory


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDPROVIDER_GMSH_HH
