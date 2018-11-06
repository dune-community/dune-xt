// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk       (2018)

#include <dune/xt/common/test/main.hxx>

#include <dune/xt/grid/gridprovider/eoc.hh>
#include <dune/xt/grid/grids.hh>


{% for name, type in config.all_grids %}

GTEST_TEST(EocProvider_{{name}}, layers)
{
  using Level DUNE_UNUSED = Dune::XT::Grid::LevelBasedEOCGridProvider<{{type}}>;
  using DD DUNE_UNUSED = Dune::XT::Grid::DdSubdomainsBasedEOCGridProvider<{{type}}>;
  using Leaf DUNE_UNUSED = Dune::XT::Grid::LeafBasedEOCGridProvider<{{type}}>;
}

{% endfor %}
