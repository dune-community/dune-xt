// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

#ifndef DUNE_XT_GRID_LAYERS_BINDINGS_HH
#define DUNE_XT_GRID_LAYERS_BINDINGS_HH

#include <dune/xt/grid/layers.hh>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <Backends backend>
struct backend_name
{
  static_assert(AlwaysFalse<typename XT::Grid::internal::backend_dependent_typename<backend>::type>::value,
                "Please add a specialization for this backend!");

  static std::string value()
  {
    return "";
  }
};


template <>
struct backend_name<Backends::view>
{
  static std::string value()
  {
    return "view";
  }
};


template <Layers layer>
struct layer_name
{
  DXT_DEPRECATED_MSG("use layer_names[layer] directly. 2018/7/2") static std::string value()
  {
    return layer_names[layer];
  }
};


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_LAYERS_BINDINGS_HH
