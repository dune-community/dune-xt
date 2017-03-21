// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_GRID_GRIDS_BINDINGS_HH
#define DUNE_XT_GRID_GRIDS_BINDINGS_HH

#include <dune/xt/common/string.hh>

#include "grids.hh"


// this is used by other headers
typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> YASP_2D_EQUIDISTANT_OFFSET;
#if HAVE_DUNE_ALUGRID
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> ALU_2D_SIMPLEX_CONFORMING;
#endif


namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {


template <class G>
struct grid_name
{
  static_assert(AlwaysFalse<G>::value, "Please add a specialization for this grid!");

  static std::string value()
  {
    return "";
  }
};


template <int dim>
struct grid_name<YaspGrid<dim, EquidistantOffsetCoordinates<double, dim>>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_cube_yaspgrid";
  }
};


#if HAVE_DUNE_ALUGRID


template <int dim, class Comm>
struct grid_name<ALUGrid<dim, dim, simplex, conforming, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_simplex_aluconformgrid";
  }
};


#endif // HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDS_BINDINGS_HH
