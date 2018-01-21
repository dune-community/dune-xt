// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)

#ifndef DUNE_XT_GRID_GRIDS_BINDINGS_HH
#define DUNE_XT_GRID_GRIDS_BINDINGS_HH

#include <dune/xt/common/string.hh>
#include <dune/xt/grid/layers.hh>

#include "grids.hh"


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

template <int dim, class Comm>
struct grid_name<ALUGrid<dim, dim, cube, nonconforming, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_cube_alunonconformgrid";
  }
};

template <int dim, class Comm>
struct grid_name<ALUGrid<dim, dim, simplex, nonconforming, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_simplex_alunonconformgrid";
  }
};

template <int dim, class Comm>
struct grid_name<ALUGrid<dim, dim, cube, conforming, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_cube_aluconformgrid";
  }
};

#endif // HAVE_DUNE_ALUGRID
#if HAVE_ALBERTA


template <int dim>
struct grid_name<Dune::AlbertaGrid<dim, dim>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_simplex_albertagrid";
  }
};


#endif // HAVE_ALBERTA
#if HAVE_DUNE_UGGRID || HAVE_UG


template <int dim>
struct grid_name<UGGrid<dim>>
{
  static std::string value()
  { // the "simplex" is due to the fact that our grid provider currently only creates ug grids with simplices
    return Common::to_string(dim) + "d_simplex_uggrid";
  }
};


#endif // HAVE_DUNE_UGGRID || HAVE_UG

} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_GRIDS_BINDINGS_HH
