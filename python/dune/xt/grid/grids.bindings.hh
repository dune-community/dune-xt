// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2017 - 2018)

#ifndef DUNE_XT_GRID_GRIDS_BINDINGS_HH
#define DUNE_XT_GRID_GRIDS_BINDINGS_HH

#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/xt/grid/grids.hh>


namespace Dune {
namespace XT {
namespace Grid {

namespace DD {
template <class G>
class SubdomainGrid;
}

namespace bindings {


template <class G>
struct grid_name
{
  static std::string value()
  {
    return Common::Typename<G>::value() + "(missing specialization of grid_name)";
  }
};


template <class G>
struct grid_name<const G>
{
  static std::string value()
  {
    return grid_name<G>::value();
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

template <class G>
struct grid_name<DD::SubdomainGrid<G>>
{
  static std::string value()
  {
    return std::string("SubdomainGrid_") + grid_name<G>::value();
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


// not optimal
template<int dim, int dimworld, ALU3dGridElementType elType, class Comm>
struct grid_name<ALU3dGrid<dim, dimworld, elType, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_alugrid";
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
