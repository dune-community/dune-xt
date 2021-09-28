// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2020)
//   René Fritze     (2017 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_GRID_GRIDS_BINDINGS_HH
#define DUNE_XT_GRID_GRIDS_BINDINGS_HH

#include <dune/xt/common/string.hh>
#include <dune/xt/common/tuple.hh>
#include <dune/xt/grid/grids.hh>


namespace Dune::XT::Grid::bindings {


template <class G>
struct grid_name
{
  static std::string value()
  {
    static_assert(AlwaysFalse<G>::value, "missing specialization of grid_name");
    return "";
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


template <>
struct grid_name<Dune::OneDGrid>
{
  static std::string value()
  {
    return "1d_simplex_onedgrid";
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

template <int dim, class Comm>
struct grid_name<ALU3dGrid<dim, dim, ALU3dGridElementType::tetra, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_tetrahedral_alugrid";
  }
};

template <int dim, class Comm>
struct grid_name<ALU3dGrid<dim, dim, ALU3dGridElementType::hexa, Comm>>
{
  static std::string value()
  {
    return Common::to_string(dim) + "d_hexahedral_alugrid";
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


/// \attention The following choices are on purpose: only two variants per dim, one cube one simplex.
///            In particular the grid_name<G>::value needs to be unique for all alugrid variants and the
///            make_...grid methods need to be more general if we extend the choice of grids here!

using Available1dGridTypes = std::tuple<ONED_1D>;

using Available2dGridTypes = std::tuple<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                        ,
                                        ALU_2D_SIMPLEX_CONFORMING
#endif
                                        >;

using Available3dGridTypes = std::tuple<YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                        ,
                                        ALU_3D_SIMPLEX_CONFORMING
#endif
                                        >;

using AvailableGridTypes = Common::tuple_cat_t<Available1dGridTypes, Available2dGridTypes, Available3dGridTypes>;


} // namespace Dune::XT::Grid::bindings

#endif // DUNE_XT_GRID_GRIDS_BINDINGS_HH
