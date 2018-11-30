// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018)

#ifndef DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_HH
#define DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_HH

#include "cube.hh"

#if DUNE_XT_WITH_PYTHON_BINDINGS


#define _DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS_DD_SUBDOMAINS_GRID(_prefix, _GRID)                         \
  _prefix Dune::XT::Grid::GridProvider<_GRID, Dune::XT::Grid::DD::SubdomainGrid<_GRID>>                                \
  Dune::XT::Grid::make_cube_dd_subdomains_grid<_GRID>(                                                                 \
      const Dune::FieldVector<typename _GRID::ctype, _GRID::dimension>&,                                               \
      const Dune::FieldVector<typename _GRID::ctype, _GRID::dimension>&,                                               \
      const std::array<unsigned int, _GRID::dimension>,                                                                \
      const unsigned int,                                                                                              \
      const std::array<unsigned int, _GRID::dimension>,                                                                \
      const std::array<unsigned int, _GRID::dimension>,                                                                \
      const size_t,                                                                                                    \
      const size_t)


#define DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS(_prefix, _GRID)                                             \
  _DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS_DD_SUBDOMAINS_GRID(_prefix, _GRID);                              \
  _prefix Dune::XT::Grid::GridProvider<_GRID> Dune::XT::Grid::make_cube_grid<_GRID>(                                   \
      const FieldVector<typename _GRID::ctype, _GRID::dimension>&,                                                     \
      const Dune::FieldVector<typename _GRID::ctype, _GRID::dimension>&,                                               \
      const std::array<unsigned int, _GRID::dimension>,                                                                \
      const unsigned int,                                                                                              \
      const std::array<unsigned int, _GRID::dimension>);                                                               \
  _prefix Dune::XT::Grid::GridProvider<_GRID> Dune::XT::Grid::make_cube_grid<_GRID>(const typename _GRID::ctype&,      \
                                                                                    const typename _GRID::ctype&,      \
                                                                                    const unsigned int,                \
                                                                                    const unsigned int,                \
                                                                                    const unsigned int);               \
  _prefix Dune::XT::Grid::GridProvider<_GRID> Dune::XT::Grid::make_cube_grid<_GRID>(                                   \
      const Dune::XT::Common::Configuration&)

#if HAVE_DUNE_ALUGRID
DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS(extern template, YASP_1D_EQUIDISTANT_OFFSET);
DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS(extern template, YASP_2D_EQUIDISTANT_OFFSET);
DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_FACTORY_METHODS(extern template, YASP_3D_EQUIDISTANT_OFFSET);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS


#endif // DUNE_XT_GRID_GRIDPROVIDER_CUBE_LIB_HH
