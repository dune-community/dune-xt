// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017, 2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_HH
#define DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_HH

#include "provider.hh"

#if DUNE_XT_WITH_PYTHON_BINDINGS


#  define _DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_LAYER(_prefix, _GRID, _type, _backend)                                   \
    _prefix                                                                                                            \
        typename Dune::XT::Grid::Layer<_GRID, Dune::XT::Grid::Layers::_type, Dune::XT::Grid::Backends::_backend>::type \
        Dune::XT::Grid::GridProvider<_GRID>::layer<Dune::XT::Grid::Layers::_type, Dune::XT::Grid::Backends::_backend>( \
            const int) const

#  define _DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_ALL_LAYERS(_prefix, _GRID)                                               \
    _DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_PART(_prefix, _GRID);                                                          \
    _DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_LAYER(_prefix, _GRID, level, view);                                            \
    _DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_LAYER(_prefix, _GRID, leaf, view);

#  define DUNE_XT_GRID_PROVIDER_PROVIDER_LIB(_prefix, _GRID)                                                           \
    _prefix class Dune::XT::Grid::GridProvider<_GRID>;                                                                 \
    _DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_ALL_LAYERS(_prefix, _GRID)

#  if HAVE_DUNE_ALUGRID
DUNE_XT_GRID_PROVIDER_PROVIDER_LIB(extern template, ALU_2D_SIMPLEX_CONFORMING);
#  endif
DUNE_XT_GRID_PROVIDER_PROVIDER_LIB(extern template, YASP_1D_EQUIDISTANT_OFFSET);
DUNE_XT_GRID_PROVIDER_PROVIDER_LIB(extern template, YASP_2D_EQUIDISTANT_OFFSET);
DUNE_XT_GRID_PROVIDER_PROVIDER_LIB(extern template, YASP_3D_EQUIDISTANT_OFFSET);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS


#endif // DUNE_XT_GRID_PROVIDER_PROVIDER_LIB_HH
