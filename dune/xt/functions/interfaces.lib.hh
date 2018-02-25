// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_INTERFACE_LIB_HH
#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_HH

#include <dune/xt/functions/interfaces.hh>

#if DUNE_XT_WITH_PYTHON_BINDINGS


#define _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_DD_SUBDOMAIN(_p, _G, _t, _b, _R, _r, _rC)                           \
  _p void Dune::XT::Functions::LocalizableFunctionInterface<                                                           \
      typename _G::template Codim<0>::Entity,                                                                          \
      typename _G::ctype,                                                                                              \
      _G::dimension,                                                                                                   \
      _R,                                                                                                              \
      _r,                                                                                                              \
      _rC>::visualize(const typename Dune::XT::Grid::Layer<_G,                                                         \
                                                           Dune::XT::Grid::Layers::_t,                                 \
                                                           Dune::XT::Grid::Backends::_b,                               \
                                                           Dune::XT::Grid::DD::SubdomainGrid<_G>>::type&,              \
                      const std::string,                                                                               \
                      const bool,                                                                                      \
                      const VTK::OutputType) const

#define _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE(_p, _G, _t, _b, _R, _r, _rC)                                        \
  _p void Dune::XT::Functions::LocalizableFunctionInterface<typename _G::template Codim<0>::Entity,                    \
                                                            typename _G::ctype,                                        \
                                                            _G::dimension,                                             \
                                                            _R,                                                        \
                                                            _r,                                                        \
                                                            _rC>::                                                     \
      visualize(                                                                                                       \
          const typename Dune::XT::Grid::Layer<_G, Dune::XT::Grid::Layers::_t, Dune::XT::Grid::Backends::_b>::type&,   \
          const std::string,                                                                                           \
          const bool,                                                                                                  \
          const VTK::OutputType) const

#if HAVE_DUNE_FEM
#define _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_PART(_p, _G, _R, _r, _rC)                                           \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE(_p, _G, adaptive_leaf, view, _R, _r, _rC);                                \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE(_p, _G, leaf, part, _R, _r, _rC);                                         \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE(_p, _G, level, part, _R, _r, _rC);                                        \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_DD_SUBDOMAIN(_p, _G, dd_subdomain, view, _R, _r, _rC)
#else
#define _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_PART(_p, _G, _R, _r, _rC)
#endif

#define _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_ALL_LAYERS(_p, _G, _R, _r, _rC)                                     \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_PART(_p, _G, _R, _r, _rC);                                                \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE(_p, _G, leaf, view, _R, _r, _rC);                                         \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE(_p, _G, level, view, _R, _r, _rC)

#define _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, _G, _R, _r, _rC)                                                          \
  _p class Dune::XT::Functions::LocalfunctionSetInterface<typename _G::template Codim<0>::Entity,                      \
                                                          typename _G::ctype,                                          \
                                                          _G::dimension,                                               \
                                                          _R,                                                          \
                                                          _r,                                                          \
                                                          _rC>;                                                        \
  _p class Dune::XT::Functions::                                                                                       \
      LocalfunctionInterface<typename _G::template Codim<0>::Entity, typename _G::ctype, _G::dimension, _R, _r, _rC>;  \
  _p class Dune::XT::Functions::LocalizableFunctionInterface<typename _G::template Codim<0>::Entity,                   \
                                                             typename _G::ctype,                                       \
                                                             _G::dimension,                                            \
                                                             _R,                                                       \
                                                             _r,                                                       \
                                                             _rC>;                                                     \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB_VISUALIZE_ALL_LAYERS(_p, _G, _R, _r, _rC)

#if HAVE_DUNE_ALUGRID
#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_ALU_2D(_p)                                                                     \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, ALU_2D_SIMPLEX_CONFORMING, double, 1, 1);                                       \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, ALU_2D_SIMPLEX_CONFORMING, double, 2, 2)
#else
#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_ALU_2D(_p)
#endif

#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_YASP_1D(_p)                                                                    \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, YASP_1D_EQUIDISTANT_OFFSET, double, 1, 1)

#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_YASP_2D(_p)                                                                    \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, YASP_2D_EQUIDISTANT_OFFSET, double, 1, 1);                                      \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, YASP_2D_EQUIDISTANT_OFFSET, double, 2, 2)

#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_YASP_3D(_p)                                                                    \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, YASP_3D_EQUIDISTANT_OFFSET, double, 1, 1);                                      \
  _DUNE_XT_FUNCTIONS_INTERFACE_LIB(_p, YASP_3D_EQUIDISTANT_OFFSET, double, 3, 3)

#define DUNE_XT_FUNCTIONS_INTERFACE_LIB_ALL_GRIDS(_p)                                                                  \
  DUNE_XT_FUNCTIONS_INTERFACE_LIB_ALU_2D(_p);                                                                          \
  DUNE_XT_FUNCTIONS_INTERFACE_LIB_YASP_1D(_p);                                                                         \
  DUNE_XT_FUNCTIONS_INTERFACE_LIB_YASP_2D(_p);                                                                         \
  DUNE_XT_FUNCTIONS_INTERFACE_LIB_YASP_3D(_p)

DUNE_XT_FUNCTIONS_INTERFACE_LIB_ALL_GRIDS(extern template);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS


#endif // DUNE_XT_FUNCTIONS_INTERFACE_LIB_HH
