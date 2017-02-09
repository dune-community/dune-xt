// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "expression.hh"

namespace Dune {
namespace XT {
namespace Functions {


#define _GRID YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>
template class ExpressionFunction<
    typename XT::Grid::Entity<
        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
    typename _GRID::ctype,
    _GRID::dimension,
    double,
    1,
    1>;
template class ExpressionFunction<
    typename XT::Grid::Entity<
        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
    typename _GRID::ctype,
    _GRID::dimension,
    double,
    2,
    1>;
// template class ExpressionFunction<
//    typename XT::Grid::Entity<
//        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
//    typename _GRID::ctype,
//    _GRID::dimension,
//    double,
//    2,
//    2>;
template class ExpressionFunction<
    typename XT::Grid::Entity<
        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
    typename _GRID::ctype,
    _GRID::dimension,
    double,
    3,
    1>;
// template class ExpressionFunction<
//    typename XT::Grid::Entity<
//        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
//    typename _GRID::ctype,
//    _GRID::dimension,
//    double,
//    3,
//    3>;
#undef _GRID

#if HAVE_ALUGRID
#define _GRID ALUGrid<2, 2, simplex, conforming>
template class ExpressionFunction<
    typename XT::Grid::Entity<
        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
    typename _GRID::ctype,
    _GRID::dimension,
    double,
    1,
    1>;
template class ExpressionFunction<
    typename XT::Grid::Entity<
        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
    typename _GRID::ctype,
    _GRID::dimension,
    double,
    2,
    1>;
// template class ExpressionFunction<
//    typename XT::Grid::Entity<
//        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
//    typename _GRID::ctype,
//    _GRID::dimension,
//    double,
//    2,
//    2>;
template class ExpressionFunction<
    typename XT::Grid::Entity<
        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
    typename _GRID::ctype,
    _GRID::dimension,
    double,
    3,
    1>;
// template class ExpressionFunction<
//    typename XT::Grid::Entity<
//        typename XT::Grid::Layer<_GRID, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>::type,
//    typename _GRID::ctype,
//    _GRID::dimension,
//    double,
//    3,
//    3>;
#undef _GRID
#endif // HAVE_ALUGRID


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
