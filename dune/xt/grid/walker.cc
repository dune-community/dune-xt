// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#include "walker.hh"

#if HAVE_DUNE_PYBINDXI

namespace Dune {
namespace XT {
namespace Grid {


template class Walker<
    typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::leaf, Backends::view>::type>;
template class Walker<
    typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::level, Backends::view>::type>;
#if HAVE_DUNE_FEM
template class Walker<
    typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::leaf, Backends::part>::type>;
template class Walker<
    typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::level, Backends::part>::type>;
#endif // HAVE_DUNE_FEM


#if HAVE_ALUGRID
template class Walker<typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::leaf, Backends::view>::type>;
template class Walker<typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::level, Backends::view>::type>;
#if HAVE_DUNE_FEM
template class Walker<typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::leaf, Backends::part>::type>;
template class Walker<typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::level, Backends::part>::type>;
#endif // HAVE_DUNE_FEM
#endif // HAVE_ALUGRID


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
