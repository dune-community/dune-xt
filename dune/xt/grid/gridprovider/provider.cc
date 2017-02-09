// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#include "provider.hh"

#if HAVE_DUNE_PYBINDXI

namespace Dune {
namespace XT {
namespace Grid {


template class GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::level, Backends::view>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::level<Backends::view>(const int) const;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::leaf, Backends::view>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::leaf<Backends::view>() const;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::leaf, Backends::view>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::layer<Layers::leaf, Backends::view>(
    const int) const;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::level, Backends::view>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::layer<Layers::level, Backends::view>(
    const int) const;

#if HAVE_DUNE_FEM

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::level, Backends::part>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::level<Backends::part>(const int lvl) const;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::leaf, Backends::part>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::leaf<Backends::part>() const;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::leaf, Backends::part>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::layer<Layers::leaf, Backends::part>(
    const int) const;

template typename Layer<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>, Layers::level, Backends::part>::type
GridProvider<YaspGrid<2, EquidistantOffsetCoordinates<double, 2>>>::layer<Layers::level, Backends::part>(
    const int) const;

#endif // HAVE_DUNE_FEM


#if HAVE_ALUGRID

template class GridProvider<ALUGrid<2, 2, simplex, conforming>>;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::level, Backends::view>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::level<Backends::view>(const int) const;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::leaf, Backends::view>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::leaf<Backends::view>() const;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::leaf, Backends::view>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::layer<Layers::leaf, Backends::view>(const int) const;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::level, Backends::view>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::layer<Layers::level, Backends::view>(const int) const;

#if HAVE_DUNE_FEM

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::level, Backends::part>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::level<Backends::part>(const int lvl) const;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::leaf, Backends::part>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::leaf<Backends::part>() const;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::leaf, Backends::part>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::layer<Layers::leaf, Backends::part>(const int) const;

template typename Layer<ALUGrid<2, 2, simplex, conforming>, Layers::level, Backends::part>::type
GridProvider<ALUGrid<2, 2, simplex, conforming>>::layer<Layers::level, Backends::part>(const int) const;

#endif // HAVE_DUNE_FEM
#endif // HAVE_ALUGRID

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
