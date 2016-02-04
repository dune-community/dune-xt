// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_TYPE_TRAITS_HH
#define DUNE_XT_GRID_TYPE_TRAITS_HH

#include <type_traits>

#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class T>
struct is_intersection : public std::false_type
{
};

template <class G, class I>
struct is_intersection<Dune::Intersection<G, I>> : public std::true_type
{
};


template <class T>
struct is_grid_view : public std::false_type
{
};

template <class T>
struct is_grid_view<Dune::GridView<T>> : public std::true_type
{
};


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_TYPE_TRAITS_HH
