// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_XT_GRID_ELEMENT_HH
#define DUNE_XT_GRID_ELEMENT_HH

#include <dune/grid/common/entity.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune::XT::Grid {


template <class E>
class SubEntityCenter
{
  static_assert(is_entity<E>::value, "");
  using G = extract_grid_t<E>;
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;

  template <int cd = d, bool anything = true>
  struct subEntity
  {
    static FieldVector<D, d> center(const E& element, const int codim, const size_t i)
    {
      if (codim == cd) {
        DUNE_THROW_IF(i >= element.subEntities(codim),
                      Common::Exceptions::index_out_of_range,
                      "element.subEntities(" << codim << ") = " << element.subEntities(codim) << "\n   i = " << i);
        return element.template subEntity<cd>(Common::numeric_cast<int>(i)).geometry().center();
      } else
        return subEntity<cd - 1>::center(element, codim, i);
    } // ... center(...)
  }; // struct subEntity

  template <bool anything>
  struct subEntity<0, anything>
  {
    static FieldVector<D, d> center(const E& element, const int codim, const size_t i)
    {
      DUNE_THROW_IF(codim != 0, Common::Exceptions::internal_error, "This must not happen");
      return element.template subEntity<0>(Common::numeric_cast<int>(i)).geometry().center();
    }
  }; // struct subEntity

public:
  static FieldVector<D, d> get(const E& element, const int codim, const size_t i)
  {
    return subEntity<>::center(element, codim, i);
  }
}; // class SubEntityCenter


template <int dim, class GridImp, template <int, int, class> class EntityImp>
FieldVector<typename GridImp::ctype, GridImp::dimension>
sub_entity_center(const Dune::Entity<0, dim, GridImp, EntityImp> element, const int codim, const size_t i)
{
  return SubEntityCenter<Dune::Entity<0, dim, GridImp, EntityImp>>::get(element, codim, i);
}


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_ELEMENT_HH
