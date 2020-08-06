// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_XT_GRID_MAPPER_HH
#define DUNE_XT_GRID_MAPPER_HH

#include <dune/grid/common/mapper.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


template <class M, class E>
class SubEntityIndex
{
  static_assert(is_entity<E>::value, "");
  using G = extract_grid_t<E>;
  using D = typename G::ctype;
  static const constexpr size_t d = G::dimension;

  template <int cd = d, bool anything = true>
  struct subEntity
  {
    static size_t index(const M& mapper, const E& element, const int codim, const size_t i)
    {
      if (codim == cd) {
        DUNE_THROW_IF(i >= element.subEntities(codim),
                      Common::Exceptions::index_out_of_range,
                      "element.subEntities(" << codim << ") = " << element.subEntities(codim) << "\n   i = " << i);
        return Common::numeric_cast<size_t>(mapper.index(element.subEntity<cd>(Common::numeric_cast<int>(i))));
      } else
        return subEntity<cd - 1>::index(mapper, element, codim, i);
    } // ... center(...)
  }; // struct subEntity

  template <bool anything>
  struct subEntity<0, anything>
  {
    static size_t index(const M& mapper, const E& element, const int codim, const size_t i)
    {
      DUNE_THROW_IF(codim != 0, Common::Exceptions::internal_error, "This must not happen");
      DUNE_THROW_IF(i >= element.subEntities(0),
                    Common::Exceptions::index_out_of_range,
                    "element.subEntities(0) = " << element.subEntities(0) << "\n   i = " << i);
      return Common::numeric_cast<size_t>(mapper.index(element.subEntity<0>(Common::numeric_cast<int>(i))));
    }
  }; // struct subEntity

public:
  static size_t get(const M& mapper, const E& element, const int codim, const size_t i)
  {
    return subEntity<>::index(mapper, element, codim, i);
  }
}; // class SubEntityIndex


template <typename G,
          typename MapperImp,
          typename IndexType,
          int dim,
          class GridImp,
          template <int, int, class>
          class EntityImp>
size_t sub_entity_index(const Dune::Mapper<G, MapperImp, IndexType>& mapper,
                        const Dune::Entity<0, dim, GridImp, EntityImp> element,
                        const int codim,
                        const size_t i)
{
  return SubEntityIndex<Dune::Mapper<G, MapperImp, IndexType>, Dune::Entity<0, dim, GridImp, EntityImp>>::get(
      mapper, element, codim, i);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_MAPPER_HH
