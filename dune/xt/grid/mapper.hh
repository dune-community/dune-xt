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

namespace Dune::XT::Grid {


template <typename G,
          typename MapperImp,
          typename IndexType,
          int dim,
          class GridImp,
          template <int, int, class>
          class EntityImp,
          int cd = GridImp::dimension>
size_t sub_entity_index(const Dune::Mapper<G, MapperImp, IndexType>& mapper,
                        const Dune::Entity<0, dim, GridImp, EntityImp> element,
                        const int codim,
                        const size_t i)
{
  if constexpr (cd == 0) {
    DUNE_THROW_IF(codim != 0, Common::Exceptions::internal_error, "This must not happen");
    DUNE_THROW_IF(i >= element.subEntities(0),
                  Common::Exceptions::index_out_of_range,
                  "element.subEntities(0) = " << element.subEntities(0) << "\n   i = " << i);
    return Common::numeric_cast<size_t>(mapper.index(element.template subEntity<0>(Common::numeric_cast<int>(i))));
  } else {
    if (codim == cd) {
      DUNE_THROW_IF(i >= element.subEntities(codim),
                    Common::Exceptions::index_out_of_range,
                    "element.subEntities(" << codim << ") = " << element.subEntities(codim) << "\n   i = " << i);
      return Common::numeric_cast<size_t>(mapper.index(element.template subEntity<cd>(Common::numeric_cast<int>(i))));
    } else
      return sub_entity_index<G, MapperImp, IndexType, dim, GridImp, EntityImp, cd - 1>(mapper, element, codim, i);
  }
}


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_MAPPER_HH
