// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_GRID_FUNCTORS_BOUNDING_BOX_HH
#define DUNE_XT_GRID_FUNCTORS_BOUNDING_BOX_HH

#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/information.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


/** \brief Functor for a \ref GridWalk calculating minima and maxima of entities' coordinates
 **/
template <class GridViewType>
struct MinMaxCoordinateFunctor : public ElementFunctor<GridViewType>
{
  typedef ElementFunctor<GridViewType> BaseType;
  typedef typename BaseType::ElementType::Geometry EntityGeometryType;
  typedef typename EntityGeometryType::ctype ctype;
  typedef FieldVector<ctype, EntityGeometryType::coorddimension> VectorType;

  MinMaxCoordinateFunctor()
    : minima_(VectorType(std::numeric_limits<ctype>::max()))
    , maxima_(VectorType(std::numeric_limits<ctype>::min()))
  {}

  void apply_local(const typename BaseType::ElementType& element) override final
  {
    const auto& geo = element.geometry();
    for (auto i : Common::value_range(geo.corners())) {
      for (auto k : Common::value_range(EntityGeometryType::coorddimension)) {
        minima_[k] = std::min(minima_[k], geo.corner(i)[k]);
        maxima_[k] = std::max(maxima_[k], geo.corner(i)[k]);
      }
    }
  }

  ElementFunctor<GridViewType>* copy() override
  {
    return new MinMaxCoordinateFunctor<GridViewType>(*this);
  }

  VectorType minima_;
  VectorType maxima_;
};


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_BOUNDING_BOX_HH
