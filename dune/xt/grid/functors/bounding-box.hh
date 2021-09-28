// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018, 2020)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_GRID_FUNCTORS_BOUNDING_BOX_HH
#define DUNE_XT_GRID_FUNCTORS_BOUNDING_BOX_HH

#include <dune/xt/common/ranges.hh>
#include <dune/xt/grid/information.hh>

#include "interfaces.hh"

namespace Dune::XT::Grid {
namespace internal {


template <class D, size_t d>
struct MinMaxCoordinateFunctorReduction
{
  using T = std::pair<FieldVector<D, d>, FieldVector<D, d>>;

  constexpr T operator()(const T& lhs, const T& rhs) const
  {
    T result;
    for (size_t ii = 0; ii < d; ++ii) {
      result.first[ii] = std::min(lhs.first[ii], rhs.first[ii]);
      result.second[ii] = std::max(lhs.second[ii], rhs.second[ii]);
    }
    return result;
  }
}; // struct MinMaxCoordinateFunctorReduction


} // namespace internal


/**
 * \brief     Functor for a \ref GridWalk calculating minima and maxima of grids' coordinates
 * \attention Broken!
 * \todo Bring up to date!
 **/
template <class GV>
class MinMaxCoordinateFunctor
  : public ElementFunctor<GV>
  , public Common::ThreadResultPropagator<
        MinMaxCoordinateFunctor<GV>,
        std::pair<FieldVector<typename GV::ctype, GV::dimension>, FieldVector<typename GV::ctype, GV::dimension>>,
        internal::MinMaxCoordinateFunctorReduction<typename GV::ctype, GV::dimension>>
{
  using BaseType = ElementFunctor<GV>;

  using D = typename GV::ctype;
  static constexpr size_t d = GV::dimension;
  using VectorType = FieldVector<D, d>;
  using Propagator = Common::ThreadResultPropagator<MinMaxCoordinateFunctor<GV>,
                                                    std::pair<VectorType, VectorType>,
                                                    internal::MinMaxCoordinateFunctorReduction<D, d>>;

public:
  MinMaxCoordinateFunctor()
    : Propagator(this)
    , bounding_box_(
          std::make_pair(VectorType(std::numeric_limits<D>::max()), VectorType(std::numeric_limits<D>::min())))
  {}

  MinMaxCoordinateFunctor(const MinMaxCoordinateFunctor& other) = default;

  void prepare() override final
  {
    bounding_box_ =
        std::make_pair(VectorType(std::numeric_limits<D>::max()), VectorType(std::numeric_limits<D>::min()));
  }

  void apply_local(const typename BaseType::ElementType& element) override final
  {
    const auto& geo = element.geometry();
    for (auto i : Common::value_range(geo.corners())) {
      const auto corner = geo.corner(i);
      for (auto k : Common::value_range(d)) {
        bounding_box_.first[k] = std::min(bounding_box_.first[k], corner[k]);
        bounding_box_.second[k] = std::max(bounding_box_.second[k], corner[k]);
      }
    }
  } // ... apply_local(...)

  std::pair<VectorType, VectorType> result() const
  {
    return bounding_box_;
  }

  void finalize() override final
  {
    Propagator::finalize_imp();
  }

  BaseType* copy() override final
  {
    return Propagator::copy_imp();
  }

protected:
  void set_result(std::pair<VectorType, VectorType> res)
  {
    bounding_box_ = res;
  }

private:
  friend Propagator;

  std::pair<VectorType, VectorType> bounding_box_;
}; // class MinMaxCoordinateFunctor


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_FUNCTORS_BOUNDING_BOX_HH
