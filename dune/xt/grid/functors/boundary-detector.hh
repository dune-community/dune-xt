// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_GRID_FUNCTORS_BOUNDARY_DETECTOR_HH
#define DUNE_XT_GRID_FUNCTORS_BOUNDARY_DETECTOR_HH


#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/common/parallel/threadstorage.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class GL>
class BoundaryDetectorFunctor : public IntersectionFunctor<GL>,
                                public Common::ThreadResultPropagator<BoundaryDetectorFunctor<GL>, size_t>
{
  using BaseType = IntersectionFunctor<GL>;
  using Propagator = Common::ThreadResultPropagator<BoundaryDetectorFunctor<GL>, size_t>;
  friend Propagator;

public:
  using ResultType = size_t;
  using typename BaseType::GridViewType;
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;

  /**
   * \attention Takes ownership of boundary_type, do not delete manually!
   */
  BoundaryDetectorFunctor(const BoundaryInfo<IntersectionType>& boundary_info, BoundaryType* boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
    , found_(0)
  {
  }

  BoundaryDetectorFunctor(const BoundaryInfo<IntersectionType>& boundary_info,
                          const std::shared_ptr<BoundaryType>& boundary_type)
    : boundary_info_(boundary_info)
    , boundary_type_(boundary_type)
    , found_(0)
  {
  }

  void prepare() override final
  {
    found_ = 0;
  }

  ResultType compute_locally(const IntersectionType& intersection,
                             const ElementType& /*inside_element*/,
                             const ElementType& /*outside_element*/)
  {
    return boundary_info_.type(intersection) == *boundary_type_;
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    found_ += compute_locally(intersection, inside_element, outside_element);
  }

  ResultType result() const
  {
    return found_;
  }

  void finalize() override
  {
    Propagator::finalize_imp();
  }

  std::shared_ptr<BaseType> copy() override
  {
    return Propagator::copy_imp();
  }

protected:
  void set_result(ResultType res)
  {
    found_ = res;
  }

private:
  const BoundaryInfo<IntersectionType>& boundary_info_;
  const std::shared_ptr<BoundaryType> boundary_type_;
  size_t found_;


}; // class BoundaryDetectorFunctor


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_FUNCTORS_BOUNDARY_DETECTOR_HH
