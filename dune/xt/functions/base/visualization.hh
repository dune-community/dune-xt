// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_FUNCTIONS_BASE_VISUALIZATION_HH
#define DUNE_XT_FUNCTIONS_BASE_VISUALIZATION_HH

#include <dune/grid/io/file/vtk/function.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class GridViewType, size_t dimRange, size_t dimRangeCols, class RangeFieldImp>
class VisualizationAdapter : public VTKFunction<GridViewType>
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");

public:
  using EntityType = XT::Grid::extract_entity_t<GridViewType>;
  using GridFunctionType = GridFunctionInterface<EntityType, dimRange, dimRangeCols, RangeFieldImp>;

private:
  using LocalFunctionType = typename GridFunctionType::LocalFunctionType;
  using DomainType = typename LocalFunctionType::DomainType;
  using RangeType = typename LocalFunctionType::RangeType;

public:
  VisualizationAdapter(const GridFunctionType& localizable_function,
                       const std::string nm = "",
                       const XT::Common::Parameter& param = {})
    : local_function_(localizable_function.local_function())
    , name_(nm.empty() ? localizable_function.name() : nm)
    , param_(param)
  {
  }

  int ncomps() const override final
  {
    return helper<>::ncomps();
  }

  std::string name() const override final
  {
    return name_;
  }

  double evaluate(int comp, const EntityType& en, const DomainType& xx) const override final
  {
    local_function_->bind(en);
    return helper<>::evaluate(comp, local_function_->evaluate(xx, param_));
  }

private:
  template <size_t r_ = dimRange, size_t rC_ = dimRangeCols, bool anything = true>
  class helper
  {
  public:
    static int ncomps()
    {
      return 1;
    }

    static double evaluate(const int& /*comp*/, const RangeType& val)
    {
      return val.frobenius_norm();
    }
  }; // class helper<...>

  template <size_t r_, bool anything>
  class helper<r_, 1, anything>
  {
  public:
    static int ncomps()
    {
      return r_;
    }

    static double evaluate(const int& comp, const RangeType& val)
    {
      assert(comp >= 0);
      assert(comp < Common::numeric_cast<int>(r_));
      return val[comp];
    }
  }; // class helper<..., 1>

  mutable std::unique_ptr<LocalFunctionType> local_function_;
  const std::string name_;
  const XT::Common::Parameter& param_;
}; // class VisualizationAdapter


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_VISUALIZATION_HH
