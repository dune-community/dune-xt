// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2016, 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_FUNCTIONS_VISUALIZATION_HH
#define DUNE_XT_FUNCTIONS_VISUALIZATION_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/grid/io/file/vtk/function.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {


template <class GridViewType, size_t dimRange, size_t dimRangeCols>
class VisualizationAdapterFunction : public VTKFunction<GridViewType>
{
public:
  using EntityType = XT::Grid::extract_entity_t<GridViewType>;

  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, double, dimRange, dimRangeCols>
      FunctionType;

  VisualizationAdapterFunction(const FunctionType& function,
                               const std::string nm = "",
                               const XT::Common::Parameter& param = {})
    : function_(function)
    , tmp_value_(0)
    , name_(nm)
    , param_(param)
  {
  }

private:
  template <size_t r, size_t rC, bool anything = true>
  class Call
  {
  public:
    static int ncomps()
    {
      return 1;
    }

    static double evaluate(const int& /*comp*/, const typename FunctionType::RangeType& val)
    {
      return val.frobenius_norm();
    }
  }; // class Call

  template <size_t r, bool anything>
  class Call<r, 1, anything>
  {
  public:
    static int ncomps()
    {
      return r;
    }

    static double evaluate(const int& comp, const typename FunctionType::RangeType& val)
    {
      return val[comp];
    }
  }; // class Call< ..., 1, ... >

public:
  virtual int ncomps() const override final
  {
    return Call<dimRange, dimRangeCols>::ncomps();
  }

  virtual std::string name() const override final
  {
    if (name_.empty())
      return function_.name();
    else
      return name_;
  }

  virtual double evaluate(int comp, const EntityType& en, const DomainType& xx) const override final
  {
    assert(comp >= 0);
    assert(comp < boost::numeric_cast<int>(dimRange));
    const auto local_func = function_.local_function(en);
    local_func->evaluate(xx, tmp_value_, param_);
    return Call<dimRange, dimRangeCols>::evaluate(comp, tmp_value_);
  }

private:
  const FunctionType& function_;
  mutable typename FunctionType::RangeType tmp_value_;
  const std::string name_;
  const XT::Common::Parameter& param_;
}; // class VisualizationAdapterFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_VISUALIZATION_HH
