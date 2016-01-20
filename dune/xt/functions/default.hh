// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2015)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_FUNCTIONS_VISUALIZATION_HH
#define DUNE_XT_FUNCTIONS_VISUALIZATION_HH

#include <boost/numeric/conversion/cast.hpp>

#if HAVE_DUNE_GRID
#include <dune/grid/io/file/vtk/function.hh>
#endif

#include <dune/xt/common/float_cmp.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {

#if HAVE_DUNE_GRID

template <class GridViewType, size_t dimRange, size_t dimRangeCols>
class VisualizationAdapter : public VTKFunction<GridViewType>
{
public:
  typedef typename GridViewType::template Codim<0>::Entity EntityType;

  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, double, dimRange, dimRangeCols>
      FunctionType;

  VisualizationAdapter(const FunctionType& function, const std::string nm = "")
    : function_(function)
    , tmp_value_(0)
    , name_(nm)
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
    local_func->evaluate(xx, tmp_value_);
    return Call<dimRange, dimRangeCols>::evaluate(comp, tmp_value_);
  }

private:
  const FunctionType& function_;
  mutable typename FunctionType::RangeType tmp_value_;
  const std::string name_;
}; // class VisualizationAdapter

#endif // HAVE_DUNE_GRID

} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_VISUALIZATION_HH
