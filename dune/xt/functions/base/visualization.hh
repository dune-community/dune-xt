// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2019)
//   René Fritze     (2013 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_BASE_VISUALIZATION_HH
#define DUNE_XT_FUNCTIONS_BASE_VISUALIZATION_HH

#include <algorithm>

#include <dune/grid/io/file/vtk/function.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


// forward
template <class Element, size_t rangeDim, size_t rangeDimCols, class RangeField>
class GridFunctionInterface;


template <size_t r, size_t rC, class R = double>
class VisualizerInterface
{
public:
  using RangeType = typename RangeTypeSelector<R, r, rC>::type;

  virtual ~VisualizerInterface(){};

  virtual int ncomps() const = 0;

  virtual double evaluate(const int& comp, const RangeType& val) const = 0;
}; // class VisualizerInterface


// visualizes all components of the function
template <size_t r, size_t rC, class R = double>
class DefaultVisualizer : public VisualizerInterface<r, rC, R>
{
  using BaseType = VisualizerInterface<r, rC, R>;

public:
  using typename BaseType::RangeType;

  int ncomps() const override final
  {
    return helper<>::ncomps();
  }

  double evaluate(const int& comp, const RangeType& val) const override final
  {
    return helper<>::evaluate(comp, val);
  }

private:
  template <size_t r_ = r, size_t rC_ = rC, bool anything = true>
  struct helper
  {
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
  struct helper<r_, 1, anything>
  {
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
  }; // struct helper<..., 1>
}; // class DefaultVisualizer

template <size_t r, size_t rC, class R = double>
const DefaultVisualizer<r, rC, R>& default_visualizer()
{
  static DefaultVisualizer<r, rC, R> default_visualizer_instance_;
  return default_visualizer_instance_;
}

template <size_t r, size_t rC = 1, class R = double>
class SumVisualizer : public VisualizerInterface<r, rC, R>
{
  static_assert(rC == 1, "Not implemented for rC > 1");
  using BaseType = VisualizerInterface<r, rC, R>;

public:
  using typename BaseType::RangeType;

  int ncomps() const override final
  {
    return 1;
  }

  double evaluate(const int& DXTC_DEBUG_ONLY(comp), const RangeType& val) const override final
  {
    assert(comp == 0);
    return std::accumulate(val.begin(), val.end(), 0.);
  }
}; // class SumVisualizer


template <size_t r, size_t rC = 1, class R = double>
class ComponentVisualizer : public VisualizerInterface<r, rC, R>
{
  static_assert(rC == 1, "Not implemented for rC > 1");
  using BaseType = VisualizerInterface<r, rC, R>;

public:
  using typename BaseType::RangeType;

  ComponentVisualizer(const int comp)
    : comp_(comp)
  {}

  int ncomps() const override final
  {
    return 1;
  }

  double evaluate(const int& comp, const RangeType& val) const override final
  {
    return val[comp];
  }

private:
  int comp_;
}; // class ComponentVisualizer


template <size_t r, size_t rC = 1, class R = double>
class GenericVisualizer : public VisualizerInterface<r, rC, R>
{
  using BaseType = VisualizerInterface<r, rC, R>;

public:
  using typename BaseType::RangeType;
  using EvalType = std::function<double(const int comp, const RangeType& val)>;

  GenericVisualizer(const int ncomps, EvalType eval)
    : ncomps_(ncomps)
    , eval_(eval)
  {}

  int ncomps() const override final
  {
    return ncomps_;
  }

  double evaluate(const int& comp, const RangeType& val) const override final
  {
    return eval_(comp, val);
  }

private:
  int ncomps_;
  EvalType eval_;
}; // class GenericVisualizer


template <class GridViewType, size_t range_dim, size_t range_dim_cols, class RangeField>
class VisualizationAdapter : public VTKFunction<GridViewType>
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");

public:
  using EntityType = XT::Grid::extract_entity_t<GridViewType>;
  using GridFunctionType = GridFunctionInterface<EntityType, range_dim, range_dim_cols, RangeField>;

private:
  using LocalFunctionType = typename GridFunctionType::LocalFunctionType;
  using DomainType = typename LocalFunctionType::DomainType;

public:
  VisualizationAdapter(const GridFunctionType& localizable_function,
                       const VisualizerInterface<range_dim, range_dim_cols, RangeField>& visualizer,
                       const std::string nm = "",
                       const XT::Common::Parameter& param = {})
    : local_function_(localizable_function.local_function())
    , visualizer_(visualizer)
    , name_(nm.empty() ? localizable_function.name() : nm)
    , param_(param)
  {}

  VisualizationAdapter(const GridFunctionType& localizable_function,
                       const std::string nm = "",
                       const XT::Common::Parameter& param = {})
    : local_function_(localizable_function.local_function())
    , visualizer_(new DefaultVisualizer<range_dim, range_dim_cols, RangeField>())
    , name_(nm.empty() ? localizable_function.name() : nm)
    , param_(param)
  {}

  int ncomps() const override final
  {
    return visualizer_.access().ncomps();
  }

  std::string name() const override final
  {
    return name_;
  }

  double evaluate(int comp, const EntityType& en, const DomainType& xx) const override final
  {
    local_function_->bind(en);
    const auto value = local_function_->evaluate(xx, param_);
    return visualizer_.access().evaluate(comp, value);
  }

private:
  mutable std::unique_ptr<LocalFunctionType> local_function_;
  const Common::ConstStorageProvider<VisualizerInterface<range_dim, range_dim_cols, RangeField>> visualizer_;
  const std::string name_;
  const XT::Common::Parameter param_;
}; // class VisualizationAdapter

template <class GridViewType, size_t range_dim, size_t range_dim_cols, class RangeField>
class GradientVisualizationAdapter : public VTKFunction<GridViewType>
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  static_assert(range_dim_cols == 1, "Not implemented!");

public:
  using EntityType = XT::Grid::extract_entity_t<GridViewType>;
  using GridFunctionType = GridFunctionInterface<EntityType, range_dim, range_dim_cols, RangeField>;
  static constexpr size_t d = GridFunctionType::d;

private:
  using LocalFunctionType = typename GridFunctionType::LocalFunctionType;
  using DomainType = typename LocalFunctionType::DomainType;

public:
  GradientVisualizationAdapter(const GridFunctionType& localizable_function,
                               const VisualizerInterface<d, 1, RangeField>& visualizer,
                               const std::string nm = "",
                               const XT::Common::Parameter& param = {})
    : local_function_(localizable_function.local_function())
    , visualizer_(visualizer)
    , name_(nm.empty() ? "grad_ " + localizable_function.name() : nm)
    , param_(param)
  {
    if (range_dim > 1)
      DUNE_THROW(Dune::NotImplemented, "Only implemented for scalar functions by now!");
  }

  GradientVisualizationAdapter(const GridFunctionType& localizable_function,
                               const std::string nm = "",
                               const XT::Common::Parameter& param = {})
    : local_function_(localizable_function.local_function())
    , visualizer_(new DefaultVisualizer<d, 1, RangeField>())
    , name_(nm.empty() ? "grad_" + localizable_function.name() : nm)
    , param_(param)
  {
    if (range_dim > 1)
      DUNE_THROW(Dune::NotImplemented, "Only implemented for scalar functions by now!");
  }

  int ncomps() const override final
  {
    return visualizer_.access().ncomps();
  }

  std::string name() const override final
  {
    return name_;
  }

  double evaluate(int comp, const EntityType& en, const DomainType& xx) const override final
  {
    local_function_->bind(en);
    const auto value = local_function_->jacobian(xx, param_);
    return visualizer_.access().evaluate(comp, value[0]);
  }

private:
  mutable std::unique_ptr<LocalFunctionType> local_function_;
  const Common::ConstStorageProvider<VisualizerInterface<d, 1, RangeField>> visualizer_;
  const std::string name_;
  const XT::Common::Parameter param_;
}; // class GradientVisualizationAdapter


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_VISUALIZATION_HH
