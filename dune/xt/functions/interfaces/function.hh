// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018, 2020)
//   René Fritze     (2018 - 2020)
//   Tim Keil        (2018)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_FUNCTION_HH

#include <memory>
#include <map>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/common.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>
#include <dune/grid/io/file/vtk/common.hh>

namespace Dune::XT::Functions {


// forwards, required in operator+-*, includes are below
namespace internal {


template <class, class, class CombinationType>
struct CombinedHelper;


} // namespace internal


template <class, class>
class DifferenceFunction;

template <class, class>
class SumFunction;

template <class, class>
class ProductFunction;

template <class, class>
class FractionFunction;


/**
 * \brief Interface for functions (in the C^\infty sense) which can thus be evaluated in global coordinates.
 *
 *        These functions do not depend on a grid, but only on the dimensions and fields of their domain and range.
 *        See in particular RangeReturnTypeSelector and DerivativeRangeReturnTypeSelector for the interpretation of a
 *        function and its derivatives, and GridFunctionInterface for functions which may have discontinuities between
 *        entities (as in: which are double-valued on intersections).
 *
 *        To turn a function into a function which is localizable w.r.t. a GridView of matching dimension (e.g.,
 *        to visualize it or to use it in a discretization scheme), use make_grid_function to obtain a wrapped version
 *        of this function:
\code
auto function = ...;
auto grid_view = ...;
using E = XT::Grid::extract_entity_t<decltype(grid_view)>;
auto grid_function = make_grid_function<E>(function);
\endcode
 *
 * \sa    RangeReturnTypeSelector
 * \sa    DerivativeRangeReturnTypeSelector
 * \sa    GridFunctionInterface
 * \sa    GridInterface
 * \sa    make_grid_function
 **/
template <size_t domainDim, size_t rangeDim = 1, size_t rangeDimCols = 1, class RangeField = double>
class FunctionInterface : public Common::ParametricInterface
{
  using ThisType = FunctionInterface;

public:
  using DomainFieldType = double;
  static constexpr size_t domain_dim = domainDim;
  using RangeFieldType = RangeField;
  static constexpr size_t range_dim = rangeDim;
  static constexpr size_t range_dim_cols = rangeDimCols;

  using D = DomainFieldType;
  static constexpr size_t d = domain_dim;
  using R = RangeFieldType;
  static constexpr size_t r = range_dim;
  static constexpr size_t rC = range_dim_cols;

  using DomainType = Dune::FieldVector<D, d>;

  using RangeSelector = RangeTypeSelector<R, r, rC>;
  using DerivativeRangeSelector = DerivativeRangeTypeSelector<d, R, r, rC>;

  using DifferenceType = Functions::DifferenceFunction<ThisType, ThisType>;
  using SumType = Functions::SumFunction<ThisType, ThisType>;

  /**
   * \name ``These types are the _standard_ types to be used.''
   * \{
   */

  using RangeReturnType = typename RangeSelector::return_type;
  using DerivativeRangeReturnType = typename DerivativeRangeSelector::return_type;
  using RowDerivativeRangeReturnType = typename DerivativeRangeSelector::row_derivative_return_type;
  using SingleDerivativeRangeReturnType = typename DerivativeRangeSelector::return_single_type;

  /**
   * \}
   * \name ``These types are used for large dimensions.''
   * \{
   */

  using DynamicRangeType = typename RangeSelector::dynamic_type;
  using DynamicRowDerivativeRangeType = typename DerivativeRangeSelector::dynamic_row_derivative_type;
  using DynamicDerivativeRangeType = typename DerivativeRangeSelector::dynamic_type;

  /// \}

  explicit FunctionInterface(const Common::ParameterType& param_type = {})
    : Common::ParametricInterface(param_type)
  {}

  FunctionInterface(const ThisType& other)
    : Common::ParametricInterface(other)
  {}

  FunctionInterface(ThisType&) = default;

  ~FunctionInterface() override = default;

  ThisType& operator=(const ThisType&) = delete;

  ThisType& operator=(ThisType&&) = delete;

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/

  /**
   * \brief Returns a (shallow) copy of the function.
   * actual implementation work is delegated to the private `copy_as_function_impl`
   * combined with hiding `copy_as_function` in derived classes, this allows us the a
   * unique_ptr with correct type at all levels of the polymorphic hierarchy
   *
   * \note This is intended to be cheap, so make sure to share resources (but in a thread-safe way)!
   */
  std::unique_ptr<ThisType> copy_as_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_function_impl());
  }

private:
  virtual ThisType* copy_as_function_impl() const = 0;

public:
  virtual int order(const XT::Common::Parameter& /*param*/ = {}) const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented to define the functionality of the function.''
   * \{
   **/

  virtual RangeReturnType evaluate(const DomainType& /*point_in_global_coordinates*/,
                                   const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This function does not provide evaluations, override the 'evaluate' method!");
  }

  virtual DerivativeRangeReturnType jacobian(const DomainType& /*point_in_global_coordinates*/,
                                             const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This function does not provide a jacobian, override the 'jacobian' method!");
  }

  virtual DerivativeRangeReturnType derivative(const std::array<size_t, d>& /*alpha*/,
                                               const DomainType& /*point_in_global_coordinates*/,
                                               const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented,
               "This function does not provide arbitrary derivatives, override the 'derivative' method!");
  }

  static std::string static_id()
  {
    return "dune.xt.functions.function";
  }

  virtual std::string name() const
  {
    return "dune.xt.functions.function";
  }

  /// \}
  /// \name Operators emulating numeric types.
  /// \{

  Functions::DifferenceFunction<ThisType, ThisType> operator-(const ThisType& other) const
  {
    return Functions::DifferenceFunction<ThisType, ThisType>(*this, other, this->name() + " - " + other.name());
  }

  Functions::SumFunction<ThisType, ThisType> operator+(const ThisType& other) const
  {
    return Functions::SumFunction<ThisType, ThisType>(*this, other, this->name() + " + " + other.name());
  }

  template <class OtherType>
  std::enable_if_t<is_function<OtherType>::value
                       && internal::CombinedHelper<ThisType, OtherType, CombinationType::product>::available,
                   Functions::ProductFunction<ThisType, as_function_interface_t<OtherType>>>
  operator*(const OtherType& other) const
  {
    return Functions::ProductFunction<ThisType, as_function_interface_t<OtherType>>(
        *this, other, "(" + this->name() + ")*(" + other.name() + ")");
  }

  template <class OtherType>
  std::enable_if_t<is_function<OtherType>::value
                       && internal::CombinedHelper<ThisType, OtherType, CombinationType::fraction>::available,
                   Functions::FractionFunction<ThisType, as_function_interface_t<OtherType>>>
  operator/(const OtherType& other) const
  {
    return Functions::FractionFunction<ThisType, as_function_interface_t<OtherType>>(
        *this, other, "(" + this->name() + ")/(" + other.name() + ")");
  }

  /**
   * \}
   * \name ´´These methods are default implemented and should be overridden to improve their performance.''
   * \{
   **/

  virtual R evaluate(const DomainType& point_in_global_coordinates,
                     const size_t row,
                     const size_t col = 0,
                     const Common::Parameter& param = {}) const
  {
    assert_correct_dims(row, col, "evaluate");
    return single_evaluate_helper<R>::call(this->evaluate(point_in_global_coordinates, param), row, col);
  }

  virtual SingleDerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                                   const size_t row,
                                                   const size_t col = 0,
                                                   const Common::Parameter& param = {}) const
  {
    assert_correct_dims(row, col, "jacobian");
    return single_derivative_helper<SingleDerivativeRangeReturnType>::call(
        this->jacobian(point_in_global_coordinates, param), row, col);
  }

  virtual SingleDerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                                     const DomainType& point_in_global_coordinates,
                                                     const size_t row,
                                                     const size_t col = 0,
                                                     const Common::Parameter& param = {}) const
  {
    assert_correct_dims(row, col, "derivative");
    return single_derivative_helper<SingleDerivativeRangeReturnType>::call(
        this->derivative(alpha, point_in_global_coordinates, param), row, col);
  }

  /**
   * \}
   * \name ´´These methods are provided for large dimensions (when RangeReturnType or DerivativeRangeReturnType do not
   *         fit on the stack) and should be overridden to improve their performance.''
   * \{
   **/

  virtual void evaluate(const DomainType& point_in_reference_element,
                        DynamicRangeType& result,
                        const Common::Parameter& param = {}) const
  {
    RangeSelector::ensure_size(result);
    RangeSelector::convert(this->evaluate(point_in_reference_element, param), result);
  }

  virtual void jacobian(const DomainType& point_in_reference_element,
                        DynamicDerivativeRangeType& result,
                        const Common::Parameter& param = {}) const
  {
    DerivativeRangeSelector::ensure_size(result);
    DerivativeRangeSelector::convert(this->jacobian(point_in_reference_element, param), result);
  }

  virtual void derivative(const std::array<size_t, d>& alpha,
                          const DomainType& point_in_reference_element,
                          DynamicDerivativeRangeType& result,
                          const Common::Parameter& param = {}) const
  {
    DerivativeRangeSelector::ensure_size(result);
    DerivativeRangeSelector::convert(this->derivative(alpha, point_in_reference_element, param), result);
  }

protected:
#ifndef DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  static void assert_correct_dims(const size_t row, const size_t col, const std::string& caller)
  {
    if (row >= r || col >= rC)
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "in " << caller << ": this function maps from u:" << d << " -> " << r << "x" << rC << "\n   "
                       << "and you are trying to access component u_" << row << "_" << col << "!");
  }
#else // DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  static void assert_correct_dims(const size_t /*row*/, const size_t /*col*/, const std::string& /*caller*/) {}
#endif

private:
  template <class SingleType, size_t _r = r, size_t _rC = rC, bool anything = true>
  struct single_evaluate_helper
  {
    template <class FullType>
    static SingleType call(const FullType& val, const size_t row, const size_t col)
    {
      return val[row][col];
    }
  }; // struct single_evaluate_helper<r, rC, ...>

  template <class SingleType, size_t _r, bool anything>
  struct single_evaluate_helper<SingleType, _r, 1, anything>
  {
    template <class FullType>
    static SingleType call(const FullType& val, const size_t row, const size_t /*col*/)
    {
      return val[row];
    }
  }; // struct single_evaluate_helper<r, 1, ...>

  template <class SingleType, size_t _r = r, size_t _rC = rC, bool anything = true>
  struct single_derivative_helper
  {
    template <class FullType>
    static SingleType call(const FullType& val, const size_t row, const size_t col)
    {
      SingleType ret;
      for (size_t dd = 0; dd < d; ++dd)
        ret[dd] = val[row][col][dd];
      return ret;
    }
  }; // struct single_derivative_helper<r, rC, ...>

  template <class SingleType, size_t _r, bool anything>
  struct single_derivative_helper<SingleType, _r, 1, anything>
  {
    template <class FullType>
    static SingleType call(const FullType& val, const size_t row, const size_t /*col*/)
    {
      return val[row];
    }
  }; // struct single_derivative_helper<r, 1, ...>
}; // class FunctionInterface


} // namespace Dune::XT::Functions


#include <dune/xt/functions/base/combined-functions.hh>

#endif // DUNE_XT_FUNCTIONS_INTERFACES_FUNCTION_HH
