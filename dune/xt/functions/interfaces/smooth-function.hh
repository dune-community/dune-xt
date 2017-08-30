// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_SMOOTH_FUNCTION_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_SMOOTH_FUNCTION_HH

#include <dune/common/fvector.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief Interface for smooth functions (in the C^\infty sense) which can thus be evaluated in global coordinates.
 *
 *        These functions do not depend on a grid, but only on the dimensions and fields of their domain and range.
 *        See in particular RangeTypeSelector and DerivativeRangeTypeSelector for the interpretation of a function and
 *        its derivatives, and LocalizableFunctionInterface for functions which may have discontinuities between entites
 *        (as in: which are double-valued on intersections).
 *
 * \sa    RangeTypeSelector
 * \sa    DerivativeRangeTypeSelector
 * \sa    LocalizableFunctionInterface
 **/
template <size_t domainDim, size_t rangeDim = 1, size_t rangeDimCols = 1, class RangeFieldImp = double>
class SmoothFunctionInterface : public Common::ParametricInterface
{
public:
  using DomainFieldType = double;
  static const constexpr size_t dimDomain = domainDim;
  using RangeFieldType = RangeFieldImp;
  static const constexpr size_t dimRange = rangeDim;
  static const constexpr size_t dimRangeCols = rangeDimCols;

  using D = DomainFieldType;
  static const constexpr size_t d = dimDomain;
  using R = RangeFieldType;
  static const constexpr size_t r = dimRange;
  static const constexpr size_t rC = dimRangeCols;

  using DomainType = FieldVector<D, d>;
  using RangeType = typename RangeTypeSelector<R, r, rC>::type;
  using DerivativeRangeType = typename DerivativeRangeTypeSelector<d, R, r, rC>::type;
  using SingleDerivativeRangeType = typename DerivativeRangeTypeSelector<d, R, r, rC>::single_type;

  virtual ~SmoothFunctionInterface() = default;

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/

  virtual int order(const XT::Common::Parameter& /*mu*/ = {}) const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented to define the functionality of the function.''
   * \{
   **/

  virtual RangeType evaluate(const DomainType& /*xx*/, const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This smooth function does not provide evaluations, override the 'evaluate' method!");
  }

  virtual DerivativeRangeType jacobian(const DomainType& /*xx*/, const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This smooth function does not provide a jacobian, override the 'jacobian' method!");
  }

  virtual DerivativeRangeType derivative(const std::array<size_t, d>& /*alpha*/,
                                         const DomainType& /*xx*/,
                                         const Common::Parameter& /*mu*/ = {}) const
  {
    DUNE_THROW(NotImplemented,
               "This smooth function does not provide arbitrary derivatives, override the 'derivative' method!");
  }

  virtual std::string type() const
  {
    return "smooth_function";
  }

  virtual std::string name() const
  {
    return "smooth_function";
  }

  /**
   * \}
   * \name ´´These methods can be overridden to improve their performance.''
   * \{
   **/

  virtual R evaluate(const DomainType& xx, const size_t row, const size_t col, const Common::Parameter& mu = {}) const
  {
    ensure_correct_dims(row, col, "evaluate");
    return single_helper<R>::call(this->evaluate(xx, mu), row, col);
  }

  virtual SingleDerivativeRangeType
  jacobian(const DomainType& xx, const size_t row, const size_t col, const Common::Parameter& mu = {}) const
  {
    ensure_correct_dims(row, col, "jacobian");
    return single_helper<SingleDerivativeRangeType>::call(this->jacobian(xx, mu), row, col);
  }

  virtual SingleDerivativeRangeType derivative(const std::array<size_t, d>& alpha,
                                               const DomainType& xx,
                                               const size_t row,
                                               const size_t col,
                                               const Common::Parameter& mu = {}) const
  {
    ensure_correct_dims(row, col, "derivative");
    return single_helper<SingleDerivativeRangeType>::call(this->derivative(alpha, xx, mu), row, col);
  }

protected:
#ifndef DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  static void ensure_correct_dims(const size_t row, const size_t col, const std::string& caller)
  {
    if (row >= r || col >= rC)
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "in " << caller << ": this function maps from u:" << d << " -> " << r << "x" << rC << "\n   "
                       << "and you are trying to access component u_"
                       << row
                       << "_"
                       << col
                       << "!");
  }
#else // DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  static void ensure_correct_dims(const size_t /*row*/, const size_t /*col*/, const std::string& /*caller*/)
  {
  }
#endif

private:
  template <class SingleType, size_t _r = r, size_t _rC = rC, bool anything = true>
  struct single_helper
  {
    template <class FullType>
    static SingleType call(const FullType& val, const size_t row, const size_t col)
    {
      SingleType ret = 0;
      for (size_t dd = 0; dd < d; ++dd)
        ret[dd] = val[row][col][dd];
      return ret;
    }
  }; // struct single_helper<r, rC, ...>

  template <class SingleType, size_t _r, bool anything>
  struct single_helper<SingleType, _r, 1, anything>
  {
    template <class FullType>
    static SingleType call(const FullType& val, const size_t row, const size_t /*col*/)
    {
      return val[row];
    }
  }; // struct single_helper<r, 1, ...>
}; // class SmoothFunctionInterface


} // namespace Functions
} // namespace XT
} // namespace Dune

#include <dune/xt/functions/base/smooth-localizable-function.hh>

#endif // DUNE_XT_FUNCTIONS_INTERFACES_SMOOTH_FUNCTION_HH
