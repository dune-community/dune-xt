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

#include <map>

#include <dune/common/fvector.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


// forward, required in SmoothFunctionInterface::as_localizable
template <class E, size_t r, size_t rC, class R>
class SmoothFunctionAsLocalizableWrapper;


/**
 * \brief Interface for smooth functions (in the C^\infty sense) which can thus be evaluated in global coordinates.
 *
 *        These functions do not depend on a grid, but only on the dimensions and fields of their domain and range.
 *        See in particular RangeTypeSelector and DerivativeRangeTypeSelector for the interpretation of a function and
 *        its derivatives, and LocalizableFunctionInterface for functions which may have discontinuities between entites
 *        (as in: which are double-valued on intersections).
 *
 *        To turn a smooth function into a function which is localizable w.r.t. a GridView of matching dimension (e.g.,
 *        to visualize it or to use it in a discretization scheme), use as_localizable to obtain a const reference to
 *        a wrapped version of this function:
\code
auto smooth_function = ...;
auto grid_view = ...;
using E = XT::Grid::extract_entity_t<decltype(grid_view)>;
const auto& localizable_function = smooth_function.template as_localizable<E>();
\endcode
 *
 * \sa    RangeTypeSelector
 * \sa    DerivativeRangeTypeSelector
 * \sa    LocalizableFunctionInterface
 **/
template <size_t domainDim, size_t rangeDim = 1, size_t rangeDimCols = 1, class RangeFieldImp = double>
class SmoothFunctionInterface : public Common::ParametricInterface
{
  using ThisType = SmoothFunctionInterface<domainDim, rangeDim, rangeDimCols, RangeFieldImp>;

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

  virtual int order(const XT::Common::Parameter& /*param*/ = {}) const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented to define the functionality of the function.''
   * \{
   **/

  virtual RangeType evaluate(const DomainType& /*point_in_global_coordinates*/,
                             const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This smooth function does not provide evaluations, override the 'evaluate' method!");
  }

  virtual DerivativeRangeType jacobian(const DomainType& /*point_in_global_coordinates*/,
                                       const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This smooth function does not provide a jacobian, override the 'jacobian' method!");
  }

  virtual DerivativeRangeType derivative(const std::array<size_t, d>& /*alpha*/,
                                         const DomainType& /*point_in_global_coordinates*/,
                                         const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented,
               "This smooth function does not provide arbitrary derivatives, override the 'derivative' method!");
  }

  virtual std::string type() const
  {
    return "dune.xt.functions.smooth";
  }

  static std::string static_id()
  {
    return "dune.xt.functions.smooth";
  }

  virtual std::string name() const
  {
    return "dune.xt.functions.smooth";
  }

  /**
   * \}
   * \name ´´These methods are default implemented and should be overridden to improve their performance.''
   * \{
   **/

  virtual R evaluate(const DomainType& point_in_global_coordinates,
                     const size_t row,
                     const size_t col,
                     const Common::Parameter& param = {}) const
  {
    ensure_correct_dims(row, col, "evaluate");
    return single_evaluate_helper<R>::call(this->evaluate(point_in_global_coordinates, param), row, col);
  }

  virtual SingleDerivativeRangeType jacobian(const DomainType& point_in_global_coordinates,
                                             const size_t row,
                                             const size_t col,
                                             const Common::Parameter& param = {}) const
  {
    ensure_correct_dims(row, col, "jacobian");
    return single_derivative_helper<SingleDerivativeRangeType>::call(
        this->jacobian(point_in_global_coordinates, param), row, col);
  }

  virtual SingleDerivativeRangeType derivative(const std::array<size_t, d>& alpha,
                                               const DomainType& point_in_global_coordinates,
                                               const size_t row,
                                               const size_t col,
                                               const Common::Parameter& param = {}) const
  {
    ensure_correct_dims(row, col, "derivative");
    return single_derivative_helper<SingleDerivativeRangeType>::call(
        this->derivative(alpha, point_in_global_coordinates, param), row, col);
  }

  /**
   * \}
   * \name ´´These methods are provided.''
   * \{
   **/

  /**
   * \note This function kepps a map of all wrappers in a local static map, to avoid temporaries.
   * \todo Check if this implementation is thread safe!
   */
  template <class E>
  const typename std::enable_if<XT::Grid::is_entity<E>::value && E::dimension == d,
                                SmoothFunctionAsLocalizableWrapper<E, r, rC, R>>::type&
  as_localizable() const
  {
    static std::map<const ThisType*, std::unique_ptr<SmoothFunctionAsLocalizableWrapper<E, r, rC, R>>> wrappers;
    if (wrappers.find(this) == wrappers.end())
      wrappers[this] = std::make_unique<SmoothFunctionAsLocalizableWrapper<E, r, rC, R>>(*this);
    return *(wrappers[this]);
  }

  /**
   * \copydoc LocalizableFunctionInterface::visualize
   */
  template <class GridLayerType>
  typename std::enable_if<Grid::is_layer<GridLayerType>::value, void>::type
  visualize(const GridLayerType& grid_layer,
            const std::string path,
            const bool subsampling = true,
            const VTK::OutputType vtk_output_type = VTK::appendedraw) const
  {
    this->as_localizable<XT::Grid::extract_entity_t<GridLayerType>>().visualize(
        grid_layer, path, subsampling, vtk_output_type);
  }

  /**
   * \}
   **/

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
}; // class SmoothFunctionInterface


} // namespace Functions
} // namespace XT
} // namespace Dune

#include <dune/xt/functions/base/smooth-localizable-function.hh>

#endif // DUNE_XT_FUNCTIONS_INTERFACES_SMOOTH_FUNCTION_HH
