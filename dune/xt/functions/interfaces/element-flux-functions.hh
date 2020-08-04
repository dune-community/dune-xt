// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2019)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017, 2019 - 2020)

#ifndef DUNE_XT_FUNCTIONS_INTERFACES_ELEMENT_FLUX_FUNCTIONS_HH
#define DUNE_XT_FUNCTIONS_INTERFACES_ELEMENT_FLUX_FUNCTIONS_HH

#include <array>
#include <sstream>
#include <vector>
#include <type_traits>

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/bound-object.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/exceptions.hh>
#include <dune/xt/functions/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \brief Interface for a set of globalvalued functions, which are evaluated with two variables, one variable in
 *element-local coordinates and one grid-independent variable. \sa    RangeTypeSelector \sa    JacobianRangeTypeSelector
 **/
template <class Element,
          size_t stateDim = 1,
          size_t rangeDim = 1,
          size_t rangeDimCols = 1,
          class RangeField = double,
          class StateField = double>
class ElementFluxFunctionSetInterface
  : public Common::ParametricInterface
  , public XT::Grid::ElementBoundObject<Element>
{
  static_assert(XT::Grid::is_entity<Element>::value, "");
  using ThisType = ElementFluxFunctionSetInterface;

public:
  using typename XT::Grid::ElementBoundObject<Element>::ElementType;
  using DomainFieldType = double;
  static const constexpr size_t domain_dim = ElementType::dimension;
  static const constexpr size_t state_dim = stateDim;
  using RangeFieldType = RangeField;
  static const constexpr size_t range_dim = rangeDim;
  static const constexpr size_t range_dim_cols = rangeDimCols;

  using E = ElementType;
  using D = DomainFieldType;
  static const constexpr size_t d = domain_dim;
  static const constexpr size_t s = state_dim;
  using R = RangeFieldType;
  using S = StateField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;

  using DomainType = Dune::FieldVector<D, d>;
  using StateType = Dune::FieldVector<S, s>;
  using DynamicStateType = Dune::DynamicVector<S>;

  using RangeSelector = RangeTypeSelector<R, r, rC>;
  using JacobianRangeSelector = DerivativeRangeTypeSelector<s, R, r, rC>;

  /**
   * \name ``These types are the _standard_ types to be used.''
   * \{
   */

  using RangeType = typename RangeSelector::type;
  using JacobianRangeType = typename JacobianRangeSelector::type;
  using SingleJacobianRangeType = typename JacobianRangeSelector::single_type;

  /**
   * \}
   * \name ``These types are used for large dimensions.''
   * \{
   */

  using DynamicRangeType = typename RangeSelector::dynamic_type;
  using DynamicJacobianRangeType = typename JacobianRangeSelector::dynamic_type;

  /// \}


  ElementFluxFunctionSetInterface(const XT::Common::ParameterType& param_type = {})
    : Common::ParametricInterface(param_type)
  {}

  ElementFluxFunctionSetInterface(const ThisType& other) = default;

  ElementFluxFunctionSetInterface(ThisType&& source) = default;

  virtual ~ElementFluxFunctionSetInterface() = default;

  /**
   * \name ´´These methods have to be implemented.''
   * \{
   **/

  /**
   * The number of functions in this set.
   *
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual size_t size(const Common::Parameter& /*param*/ = {}) const = 0;

  /**
   * The maximum number of functions in this set over all elements. The following has to hold true:
\code
auto local_function_set = ...;
size_t max_set_size = 0;
for (auto&& element : grid_view)
  max_set_size = std::max(max_set_size, local_function_set.bind(element).size(param));
assert(max_set_size <= local_function_set.max_size());
\endcode
   **/
  virtual size_t max_size(const Common::Parameter& /*param*/ = {}) const = 0;

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual int order(const Common::Parameter& /*param*/ = {}) const = 0;

  /**
   * \}
   * \name ´´These methods should be implemented to define the functionality of the set of functions.''
   * \{
   **/

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const DomainType& /*point_in_reference_element*/,
                        const StateType& /*u*/,
                        std::vector<RangeType>& /*result*/,
                        const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented,
               "This set of element functions does not provide evaluations, override the 'evaluate' method!");
  }

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void jacobian(const DomainType& /*point_in_reference_element*/,
                        const StateType& /*u*/,
                        std::vector<JacobianRangeType>& /*result*/,
                        const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented,
               "This set of element functions does not provide jacobian, override the 'jacobian' method!");
  }

  /**
   * \{
   * \name ´´These methods are provided for convenience and should not be used within library code.''
   * \{
   **/

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual std::vector<RangeType> evaluate_set(const DomainType& point_in_reference_element,
                                              const StateType& u,
                                              const Common::Parameter& param = {}) const
  {
    std::vector<RangeType> result(this->size(param));
    this->evaluate(point_in_reference_element, u, result, param);
    return result;
  }

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual std::vector<JacobianRangeType> jacobian_of_set(const DomainType& point_in_reference_element,
                                                         const StateType& u,
                                                         const Common::Parameter& param = {}) const
  {
    std::vector<JacobianRangeType> result(this->size(param));
    this->jacobian(point_in_reference_element, u, result, param);
    return result;
  }

  /**
   * \{
   * \name ´´These methods are used to access individual range dimensions and should be overridden to improve their
   *         performance.''
   * \{
   **/

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const DomainType& point_in_reference_element,
                        const StateType& u,
                        std::vector<R>& result,
                        const size_t row,
                        const size_t col = 0,
                        const Common::Parameter& param = {}) const
  {
    assert_correct_dims(row, col, "evaluate");
    const auto tmp_values = this->evaluate_set(point_in_reference_element, u, param);
    if (result.size() < tmp_values.size())
      result.resize(tmp_values.size());
    if constexpr (rC == 1) {
      for (size_t ii = 0; ii < tmp_values.size(); ++ii)
        result[ii] = tmp_values[ii][row];
    } else {
      for (size_t ii = 0; ii < tmp_values.size(); ++ii)
        result[ii] = tmp_values[ii][row][col];
    }
  }

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void jacobian(const DomainType& point_in_reference_element,
                        const StateType& u,
                        std::vector<SingleJacobianRangeType>& result,
                        const size_t row,
                        const size_t col = 0,
                        const Common::Parameter& param = {}) const
  {
    assert_correct_dims(row, col, "jacobian");
    const auto tmp_values = this->jacobian_of_set(point_in_reference_element, u, param);
    if (result.size() < tmp_values.size())
      result.resize(tmp_values.size());
    if constexpr (rC == 1) {
      for (size_t ii = 0; ii < tmp_values.size(); ++ii)
        result[ii] = tmp_values[ii][row];
    } else {
      for (size_t ii = 0; ii < tmp_values.size(); ++ii)
        for (size_t dd = 0; dd < d; ++dd)
          result[ii][dd] = tmp_values[ii][row][col][dd];
    }
  }

  /**
   * \{
   * \name ´´These methods are provided for large dimensions (when RangeType or JacobianRangeType do not fit on the
   *         stack) and should be overridden to improve their performance.''
   * \{
   **/

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void evaluate(const DomainType& point_in_reference_element,
                        const StateType& u,
                        std::vector<DynamicRangeType>& result,
                        const Common::Parameter& param = {}) const
  {
    // prepare result
    const auto sz = this->size();
    if (result.size() < sz)
      result.resize(sz);
    for (size_t ii = 0; ii < sz; ++ii)
      RangeSelector::ensure_size(result[ii]);
    // call actual evaluate
    auto tmp_result = std::make_unique<std::vector<RangeType>>(sz);
    this->evaluate(point_in_reference_element, u, *tmp_result, param);
    // convert
    for (size_t ii = 0; ii < sz; ++ii)
      RangeSelector::convert((*tmp_result)[ii], result[ii]);
  } // ... evaluate(...)

  /**
   * \note Will throw Exceptions::not_bound_to_an_element_yet error if not bound yet!
   **/
  virtual void jacobian(const DomainType& point_in_reference_element,
                        const StateType& u,
                        std::vector<DynamicJacobianRangeType>& result,
                        const Common::Parameter& param = {}) const
  {
    // prepare result
    const auto sz = this->size();
    if (result.size() < sz)
      result.resize(sz);
    for (size_t ii = 0; ii < sz; ++ii)
      JacobianRangeSelector::ensure_size(result[ii]);
    // call actual jacobian
    auto tmp_result = std::make_unique<std::vector<JacobianRangeType>>(sz);
    this->jacobian(point_in_reference_element, u, *tmp_result, param);
    // convert
    for (size_t ii = 0; ii < sz; ++ii)
      JacobianRangeSelector::convert((*tmp_result)[ii], result[ii]);
  } // ... jacobian(...)

  /**
   * \}
   **/
protected:
#ifndef DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  void assert_inside_reference_element(const DomainType& point_in_reference_element) const
  {
    if (!ReferenceElements<D, d>::general(this->element().type()).checkInside(point_in_reference_element)) {
      std::stringstream error_message;
      error_message << "Given point_in_reference_element is not inside the current elements reference element!\n\n";
      XT::Grid::print_entity(this->element(), XT::Common::Typename<E>::value(), error_message, "   ");
      error_message << "\n   point_in_reference_element = " << point_in_reference_element << std::endl;
      DUNE_THROW(XT::Functions::Exceptions::wrong_input_given, error_message.str());
    }
  }
#else // DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  void assert_inside_reference_element(const DomainType& /*point_in_reference_element*/) const {}
#endif

#ifndef DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  static void assert_correct_dims(const size_t row, const size_t col, const std::string& caller)
  {
    if (row >= r || col >= rC)
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "in " << caller << ": this function maps from " << d << " -> " << r << "x" << rC << "\n   "
                       << "and you are trying to access the entry in row " << row << " and column " << col
                       << " of its range!");
  }
#else // DUNE_XT_FUNCTIONS_DISABLE_CHECKS
  static void assert_correct_dims(const size_t /*row*/, const size_t /*col*/, const std::string& /*caller*/) {}
#endif
}; // class ElementFluxFunctionSetInterface


/**
 *  \brief  Interface for a globalvalued function, which can be evaluated locally on one Element.
 */
template <class Element,
          size_t state_dim = 1,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double>
class ElementFluxFunctionInterface
  : public ElementFluxFunctionSetInterface<Element, state_dim, range_dim, range_dim_cols, RangeField>
{
  using BaseType = ElementFluxFunctionSetInterface<Element, state_dim, range_dim, range_dim_cols, RangeField>;
  using ThisType = ElementFluxFunctionInterface;

public:
  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using BaseType::s;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicStateType;
  using typename BaseType::ElementType;
  using typename BaseType::R;
  using typename BaseType::S;
  using typename BaseType::StateType;

  using typename BaseType::JacobianRangeSelector;
  using typename BaseType::RangeSelector;

  /**
   * \name ``These types are the _standard_ types used by the LocalFunctionSetInterface.''
   * \{
   */

  using typename BaseType::JacobianRangeType;
  using typename BaseType::RangeType;
  using typename BaseType::SingleJacobianRangeType;

  /**
   * \}
   * \name ``These types are used by the LocalFunctionSetInterface for large dimensions.''
   * \{
   */

  using typename BaseType::DynamicJacobianRangeType;
  using typename BaseType::DynamicRangeType;

  /**
   * \name ``These types are the _standard_ types to be used for LocalFunctionInterface.''
   * \{
   */

  using RangeReturnType = typename RangeSelector::return_type;
  using JacobianRangeReturnType = typename JacobianRangeSelector::return_type;
  using SingleJacobianRangeReturnType = typename JacobianRangeSelector::return_single_type;

  /// \}

  ElementFluxFunctionInterface(const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
  {}

  ElementFluxFunctionInterface(const ThisType& other) = default;
  ElementFluxFunctionInterface(ThisType&& source) = default;

  virtual ~ElementFluxFunctionInterface() = default;

  ThisType& operator=(const ThisType& other) = default;
  ThisType& operator=(ThisType&& source) = default;

  using BaseType::evaluate;
  using BaseType::jacobian;

  /**
   * \}
   * \name ´´These methods should be implemented to define the functionality of the function.''
   * \{
   **/

  virtual RangeReturnType evaluate(const DomainType& /*point_in_reference_element*/,
                                   const StateType& /*u*/,
                                   const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This local function does not provide evaluations, override the 'evaluate' method!");
  }

  // This is df/du, so the jacobian with respect to the state variable u. We do not yet provide derivatives in x (feel
  // free to implement them if you need x-derivatives).
  virtual JacobianRangeReturnType jacobian(const DomainType& /*point_in_reference_element*/,
                                           const StateType& /*u*/,
                                           const Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "This local function does not provide a jacobian, override the 'jacobian' method!");
  }

  /**
   * \}
   * \name ´´These methods are used to access individual range dimensions and should be overridden to improve their
   *         performance.''
   * \{
   **/

  virtual R evaluate(const DomainType& point_in_reference_element,
                     const StateType& u,
                     const size_t row,
                     const size_t col = 0,
                     const Common::Parameter& param = {}) const
  {
    this->assert_correct_dims(row, col, "evaluate");

    const auto val = this->evaluate(point_in_reference_element, u, param);
    if constexpr (rC == 1) {
      return val[row];
    } else {
      return val[row][col];
    }
  }

  virtual SingleJacobianRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                                 const StateType& u,
                                                 const size_t row,
                                                 const size_t col = 0,
                                                 const Common::Parameter& param = {}) const
  {
    this->assert_correct_dims(row, col, "jacobian");
    const auto val = this->jacobian(point_in_reference_element, u, param);
    if constexpr (rC == 1) {
      return val[row];
    } else {
      SingleJacobianRangeReturnType ret;
      for (size_t dd = 0; dd < d; ++dd)
        ret[dd] = val[row][col][dd];
      return ret;
    }
  }

  /**
   * \}
   * \name ´´These methods are provided for large dimensions (when RangeType or JacobianRangeType do not fit on the
   *         stack) and should be overridden to improve their performance.''
   * \{
   **/

  virtual void evaluate(const DomainType& point_in_reference_element,
                        const StateType& u,
                        DynamicRangeType& result,
                        const Common::Parameter& param = {}) const
  {
    // prepare result
    RangeSelector::ensure_size(result);
    // call actual evaluate
    auto tmp_result = std::make_unique<std::vector<RangeType>>(1);
    this->evaluate(point_in_reference_element, u, *tmp_result, param);
    // convert
    RangeSelector::convert((*tmp_result)[0], result);
  } // ... evaluate(...)

  virtual void jacobian(const DomainType& point_in_reference_element,
                        const StateType& u,
                        DynamicJacobianRangeType& result,
                        const Common::Parameter& param = {}) const
  {
    // prepare result
    JacobianRangeSelector::ensure_size(result);
    // call actual jacobian
    auto tmp_result = std::make_unique<std::vector<JacobianRangeType>>(1);
    this->jacobian(point_in_reference_element, u, *tmp_result, param);
    // convert
    JacobianRangeSelector::convert((*tmp_result)[0], result);
  } // ... jacobian(...)

  /**
   * \{
   * \name ´´These methods are required by ElementFluxFunctionSetInterface and are provided by this interface.''
   * \{
   **/

  virtual size_t size(const Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  virtual size_t max_size(const Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  virtual void evaluate(const DomainType& point_in_reference_element,
                        const StateType& u,
                        std::vector<RangeType>& result,
                        const Common::Parameter& param = {}) const override
  {
    if (result.size() < 1)
      result.resize(1);
    result[0] = this->evaluate(point_in_reference_element, u, param);
  }

  virtual void jacobian(const DomainType& point_in_reference_element,
                        const StateType& u,
                        std::vector<JacobianRangeType>& result,
                        const Common::Parameter& param = {}) const override
  {
    if (result.size() < 1)
      result.resize(1);
    result[0] = this->jacobian(point_in_reference_element, u, param);
  }
}; // class ElementFluxFunctionInterface


} // namespace Functions
} // namespace XT
} // namespace Dune

#include <dune/xt/functions/base/combined-element-functions.hh>

#endif // DUNE_XT_FUNCTIONS_INTERFACES_ELEMENT_FLUX_FUNCTIONS_HH
