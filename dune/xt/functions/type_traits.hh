// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)
//   Tim Keil        (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH
#define DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * Given a function u: R^d -> R^{r x rC}, we interpret it as a r functions u_1, ..., u_r in the sense of
 * u(x) = (u_1(x), ..., u_r(x))T, where each u_s: R^d -> R^rC for 1 \leq s \leq r, and thus ultimately
 *
 *        | (u_1(x))_1, (u_1(x))_2, ..., (u_1(x))_rC |
 *        | (u_2(x))_1,    .        ...      .       |
 * u(x) = |      .         .        ...      .       |
 *        |      .         .        ...      .       |
 *        | (u_r(x))_1, (u_r(x))_2, ..., (u_r(x))_rC |
 *
 * (thus u_s corresponds to one row of u, for s = 1, ..., r).
 * Here, (u_s(x))_i can be interpreted as the evaluation of a scalar function u_s_i: R^d -> R
 * for 1 \leq i \leq rC, and we can thus identify
 *
 * u = ( u_s_i )_{s = 1, ..., r | i = 1, ..., rC }
 *
 * This interpretation of u is modeled by RangeTypeSelector, with a special case if rC = 1, to simplify manners.
 **/
template <class R, size_t r, size_t rC>
struct RangeTypeSelector
{
  using type = FieldMatrix<R, r, rC>;
  using return_type = XT::Common::FieldMatrix<R, r, rC>;
  using dynamic_type = DynamicMatrix<R>;

  static void ensure_size(dynamic_type& arg)
  {
    if (arg.rows() < r || arg.cols() < rC)
      arg.resize(r, rC);
  }

  static void convert(const type& in, dynamic_type& out)
  {
    assert(out.rows() >= r);
    assert(out.cols() >= rC);
    for (size_t ii = 0; ii < r; ++ii)
      for (size_t jj = 0; jj < rC; ++jj)
        out[ii][jj] = in[ii][jj];
  }
}; // struct RangeTypeSelector<...>

template <class R, size_t r>
struct RangeTypeSelector<R, r, 1>
{
  using type = FieldVector<R, r>;
  using return_type = XT::Common::FieldVector<R, r>;
  using dynamic_type = DynamicVector<R>;

  static void ensure_size(dynamic_type& arg)
  {
    if (arg.size() < r)
      arg.resize(r);
  }

  static void convert(const type& in, dynamic_type& out)
  {
    assert(out.size() >= r);
    for (size_t ii = 0; ii < r; ++ii)
      out[ii] = in[ii];
  }
}; // struct RangeTypeSelector<..., 1>


/**
 * \sa RangeTypeSelector
 *
 *     Since one is often interested in the derivatives of the individual functions u_s (see RangeTypeSelector),
 *     it makes sense to consider the jacobian u: R^d -> R^{r x {rC x d}} as a vector of the jacobians of the
 *     u_s (for s = 1, ..., r) in the sense of
 *
 *     jacobian u = (jacobian u_1, ..., jacobian u_r)
 *
 *     where
 *
 *                       |  (d u_s_1 / d x_1)(x),  (d u_s_1 / d x_2)(x), ...,  (d u_s_1 / d x_d)(x) |
 *                       |  (d u_s_2 / d x_1)(x),           .            ...            .           |
 *     jacobian u_s(x) = |           .                      .            ...            .           |
 *                       |           .                      .            ...            .           |
 *                       | (d u_s_rC / d x_1)(x), (d u_s_rC / d x_2)(x), ..., (d u_s_rC / d x_d)(x) |
 *
 *     for s = 1, ..., r.
 *     Again, (d u_s_i / d x_j)(x) can be interpreted as the evaluation of a scalar function
 *     (jacobian u)_s_i_j: R^d -> R, and we can thus identify
 *
 *     (jacobian u)_s = ( d u_s_i / d x_j )_{i = 1, ..., rC | j = 1, ..., d}
 *
 *     This interpretation of the jacobian of u is modeled by DerivativeRangeTypeSelector, with a special case if
 *     rC = 1, to simplify manners. To access the jacobian of u_s, simply use
 *     \code
 *     u.jacobian(x)[s];
 *     \endcode
 *
 *     In addition, one is sometimes interested in evaluating the jacobian of the scalar function u_s_i (see
 *     RangeTypeSelector), jacobian u_s_i: R^d -> R^d. This is modeled by the single_type.
 **/
template <size_t d, class R, size_t r, size_t rC>
struct DerivativeRangeTypeSelector
{
  using single_type = XT::Common::FieldVector<R, d>;
  using return_single_type = FieldVector<R, d>;
  using dynamic_single_type = DynamicVector<R>;

  using type = FieldVector<FieldMatrix<R, rC, d>, r>;
  using return_type = XT::Common::FieldVector<XT::Common::FieldMatrix<R, rC, d>, r>;
  using dynamic_type = DynamicVector<DynamicMatrix<R>>;

  static void ensure_size(dynamic_type& arg)
  {
    if (arg.size() < r)
      arg.resize(r);
    for (size_t ii = 0; ii < r; ++ii)
      if (arg[ii].rows() < rC || arg[ii].cols() < d)
        arg[ii].resize(rC, d);
  }

  static void convert(const type& in, dynamic_type& out)
  {
    assert(out.size() >= r);
    for (size_t ii = 0; ii < r; ++ii) {
      assert(out[ii].rows() >= rC);
      assert(out[ii].cols() >= d);
      for (size_t jj = 0; jj < rC; ++jj)
        for (size_t kk = 0; kk < d; ++kk)
          out[ii][jj][kk] = in[ii][jj][kk];
    }
  } // ... convert(...)
}; // struct DerivativeRangeTypeSelector<...>

template <size_t d, class R, size_t r>
struct DerivativeRangeTypeSelector<d, R, r, 1>
{
  using single_type = FieldVector<R, d>;
  using return_single_type = XT::Common::FieldVector<R, d>;
  using dynamic_single_type = DynamicVector<R>;

  using type = FieldMatrix<R, r, d>;
  using return_type = XT::Common::FieldMatrix<R, r, d>;
  using dynamic_type = DynamicMatrix<R>;

  static void ensure_size(dynamic_type& arg)
  {
    if (arg.rows() < r || arg.cols() < d)
      arg.resize(r, d);
  }

  static void convert(const type& in, dynamic_type& out)
  {
    assert(out.rows() >= r);
    assert(out.cols() >= d);
    for (size_t ii = 0; ii < r; ++ii)
      for (size_t jj = 0; jj < d; ++jj)
        out[ii][jj] = in[ii][jj];
  }
}; // struct DerivativeRangeTypeSelector<..., 1>


// forwards
template <class E, size_t r, size_t rC, class R>
class ElementFunctionSetInterface;

template <class E, size_t r, size_t rC, class R>
class ElementFunctionInterface;

template <class E, size_t r, size_t rC, class R>
class GridFunctionInterface;

// smooth

// template <class E, class D, size_t d, class U_, size_t s_, class R, size_t r, size_t rC>
// class LocalizableFluxFunctionInterface;


namespace internal {


template <class Tt>
struct is_grid_function_helper
{
  DXTC_has_typedef_initialize_once(E);
  DXTC_has_typedef_initialize_once(R);
  DXTC_has_static_member_initialize_once(r);
  DXTC_has_static_member_initialize_once(rC);
  static const bool is_candidate = DXTC_has_typedef(E)<Tt>::value && DXTC_has_typedef(R)<Tt>::value
                                   && DXTC_has_static_member(r)<Tt>::value && DXTC_has_static_member(rC)<Tt>::value;
}; // struct is_grid_function_helper


} // namespace internal


template <class T, bool is_candidate = internal::is_grid_function_helper<T>::is_candidate>
struct is_grid_function;

template <class T>
struct is_grid_function<T, false> : public std::false_type
{};

template <class T>
struct is_grid_function<T, true> : std::is_base_of<GridFunctionInterface<typename T::E, T::r, T::rC, typename T::R>, T>
{};


// template <class T, bool candidate = internal::is_localfunction_set_helper<T>::is_candidate>
// struct is_localfunction_set : public std::is_base_of<LocalfunctionSetInterface<typename T::EntityType,
//                                                                               typename T::DomainFieldType,
//                                                                               T::dimDomain,
//                                                                               typename T::RangeFieldType,
//                                                                               T::dimRange,
//                                                                               T::dimRangeCols>,
//                                                     T>
//{
//};

// template <class T>
// struct is_localfunction_set<T, false> : public std::false_type
//{
//};


// template <class T, bool candidate = internal::is_localfunction_set_helper<T>::is_candidate>
// struct is_localfunction : public std::is_base_of<LocalfunctionInterface<typename T::EntityType,
//                                                                        typename T::DomainFieldType,
//                                                                        T::dimDomain,
//                                                                        typename T::RangeFieldType,
//                                                                        T::dimRange,
//                                                                        T::dimRangeCols>,
//                                                 T>
//{
//};

// template <class T>
// struct is_localfunction<T, false> : public std::false_type
//{
//};


// template <class T, bool candidate = internal::is_localfunction_set_helper<T>::is_candidate>
// struct is_grid_function : public std::is_base_of<GridFunctionInterface<typename T::EntityType,
//                                                                                     typename T::DomainFieldType,
//                                                                                     T::dimDomain,
//                                                                                     typename T::RangeFieldType,
//                                                                                     T::dimRange,
//                                                                                     T::dimRangeCols>,
//                                                        T>
//{
//};

// template <class T>
// struct is_grid_function<T, false> : public std::false_type
//{
//};


// template <class T, bool candidate = internal::is_localizable_flux_function_helper<T>::is_candidate>
// struct is_localizable_flux_function
//    : public std::is_base_of<LocalizableFluxFunctionInterface<typename T::EntityType,
//                                                              typename T::DomainFieldType,
//                                                              T::dimDomain,
//                                                              typename T::StateType,
//                                                              T::stateDerivativeOrder,
//                                                              typename T::RangeFieldType,
//                                                              T::dimRange,
//                                                              T::dimRangeCols>,
//                             T>
//{
//};

// template <class T>
// struct is_localizable_flux_function<T, false> : public std::false_type
//{
//};


////! Utility to generate a complete Function Type from an existing one and a template
// template <class FunctionImp, template <class, class, size_t, class, size_t, size_t> class OutTemplate>
// struct FunctionTypeGenerator
//{
//  typedef OutTemplate<typename FunctionImp::EntityType,
//                      typename FunctionImp::DomainFieldType,
//                      FunctionImp::dimDomain,
//                      typename FunctionImp::RangeFieldType,
//                      FunctionImp::dimRange,
//                      FunctionImp::dimRangeCols>
//      type;
//};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH
