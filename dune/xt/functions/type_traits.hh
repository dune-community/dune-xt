// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH
#define DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH

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
  using type = XT::Common::FieldMatrix<R, r, rC>;
};

template <class R, size_t r>
struct RangeTypeSelector<R, r, 1>
{
  using type = XT::Common::FieldVector<R, r>;
};


/**
 * \sa RangeTypeSelector
 *
 *     Since one is often interested in the derivatives of the individual functions u_s, it makes sense to consider the
 *     jacobian u: R^d -> R^{r x {d x rC}}, where
 *
 *                         |  (d u_s_1 / d x_1)(x),  (d u_s_1 / d x_2)(x), ...,  (d u_s_1 / d x_d)(x) |
 *                         |  (d u_s_2 / d x_1)(x),           .             ...           .           |
 *     (jacobian u)_s(x) = |           .                      .             ...           .           |
 *                         |           .                      .             ...           .           |
 *                         | (d u_s_rC / d x_1)(x), (d u_s_rC / d x_2)(x), ..., (d u_s_rC / d x_d)(x) |
 *
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
  using type = XT::Common::FieldVector<Dune::XT::Common::FieldMatrix<R, rC, d>, r>;
};

template <size_t d, class R, size_t r>
struct DerivativeRangeTypeSelector<d, R, r, 1>
{
  using single_type = XT::Common::FieldVector<R, d>;
  using type = XT::Common::FieldMatrix<R, r, d>;
};


//// forwards
// template <class E, class D, size_t d, class R, size_t r, size_t rC>
// class LocalfunctionSetInterface;

// template <class E, class D, size_t d, class R, size_t r, size_t rC>
// class LocalfunctionInterface;

// template <class E, class D, size_t d, class R, size_t r, size_t rC>
// class LocalizableFunctionInterface;

// template <class E, class D, size_t d, class U_, size_t s_, class R, size_t r, size_t rC>
// class LocalizableFluxFunctionInterface;


// namespace internal {


// template <class Tt>
// struct is_localfunction_set_helper
//{
//  DXTC_has_typedef_initialize_once(EntityType);
//  DXTC_has_typedef_initialize_once(DomainFieldType);
//  DXTC_has_typedef_initialize_once(RangeFieldType);
//  DXTC_has_static_member_initialize_once(dimDomain);
//  DXTC_has_static_member_initialize_once(dimRange);
//  DXTC_has_static_member_initialize_once(dimRangeCols);
//  static const bool is_candidate =
//      DXTC_has_typedef(EntityType)<Tt>::value && DXTC_has_typedef(DomainFieldType)<Tt>::value
//      && DXTC_has_typedef(RangeFieldType)<Tt>::value && DXTC_has_static_member(dimDomain)<Tt>::value
//      && DXTC_has_static_member(dimRange)<Tt>::value && DXTC_has_static_member(dimRangeCols)<Tt>::value;
//}; // struct is_localfunction_set_helper

// template <class Tt>
// struct is_localizable_flux_function_helper
//{
//  DXTC_has_typedef_initialize_once(EntityType);
//  DXTC_has_typedef_initialize_once(DomainFieldType);
//  DXTC_has_typedef_initialize_once(StateType);
//  DXTC_has_typedef_initialize_once(RangeFieldType);
//  DXTC_has_static_member_initialize_once(dimDomain);
//  DXTC_has_static_member_initialize_once(dimRange);
//  DXTC_has_static_member_initialize_once(dimRangeCols);
//  DXTC_has_static_member_initialize_once(stateDerivativeOrder);
//  static const bool is_candidate =
//      DXTC_has_typedef(EntityType)<Tt>::value && DXTC_has_typedef(DomainFieldType)<Tt>::value
//      && DXTC_has_typedef(StateType)<Tt>::value && DXTC_has_typedef(RangeFieldType)<Tt>::value
//      && DXTC_has_static_member(dimDomain)<Tt>::value && DXTC_has_static_member(dimRange)<Tt>::value
//      && DXTC_has_static_member(dimRangeCols)<Tt>::value && DXTC_has_static_member(stateDerivativeOrder)<Tt>::value;
//}; // struct is_localfunction_set_helper


//} // namespace internal


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
// struct is_localizable_function : public std::is_base_of<LocalizableFunctionInterface<typename T::EntityType,
//                                                                                     typename T::DomainFieldType,
//                                                                                     T::dimDomain,
//                                                                                     typename T::RangeFieldType,
//                                                                                     T::dimRange,
//                                                                                     T::dimRangeCols>,
//                                                        T>
//{
//};

// template <class T>
// struct is_localizable_function<T, false> : public std::false_type
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
