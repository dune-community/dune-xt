// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Tobias Leibner  (2017)

#ifndef DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH
#define DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace XT {
namespace Functions {


// forwards
template <class E, class D, size_t d, class R, size_t r, size_t rC>
class LocalfunctionSetInterface;

template <class E, class D, size_t d, class R, size_t r, size_t rC>
class LocalfunctionInterface;

template <class E, class D, size_t d, class R, size_t r, size_t rC>
class LocalizableFunctionInterface;

template <class E, class D, size_t d, class U_, size_t s_, class R, size_t r, size_t rC>
class LocalizableFluxFunctionInterface;


namespace internal {


template <class Tt>
struct is_localfunction_set_helper
{
  DXTC_has_typedef_initialize_once(EntityType);
  DXTC_has_typedef_initialize_once(DomainFieldType);
  DXTC_has_typedef_initialize_once(RangeFieldType);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);
  DXTC_has_static_member_initialize_once(dimRangeCols);
  static const bool is_candidate =
      DXTC_has_typedef(EntityType)<Tt>::value && DXTC_has_typedef(DomainFieldType)<Tt>::value
      && DXTC_has_typedef(RangeFieldType)<Tt>::value && DXTC_has_static_member(dimDomain)<Tt>::value
      && DXTC_has_static_member(dimRange)<Tt>::value && DXTC_has_static_member(dimRangeCols)<Tt>::value;
}; // struct is_localfunction_set_helper

template <class Tt>
struct is_localizable_flux_function_helper
{
  DXTC_has_typedef_initialize_once(EntityType);
  DXTC_has_typedef_initialize_once(DomainFieldType);
  DXTC_has_typedef_initialize_once(StateType);
  DXTC_has_typedef_initialize_once(RangeFieldType);
  DXTC_has_static_member_initialize_once(dimDomain);
  DXTC_has_static_member_initialize_once(dimRange);
  DXTC_has_static_member_initialize_once(dimRangeCols);
  DXTC_has_static_member_initialize_once(stateDerivativeOrder);
  static const bool is_candidate =
      DXTC_has_typedef(EntityType)<Tt>::value && DXTC_has_typedef(DomainFieldType)<Tt>::value
      && DXTC_has_typedef(StateType)<Tt>::value && DXTC_has_typedef(RangeFieldType)<Tt>::value
      && DXTC_has_static_member(dimDomain)<Tt>::value && DXTC_has_static_member(dimRange)<Tt>::value
      && DXTC_has_static_member(dimRangeCols)<Tt>::value && DXTC_has_static_member(stateDerivativeOrder)<Tt>::value;
}; // struct is_localfunction_set_helper


} // namespace internal


template <class T, bool candidate = internal::is_localfunction_set_helper<T>::is_candidate>
struct is_localfunction_set : public std::is_base_of<LocalfunctionSetInterface<typename T::EntityType,
                                                                               typename T::DomainFieldType,
                                                                               T::dimDomain,
                                                                               typename T::RangeFieldType,
                                                                               T::dimRange,
                                                                               T::dimRangeCols>,
                                                     T>
{
};

template <class T>
struct is_localfunction_set<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_localfunction_set_helper<T>::is_candidate>
struct is_localfunction : public std::is_base_of<LocalfunctionInterface<typename T::EntityType,
                                                                        typename T::DomainFieldType,
                                                                        T::dimDomain,
                                                                        typename T::RangeFieldType,
                                                                        T::dimRange,
                                                                        T::dimRangeCols>,
                                                 T>
{
};

template <class T>
struct is_localfunction<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_localfunction_set_helper<T>::is_candidate>
struct is_localizable_function : public std::is_base_of<LocalizableFunctionInterface<typename T::EntityType,
                                                                                     typename T::DomainFieldType,
                                                                                     T::dimDomain,
                                                                                     typename T::RangeFieldType,
                                                                                     T::dimRange,
                                                                                     T::dimRangeCols>,
                                                        T>
{
};

template <class T>
struct is_localizable_function<T, false> : public std::false_type
{
};


template <class T, bool candidate = internal::is_localizable_flux_function_helper<T>::is_candidate>
struct is_localizable_flux_function
    : public std::is_base_of<LocalizableFluxFunctionInterface<typename T::EntityType,
                                                              typename T::DomainFieldType,
                                                              T::dimDomain,
                                                              typename T::StateType,
                                                              T::stateDerivativeOrder,
                                                              typename T::RangeFieldType,
                                                              T::dimRange,
                                                              T::dimRangeCols>,
                             T>
{
};

template <class T>
struct is_localizable_flux_function<T, false> : public std::false_type
{
};


//! Utility to generate a complete Function Type from an existing one and a template
template <class FunctionImp, template <class, class, size_t, class, size_t, size_t> class OutTemplate>
struct FunctionTypeGenerator
{
  typedef OutTemplate<typename FunctionImp::EntityType,
                      typename FunctionImp::DomainFieldType,
                      FunctionImp::dimDomain,
                      typename FunctionImp::RangeFieldType,
                      FunctionImp::dimRange,
                      FunctionImp::dimRangeCols>
      type;
};

//! Utility for uniform handling of RangeType for functions with dimRangeCols == 1
//! and functions with dimRangeCols > 1.
template <size_t dimRange, size_t dimRangeCols>
struct RangeTypeConverter
{
  typedef Dune::FieldMatrix<double, dimRange, dimRangeCols> UnifiedRangeType;
  typedef UnifiedRangeType RangeType;

  static UnifiedRangeType convert(const UnifiedRangeType& in)
  {
    return in;
  }

  static RangeType convert_back(const UnifiedRangeType& in)
  {
    return in;
  }
};

template <size_t dimRange>
struct RangeTypeConverter<dimRange, 1>
{
  typedef Dune::FieldMatrix<double, dimRange, 1> UnifiedRangeType;
  typedef Dune::FieldVector<double, dimRange> RangeType;

  static UnifiedRangeType convert(const RangeType& in)
  {
    UnifiedRangeType out;
    for (size_t ii = 0; ii < dimRange; ++ii)
      out[ii][0] = in[ii];
    return out;
  }

  static RangeType convert_back(const UnifiedRangeType& in)
  {
    RangeType out;
    for (size_t ii = 0; ii < dimRange; ++ii)
      out[ii] = in[ii][0];
    return out;
  }
};

//! Utility for uniform handling of JacobianRangeType for functions with dimRangeCols == 1
//! and functions with dimRangeCols > 1.
template <size_t dimDomain, size_t dimRange, size_t dimRangeCols>
struct JacobianRangeTypeConverter
{
  typedef typename Dune::FieldVector<Dune::FieldMatrix<double, dimRange, dimDomain>, dimRangeCols> JacobianRangeType;
  typedef JacobianRangeType UnifiedJacobianRangeType;

  static UnifiedJacobianRangeType convert(const JacobianRangeType& ret)
  {
    return ret;
  }

  static JacobianRangeType convert_back(const UnifiedJacobianRangeType& ret)
  {
    return ret;
  }
};

template <size_t dimDomain, size_t dimRange>
struct JacobianRangeTypeConverter<dimDomain, dimRange, 1>
{
  typedef typename Dune::FieldMatrix<double, dimRange, dimDomain> JacobianRangeType;
  typedef typename Dune::FieldVector<JacobianRangeType, 1> UnifiedJacobianRangeType;

  static UnifiedJacobianRangeType convert(const JacobianRangeType& in)
  {
    return UnifiedJacobianRangeType(in);
  }

  static JacobianRangeType convert_back(const UnifiedJacobianRangeType& in)
  {
    return in[0];
  }
};


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_TYPE_TRAITS_HH
