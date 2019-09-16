// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)

#ifndef DUNE_XT_COMMON_FVECTOR_PBH
#define DUNE_XT_COMMON_FVECTOR_PBH

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/cast.h>

#include <dune/xt/common/fvector.hh>

NAMESPACE_BEGIN(pybind11)
NAMESPACE_BEGIN(detail)


template <class FieldVectorImp,
          bool is_number = (Dune::XT::Common::is_arithmetic<typename FieldVectorImp::value_type>::value
                            || Dune::XT::Common::is_complex<typename FieldVectorImp::value_type>::value)>
struct FieldVector_type_caster
{
  using type = FieldVectorImp;
  typedef typename type::value_type K;
  static const int SZ = type::dimension;
  using value_conv = make_caster<K>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<sequence>(src))
      return false;
    auto s = reinterpret_borrow<sequence>(src);
    if (s.size() != SZ)
      return false;
    value_conv conv;
    value *= K(0.0);
    size_t ii = 0;
    for (auto it : s) {
      if (ii >= SZ)
        return false;
      if (!conv.load(it, convert))
        return false;
      value[ii++] = cast_op<K>(conv);
    }
    return true;
  } // ... load(...)

  static handle cast(const type& src, return_value_policy policy, handle parent)
  {
    list l(SZ);
    for (size_t ii = 0; ii < src.size(); ++ii) {
      object val = reinterpret_steal<object>(value_conv::cast(src[ii], policy, parent));
      if (!val)
        return handle();
      PyList_SET_ITEM(l.ptr(), ii, val.release().ptr()); // steals a reference
    }
    return l.release();
  } // ... cast(...)

  PYBIND11_TYPE_CASTER(type, _("List[") + value_conv::name + _("[") + _<SZ>() + _("]]"));
}; // struct FieldVector_type_caster

/**
 * This specialization is needed because we also want to have std::string in FieldVector.
 * For Dune::XT::Common::FieldMatrix we already have the correct default value, so no string version of value *= K(0.0)
 * is required
 **/
template <class FieldVectorImp>
struct FieldVector_type_caster<FieldVectorImp, false>
{
  using type = FieldVectorImp;
  typedef typename type::value_type K;
  static const int SZ = type::dimension;
  using value_conv = make_caster<K>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<sequence>(src))
      return false;
    auto s = reinterpret_borrow<sequence>(src);
    if (s.size() != SZ)
      return false;
    value_conv conv;
    size_t ii = 0;
    for (auto it : s) {
      if (ii >= SZ)
        return false;
      if (!conv.load(it, convert))
        return false;
      value[ii++] = cast_op<K>(conv);
    }
    return true;
  } // ... load(...)

  static handle cast(const type& src, return_value_policy policy, handle parent)
  {
    list l(SZ);
    for (size_t ii = 0; ii < src.size(); ++ii) {
      object val = reinterpret_steal<object>(value_conv::cast(src[ii], policy, parent));
      if (!val)
        return handle();
      PyList_SET_ITEM(l.ptr(), ii, val.release().ptr()); // steals a reference
    }
    return l.release();
  } // ... cast(...)

  PYBIND11_TYPE_CASTER(type, _("List[") + value_conv::name + _("[") + _<SZ>() + _("]]"));
}; // struct FieldVector_type_caster


template <class K, int SIZE>
struct type_caster<Dune::FieldVector<K, SIZE>> : public FieldVector_type_caster<Dune::FieldVector<K, SIZE>>
{};

template <class K, int SIZE>
struct type_caster<Dune::XT::Common::FieldVector<K, SIZE>>
  : public FieldVector_type_caster<Dune::XT::Common::FieldVector<K, SIZE>>
{};


NAMESPACE_END(detail)
NAMESPACE_END(pybind11)

#endif // DUNE_XT_COMMON_FVECTOR_PBH
