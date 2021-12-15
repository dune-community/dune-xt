// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   René Fritze     (2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2020 - 2021)

#ifndef DUNE_XT_COMMON_FMATRIX_PBH
#define DUNE_XT_COMMON_FMATRIX_PBH

#include <dune/pybindxi/complex.h>
#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/cast.h>

#include <dune/xt/common/fmatrix.hh>

#include "fvector.hh"

PYBIND11_NAMESPACE_BEGIN(pybind11)
PYBIND11_NAMESPACE_BEGIN(detail)


template <class FieldMatrixImp>
struct FieldMatrix_type_caster
{
  using type = FieldMatrixImp;
  using K = typename type::value_type;
  using row_type = typename type::row_type;
  static constexpr int ROWS = type::rows;
  static constexpr int COLS = type::cols;
  using value_conv = make_caster<K>;
  using row_conv = make_caster<row_type>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<sequence>(src))
      return false;
    auto s = reinterpret_borrow<sequence>(src);
    if (s.size() != ROWS)
      return false;
    row_conv conv;
    value *= K(0.0);
    size_t ii = 0;
    for (auto it : s) {
      if (ii >= ROWS)
        return false;
      if (!conv.load(it, convert))
        return false;
      value[ii++] = cast_op<row_type>(conv);
    }
    return true;
  } // ... load(...)

  static handle cast(const type& src, return_value_policy policy, handle parent)
  {
    list l(ROWS);
    for (size_t ii = 0; ii < src.size(); ++ii) {
      auto val = reinterpret_steal<object>(row_conv::cast(src[ii], policy, parent));
      if (!val)
        return handle();
      PyList_SET_ITEM(l.ptr(), ii, val.release().ptr()); // steals a reference
    }
    return l.release();
    return list(0).release();
  } // ... cast(...)

  PYBIND11_TYPE_CASTER(type,
                       _("List[List[") + value_conv::name + _("[") + _<COLS>() + _("]]") + _("[") + _<ROWS>()
                           + _("]]"));
}; // struct FieldMatrix_type_caster


template <class K, int N, int M>
struct type_caster<Dune::FieldMatrix<K, N, M>> : public FieldMatrix_type_caster<Dune::FieldMatrix<K, N, M>>
{};

/**
 * Dune::XT::Common::FieldMatrix is already able to handle std::string as we already provide the correct default value
 * (see type_traits::suitable_default). Thus, we do not need to initalize with value *= K(0.0)
 *
 * Attention: If this is not enough you need to specialize for std::string. (see fvector.hh)
 */
template <class K, int N, int M>
struct type_caster<Dune::XT::Common::FieldMatrix<K, N, M>>
  : public FieldMatrix_type_caster<Dune::XT::Common::FieldMatrix<K, N, M>>
{
  using type = Dune::XT::Common::FieldMatrix<K, N, M>;
  using row_type = typename type::row_type;
  static constexpr int ROWS = type::rows;
  static constexpr int COLS = type::cols;
  using value_conv = make_caster<K>;
  using row_conv = make_caster<row_type>;

  bool load(handle src, bool convert)
  {
    if (!isinstance<sequence>(src))
      return false;
    auto s = reinterpret_borrow<sequence>(src);
    if (s.size() != ROWS)
      return false;
    row_conv conv;
    size_t ii = 0;
    for (auto it : s) {
      if (ii >= ROWS)
        return false;
      if (!conv.load(it, convert))
        return false;
      value[ii++] = cast_op<row_type>(conv);
    }
    return true;
  } // ... load(...)

  PYBIND11_TYPE_CASTER(type,
                       _("List[List[") + value_conv::name + _("[") + _<COLS>() + _("]]") + _("[") + _<ROWS>()
                           + _("]]"));
};


PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(pybind11)

#endif // DUNE_XT_COMMON_FMATRIX_PBH
