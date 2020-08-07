// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze    (2019)
//   Tobias Leibner (2018 - 2020)

#ifndef DUNE_XT_COMMON_VECTOR_LESS_HH
#define DUNE_XT_COMMON_VECTOR_LESS_HH

#include <dune/xt/common/vector.hh>
#include <dune/xt/common/float_cmp.hh>

namespace Dune::XT::Common {


//! struct to be used as comparison function e.g. in a std::map<FieldVector<...>, ..., VectorLess>
struct VectorLess
{
  template <class FirstVectorType, class SecondVectorType>
  bool operator()(const FirstVectorType& a, const SecondVectorType& b) const
  {
    using V1 = VectorAbstraction<FirstVectorType>;
    using V2 = VectorAbstraction<SecondVectorType>;
    for (size_t dd = 0; dd < a.size(); ++dd) {
      if (V1::get_entry(a, dd) < V2::get_entry(b, dd))
        return true;
      else if (V1::get_entry(a, dd) > V2::get_entry(b, dd))
        return false;
    }
    return false;
  }
};

//! like VectorLess but using FloatCmp::lt instead of "<"
struct VectorFloatLess
{
  template <class FirstVectorType, class SecondVectorType>
  bool operator()(const FirstVectorType& a, const SecondVectorType& b) const
  {
    using V1 = VectorAbstraction<FirstVectorType>;
    using V2 = VectorAbstraction<SecondVectorType>;
    for (size_t dd = 0; dd < a.size(); ++dd) {
      if (XT::Common::FloatCmp::lt(V1::get_entry(a, dd), V2::get_entry(b, dd)))
        return true;
      else if (XT::Common::FloatCmp::gt(V1::get_entry(a, dd), V2::get_entry(b, dd)))
        return false;
    }
    return false;
  }
};


} // namespace Dune::XT::Common

#endif // DUNE_XT_COMMON_VECTOR_LESS_HH
