// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_INTERFACE_HH
#define DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_INTERFACE_HH

#include <dune/xt/common/crtp.hh>

#include <dune/xt/la/container/vector-interface.hh>

namespace Dune {
namespace XT {
namespace LA {


/**
 * \brief Models a collection of length() number of vectors of size dim() each.
 */
template <class TraitsImp>
class VectorArrayInterface : public Common::CRTPInterface<VectorArrayInterface<TraitsImp>, TraitsImp>
{
  using BaseType = Common::CRTPInterface<VectorArrayInterface<TraitsImp>, TraitsImp>;

public:
  using Traits = typename BaseType::Traits;
  using derived_type = typename BaseType::derived_type;

  using VectorType = typename Traits::VectorType;

private:
  static_assert(is_vector<VectorType>::value, "");

public:
  using const_iterator = typename Traits::const_iterator;
  using iterator = typename Traits::iterator;

  size_t dim() const
  {
    CHECK_CRTP(this->as_imp().dim());
    return this->as_imp().dim();
  }

  size_t length() const
  {
    CHECK_CRTP(this->as_imp().length());
    return this->as_imp().length();
  }

  const VectorType& operator[](const size_t ii) const
  {
    CHECK_CRTP(this->as_imp()[ii]);
    return this->as_imp()[ii];
  }

  VectorType& operator[](const size_t ii)
  {
    CHECK_CRTP(this->as_imp()[ii]);
    return this->as_imp()[ii];
  }

  const_iterator begin() const
  {
    CHECK_CRTP(this->as_imp().begin());
    return this->as_imp().begin();
  }

  iterator begin()
  {
    CHECK_CRTP(this->as_imp().begin());
    return this->as_imp().begin();
  }

  const_iterator end() const
  {
    CHECK_CRTP(this->as_imp().end());
    return this->as_imp().end();
  }

  iterator end()
  {
    CHECK_CRTP(this->as_imp().end());
    return this->as_imp().end();
  }

  template <class F>
  void set_all(const F& scalar)
  {
    for (auto& vec : *this)
      vec.set_all(scalar);
  }
}; // class VectorArrayInterface


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_INTERFACE_HH
