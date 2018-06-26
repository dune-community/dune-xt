// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH
#define DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH

#include <vector>

#include "interface.hh"

namespace Dune {
namespace XT {
namespace LA {


// forward, needed for the Traits
template <class Vector>
class ListVectorArray;


namespace internal {


template <class Vector>
class ListVectorArrayTraits
{
  static_assert(is_vector<Vector>::value, "");

public:
  using derived_type = ListVectorArray<Vector>;
  using VectorType = Vector;
  using const_iterator = typename std::vector<Vector>::const_iterator;
  using iterator = typename std::vector<Vector>::iterator;
}; // class ListVectorArrayTraits


} // namespace internal


/**
 * \brief Implementation of VectorArrayInterface as an array of vectors derived from XT::LA::VectorInterface
 */
template <class Vector>
class ListVectorArray : public VectorArrayInterface<internal::ListVectorArrayTraits<Vector>>
{
  using BaseType = VectorArrayInterface<internal::ListVectorArrayTraits<Vector>>;

public:
  using Traits = typename BaseType::Traits;
  using derived_type = typename BaseType::derived_type;

  using typename BaseType::VectorType;
  using typename BaseType::const_iterator;
  using typename BaseType::iterator;

  ListVectorArray(const size_t dm, const size_t lngth = 0, const size_t resrv = 0)
    : dim_(dm)
    , array_(lngth)
  {
    for (auto& v : array_)
      v = VectorType(dm, 0.);
    array_.reserve(resrv);
  }

  size_t dim() const
  {
    return dim_;
  }

  size_t length() const
  {
    return array_.size();
  }

  void reserve(const size_t len)
  {
    array_.reserve(len);
  }

  void append(const VectorType& vec)
  {
    array_.emplace_back(vec);
  }

  void append(VectorType&& vec)
  {
    array_.emplace_back(std::move(vec));
  }

  const VectorType& operator[](const size_t ii) const
  {
    return array_[ii];
  }

  VectorType& operator[](const size_t ii)
  {
    return array_[ii];
  }

  const_iterator begin() const
  {
    return array_.begin();
  }

  iterator begin()
  {
    return array_.begin();
  }

  const_iterator end() const
  {
    return array_.end();
  }

  iterator end()
  {
    return array_.end();
  }

private:
  const size_t dim_;
  std::vector<Vector> array_;
}; // class ListVectorArray


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH
