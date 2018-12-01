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

#include <dune/xt/common/parameter.hh>
#include <dune/xt/la/exceptions.hh>

//#include "interface.hh"

namespace Dune {
namespace XT {
namespace LA {


// forward, needed for the Traits
template <class Vector>
class ListVectorArray;


#if 0
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
#endif // 0


/**
 * \brief Implementation of VectorArrayInterface as an array of vectors
 */
template <class Vector>
class ListVectorArray /*: public VectorArrayInterface<internal::ListVectorArrayTraits<Vector>>*/
{
  //  using BaseType = VectorArrayInterface<internal::ListVectorArrayTraits<Vector>>;
  using ThisType = ListVectorArray<Vector>;

public:
  //  using Traits = typename BaseType::Traits;
  //  using derived_type = typename BaseType::derived_type;

  //  using typename BaseType::VectorType;
  //  using typename BaseType::const_iterator;
  //  using typename BaseType::iterator;

  using VectorType = Vector;

private:
  class AnnotatedVector
  {
  public:
    AnnotatedVector(std::vector<Vector>& vectors, std::vector<Common::Parameter>& notes, const size_t ii)
      : vectors_(vectors)
      , notes_(notes)
      , ii_(ii)
    {}

    const VectorType& vector() const
    {
      DUNE_THROW_IF(ii_ >= vectors_.size(),
                    InvalidStateException,
                    "This should not happen: ii_ = " << ii_ << "\n   vectors_.size() = " << vectors_.size());
      return vectors_[ii_];
    }

    VectorType& vector()
    {
      DUNE_THROW_IF(ii_ >= vectors_.size(),
                    InvalidStateException,
                    "This should not happen: ii_ = " << ii_ << "\n   vectors_.size() = " << vectors_.size());
      return vectors_[ii_];
    }

    const Common::Parameter& note() const
    {
      DUNE_THROW_IF(ii_ >= vectors_.size(),
                    InvalidStateException,
                    "This should not happen: ii_ = " << ii_ << "\n   notes_.size() = " << notes_.size());
      return notes_[ii_];
    }

    Common::Parameter& note()
    {
      DUNE_THROW_IF(ii_ >= vectors_.size(),
                    InvalidStateException,
                    "This should not happen: ii_ = " << ii_ << "\n   notes_.size() = " << notes_.size());
      return notes_[ii_];
    }

  private:
    std::vector<Vector>& vectors_;
    std::vector<Common::Parameter>& notes_;
    const size_t ii_;
  }; // class AnnotatedVector

public:
  using iterator = typename std::vector<AnnotatedVector>::iterator;
  using const_iterator = typename std::vector<AnnotatedVector>::const_iterator;
  using reverse_iterator = typename std::vector<AnnotatedVector>::reverse_iterator;
  using const_reverse_iterator = typename std::vector<AnnotatedVector>::const_reverse_iterator;
  using reference = typename std::vector<AnnotatedVector>::reference;
  using const_reference = typename std::vector<AnnotatedVector>::const_reference;

  ListVectorArray(const size_t dm, const size_t lngth = 0, const size_t resrv = 0)
    : dim_(dm)
    , len_(lngth)
    , vectors_()
    , notes_()
    , pairs_()
  {
    vectors_.reserve(resrv);
    notes_.reserve(resrv);
    pairs_.reserve(resrv);
    for (size_t ii = 0; ii < len_; ++ii) {
      vectors_.emplace_back(dim_, 0);
      notes_.emplace_back();
      pairs_.emplace_back(vectors_, notes_, ii);
    }
  } // ListVectorArray(...)

  ListVectorArray(const ThisType& other)
    : dim_(other.dim_)
    , len_(other.len_)
    , vectors_(other.vectors_)
    , notes_(other.notes_)
    , pairs_()
  {
    pairs_.reserve(len_);
    for (size_t ii = 0; ii < len_; ++ii)
      pairs_.emplace_back(vectors_, notes_, ii);
  } // ListVectorArray(...)

  ListVectorArray(ThisType&& source)
    : dim_(source.dim_)
    , len_(source.len_)
    , vectors_(std::move(source.vectors_))
    , notes_(std::move(source.notes_))
    , pairs_()
  {
    pairs_.reserve(len_);
    for (size_t ii = 0; ii < len_; ++ii)
      pairs_.emplace_back(vectors_, notes_, ii);
  } // ListVectorArray(...)

  size_t dim() const
  {
    return dim_;
  }

  size_t length() const
  {
    return len_;
  }

  const std::vector<VectorType>& vectors() const
  {
    return vectors_;
  }

  const std::vector<Common::Parameter>& notes() const
  {
    return notes_;
  }

  void reserve(const size_t len)
  {
    vectors_.reserve(len);
    notes_.reserve(len);
    pairs_.reserve(len);
  }

  void append(const VectorType& vec, const Common::Parameter& note = {})
  {
    vectors_.emplace_back(vec);
    notes_.emplace_back(note);
    pairs_.emplace_back(vectors_, notes_, len_);
    ++len_;
  }

  void append(VectorType&& vec, const Common::Parameter& note = {})
  {
    vectors_.emplace_back(std::move(vec));
    notes_.emplace_back(note);
    pairs_.emplace_back(vectors_, notes_, len_);
    ++len_;
  }

  const AnnotatedVector& operator[](const size_t ii) const
  {
    DUNE_THROW_IF(ii >= len_, Common::Exceptions::index_out_of_range, "ii = " << ii << "\n   len_ = " << len_);
    return pairs_[ii];
  }

  AnnotatedVector& operator[](const size_t ii)
  {
    DUNE_THROW_IF(ii >= len_, Common::Exceptions::index_out_of_range, "ii = " << ii << "\n   len_ = " << len_);
    return pairs_[ii];
  }

  reference front()
  {
    return pairs_.front();
  }

  const_reference front() const
  {
    return pairs_.front();
  }

  reference back()
  {
    return pairs_.back();
  }

  const_reference back() const
  {
    return pairs_.back();
  }

  iterator begin()
  {
    return pairs_.begin();
  }

  const_iterator begin() const
  {
    return pairs_.begin();
  }

  const_iterator cbegin() const
  {
    return pairs_.cbegin();
  }

  iterator end()
  {
    return pairs_.end();
  }

  const_iterator end() const
  {
    return pairs_.end();
  }

  const_iterator cend() const
  {
    return pairs_.cend();
  }

  reverse_iterator rbegin()
  {
    return pairs_.rbegin();
  }

  const_reverse_iterator rbegin() const
  {
    return pairs_.rbegin();
  }

  const_reverse_iterator crbegin() const
  {
    return pairs_.crbegin();
  }

  reverse_iterator rend()
  {
    return pairs_.rend();
  }

  const_reverse_iterator rend() const
  {
    return pairs_.end();
  }

  const_reverse_iterator crend() const
  {
    return pairs_.crend();
  }

private:
  const size_t dim_;
  size_t len_;
  std::vector<Vector> vectors_;
  std::vector<Common::Parameter> notes_;
  std::vector<AnnotatedVector> pairs_;
}; // class ListVectorArray


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH
