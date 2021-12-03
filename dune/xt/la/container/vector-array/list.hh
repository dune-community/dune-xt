// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018 - 2020)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH
#define DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH

#include <vector>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class Vector>
class ListVectorArray
{
  static_assert(is_vector<Vector>::value);

public:
  using ThisType = ListVectorArray;

  using VectorType = Vector;
  using VectorArrayType = std::vector<Common::StorageProvider<VectorType>>;

private:
  class AnnotatedVector
  {
  public:
    AnnotatedVector(VectorArrayType& vectors, std::vector<Common::Parameter>& notes, const size_t ii)
      : vectors_(vectors)
      , notes_(notes)
      , ii_(ii)
    {}

    const VectorType& vector() const
    {
      DUNE_THROW_IF(ii_ >= vectors_.size(),
                    InvalidStateException,
                    "This should not happen: ii_ = " << ii_ << "\n   vectors_.size() = " << vectors_.size());
      return vectors_[ii_].access();
    }

    VectorType& vector()
    {
      DUNE_THROW_IF(ii_ >= vectors_.size(),
                    InvalidStateException,
                    "This should not happen: ii_ = " << ii_ << "\n   vectors_.size() = " << vectors_.size());
      return vectors_[ii_].access();
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
    VectorArrayType& vectors_;
    std::vector<Common::Parameter>& notes_;
    const size_t ii_;
  }; // class AnnotatedVector

  using AnnotatedVectorArrayType = std::vector<AnnotatedVector>;

public:
  using iterator = typename AnnotatedVectorArrayType::iterator;
  using const_iterator = typename AnnotatedVectorArrayType::const_iterator;
  using reverse_iterator = typename AnnotatedVectorArrayType::reverse_iterator;
  using const_reverse_iterator = typename AnnotatedVectorArrayType::const_reverse_iterator;
  using reference = typename AnnotatedVectorArrayType::reference;
  using const_reference = typename AnnotatedVectorArrayType::const_reference;

  ListVectorArray(const size_t dm, const size_t lngth = 0, const size_t resrv = 0)
    : dim_(dm)
    , len_(lngth)
    , vectors_()
    , pairs_()
  {
    vectors_.reserve(resrv);
    notes_.reserve(resrv);
    pairs_.reserve(resrv);
    for (size_t ii = 0; ii < len_; ++ii) {
      vectors_.emplace_back(new VectorType(dim_, 0.));
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

  ListVectorArray(ThisType&& source) noexcept
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

  const VectorArrayType& vectors() const
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

  // careful, makes a copy
  void append(const VectorType& vec, const Common::Parameter& note = {})
  {
    vectors_.emplace_back(new VectorType(vec));
    notes_.emplace_back(note);
    pairs_.emplace_back(vectors_, notes_, len_);
    ++len_;
  }

  void append(VectorType& vec, const Common::Parameter& note = {})
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
  VectorArrayType vectors_;
  std::vector<Common::Parameter> notes_;
  AnnotatedVectorArrayType pairs_;
}; // class ListVectorArray


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_VECTOR_ARRAY_LIST_HH
