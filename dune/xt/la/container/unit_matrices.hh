// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Barbara Verfürth (2015)
//   Felix Schindler  (2013 - 2017)
//   Rene Milk        (2014 - 2016)
//   Tobias Leibner   (2014, 2016 - 2017)

#ifndef DUNE_XT_LA_CONTAINER_UNIT_MATRICES_HH
#define DUNE_XT_LA_CONTAINER_UNIT_MATRICES_HH

#include <memory>

#include <dune/common/fmatrix.hh>

#include <dune/xt/common/memory.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/pattern.hh>

namespace Dune {
namespace XT {
namespace LA {


template <class MatrixType, class Enable = void>
struct UnitMatrix;

template <class MatrixType>
struct UnitMatrix<MatrixType,
                  std::enable_if_t<std::is_base_of<Dune::FieldMatrix<typename MatrixType::value_type,
                                                                     MatrixType::rows,
                                                                     MatrixType::rows>,
                                                   MatrixType>::value,
                                   void>>
{
  static std::unique_ptr<MatrixType> get(const size_t /*size*/ = MatrixType::rows, const size_t /*num_mutexes*/ = 1)
  {
    auto ret = XT::Common::make_unique<MatrixType>(0.);
    for (size_t ii = 0; ii < MatrixType::rows; ++ii)
      (*ret)[ii][ii] = 1.;
    return ret;
  }
};

template <class MatrixType>
struct UnitMatrix<MatrixType,
                  std::enable_if_t<std::is_base_of<Dune::DynamicMatrix<typename MatrixType::value_type>,
                                                   MatrixType>::value,
                                   void>>
{
  static std::unique_ptr<MatrixType> get(const size_t size, const size_t /*num_mutexes*/ = 1)
  {
    auto ret = XT::Common::make_unique<MatrixType>(size, size, 0.);
    for (size_t ii = 0; ii < size; ++ii)
      (*ret)[ii][ii] = 1.;
    return ret;
  }
};

template <class MatrixVectorType>
struct UnitMatrix<MatrixVectorType,
                  std::enable_if_t<std::is_base_of<Dune::FieldVector<Dune::FieldMatrix<
                                                                         typename MatrixVectorType::value_type::
                                                                             value_type,
                                                                         MatrixVectorType::value_type::rows,
                                                                         MatrixVectorType::value_type::rows>,
                                                                     MatrixVectorType::dimension>,
                                                   MatrixVectorType>::value,
                                   void>>
{
  static std::unique_ptr<MatrixVectorType> get(const size_t size, const size_t num_mutexes = 1)
  {
    auto ret = XT::Common::make_unique<MatrixVectorType>(0.);
    for (size_t jj = 0; jj < size; ++jj)
      for (size_t ii = 0; ii < MatrixVectorType::value_type::rows; ++ii)
        (*ret)[jj][ii][ii] = 1.;
    return ret;
  }
};

template <class MatrixType>
struct UnitMatrix<MatrixType, std::enable_if_t<XT::LA::is_matrix<MatrixType>::value, void>>
{
  static std::unique_ptr<MatrixType> get(const size_t size, const size_t num_mutexes = 1)
  {
    XT::LA::SparsityPatternDefault pattern(size);
    for (size_t ii = 0; ii < size; ++ii)
      pattern.insert(ii, ii);
    auto ret = XT::Common::make_unique<MatrixType>(size, size, pattern, num_mutexes);
    for (size_t ii = 0; ii < size; ++ii)
      ret->set_entry(ii, ii, 1.);
    return ret;
  }
};


} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_CONTAINER_UNIT_MATRICES_HH
