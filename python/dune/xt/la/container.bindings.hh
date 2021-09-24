// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2020)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_LA_CONTAINER_BINDINGS_HH
#define DUNE_XT_LA_CONTAINER_BINDINGS_HH

#include <dune/xt/la/container.hh>


// this is used in other headers
using COMMON_DENSE_VECTOR = Dune::XT::LA::CommonDenseVector<double>;
using COMMON_DENSE_MATRIX = Dune::XT::LA::CommonDenseMatrix<double>;
using COMMON_SPARSE_MATRIX = Dune::XT::LA::CommonSparseMatrix<double>;
#if HAVE_EIGEN
using EIGEN_DENSE_VECTOR = Dune::XT::LA::EigenDenseVector<double>;
using EIGEN_DENSE_MATRIX = Dune::XT::LA::EigenDenseMatrix<double>;
using EIGEN_SPARSE_MATRIX = Dune::XT::LA::EigenRowMajorSparseMatrix<double>;
#endif
using ISTL_DENSE_VECTOR = Dune::XT::LA::IstlDenseVector<double>;
using ISTL_SPARSE_MATRIX = Dune::XT::LA::IstlRowMajorSparseMatrix<double>;

namespace Dune {
namespace XT {
namespace LA {
namespace bindings {


template <Backends bb>
struct backend_name
{
  static_assert(AlwaysFalse<typename Container<double, bb>::VectorType1>::value,
                "Please add a specialization for this backend!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct backend_name<Backends::common_dense>
{
  static std::string value()
  {
    return "common_dense";
  }
};

template <>
struct backend_name<Backends::common_sparse>
{
  static std::string value()
  {
    return "common_sparse";
  }
};

template <>
struct backend_name<Backends::istl_dense>
{
  static std::string value()
  {
    return "istl_dense";
  }
};

template <>
struct backend_name<Backends::istl_sparse>
{
  static std::string value()
  {
    return "istl_sparse";
  }
};

template <>
struct backend_name<Backends::eigen_dense>
{
  static std::string value()
  {
    return "eigen_dense";
  }
};

template <>
struct backend_name<Backends::eigen_sparse>
{
  static std::string value()
  {
    return "eigen_sparse";
  }
};


template <class T>
struct container_name
{
  static_assert(AlwaysFalse<T>::value, "Please add a specialization for this container!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct container_name<CommonDenseVector<double>>
{
  static std::string value()
  {
    return "common_vector";
  }
};

template <>
struct container_name<CommonDenseVector<size_t>>
{
  static std::string value()
  {
    return "common_vector_size_t";
  }
};

template <>
struct container_name<CommonDenseMatrix<double>>
{
  static std::string value()
  {
    return "common_dense_matrix";
  }
};

template <>
struct container_name<CommonSparseMatrix<double>>
{
  static std::string value()
  {
    return "common_sparse_matrix";
  }
};

#if HAVE_EIGEN

template <>
struct container_name<EigenDenseVector<double>>
{
  static std::string value()
  {
    return "eigen_vector";
  }
};

template <>
struct container_name<EigenDenseMatrix<double>>
{
  static std::string value()
  {
    return "eigen_dense_matrix";
  }
};

template <>
struct container_name<EigenRowMajorSparseMatrix<double>>
{
  static std::string value()
  {
    return "eigen_sparse_matrix";
  }
};

#endif // HAVE_EIGEN

template <>
struct container_name<IstlDenseVector<double>>
{
  static std::string value()
  {
    return "istl_vector";
  }
};

template <>
struct container_name<IstlRowMajorSparseMatrix<double>>
{
  static std::string value()
  {
    return "istl_sparse_matrix";
  }
};


} // namespace bindings
} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_CONTAINER_BINDINGS_HH
