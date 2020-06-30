// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_LA_CONTAINER_BINDINGS_HH
#define DUNE_XT_LA_CONTAINER_BINDINGS_HH

#include <dune/xt/la/container.hh>


// this is used in other headers
typedef Dune::XT::LA::CommonDenseVector<double> COMMON_DENSE_VECTOR;
typedef Dune::XT::LA::CommonDenseMatrix<double> COMMON_DENSE_MATRIX;
typedef Dune::XT::LA::CommonSparseMatrix<double> COMMON_SPARSE_MATRIX;
#if HAVE_EIGEN
typedef Dune::XT::LA::EigenDenseVector<double> EIGEN_DENSE_VECTOR;
typedef Dune::XT::LA::EigenDenseMatrix<double> EIGEN_DENSE_MATRIX;
typedef Dune::XT::LA::EigenRowMajorSparseMatrix<double> EIGEN_SPARSE_MATRIX;
#endif
typedef Dune::XT::LA::IstlDenseVector<double> ISTL_DENSE_VECTOR;
typedef Dune::XT::LA::IstlRowMajorSparseMatrix<double> ISTL_SPARSE_MATRIX;

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
