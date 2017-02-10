// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_CONTAINER_BINDINGS_HH
#define DUNE_XT_LA_CONTAINER_BINDINGS_HH

#include "container.hh"

namespace Dune {


// this is used in other headers
typedef XT::LA::CommonDenseVector<double> COMMON_DENSE_VECTOR;
typedef XT::LA::CommonDenseMatrix<double> COMMON_DENSE_MATRIX;
typedef XT::LA::CommonSparseMatrix<double> COMMON_SPARSE_MATRIX;
#if HAVE_EIGEN
typedef XT::LA::EigenDenseVector<double> EIGEN_DENSE_VECTOR;
typedef XT::LA::EigenDenseMatrix<double> EIGEN_DENSE_MATRIX;
typedef XT::LA::EigenRowMajorSparseMatrix<double> EIGEN_SPARSE_MATRIX;
#endif
#if HAVE_DUNE_ISTL
typedef XT::LA::IstlDenseVector<double> ISTL_DENSE_VECTOR;
typedef XT::LA::IstlRowMajorSparseMatrix<double> ISTL_SPARSE_MATRIX;
#endif


} // namespace Dune


#endif // DUNE_XT_LA_CONTAINER_BINDINGS_HH
