// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <string>
#include <vector>

#include <dune/pybindxi/pybind11.h>
//#include <dune/pybindxi/stl_bind.h> // <- see dune/xt/common/bindings.cc

#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include "container/container-interface.pbh"
#include "container/vector-interface.pbh"
#include "container/pattern.pbh"
#include "container/matrix-interface.pbh"

#include "container.hh"


namespace py = pybind11;


PYBIND11_PLUGIN(la)
{
  py::module m("la", "dune-xt-la");

  py::module::import("common");

  Dune::XT::LA::bind_Backends(m);

#define BIND_VECTOR(C, c, id)                                                                                          \
  auto c = Dune::XT::LA::bind_Vector<C>(m, id);                                                                        \
  Dune::XT::LA::addbind_ProvidesBackend(c);                                                                            \
  Dune::XT::LA::addbind_ProvidesDataAccess(c)

  BIND_VECTOR(Dune::XT::LA::CommonDenseVector<double>, common_dense_vector_double, "CommonDenseVector_double");
#if HAVE_DUNE_ISTL
  BIND_VECTOR(Dune::XT::LA::IstlDenseVector<double>, istl_dense_vector_double, "IstlDenseVector_double");
#endif
#if HAVE_EIGEN
  BIND_VECTOR(Dune::XT::LA::EigenDenseVector<double>, eigen_dense_vector_double, "EigenDenseVector_double");
#endif
#undef BIND_VECTOR

  Dune::XT::LA::bind_SparsityPatternDefault(m);

#define BIND_MATRIX(C, s, c, id)                                                                                       \
  auto c = Dune::XT::LA::bind_Matrix<C, s>(m, id);                                                                     \
  Dune::XT::LA::addbind_ProvidesBackend(c);

  BIND_MATRIX(Dune::XT::LA::CommonDenseMatrix<double>, false, common_dense_matrix_double, "CommonDenseMatrix_double");
  BIND_MATRIX(Dune::XT::LA::CommonSparseMatrix<double>, true, common_sparse_matrix_double, "CommonSparseMatrix_double");
#if HAVE_DUNE_ISTL
  BIND_MATRIX(Dune::XT::LA::IstlRowMajorSparseMatrix<double>,
              true,
              istl_row_major_sparse_matrix_double,
              "IstlRowMajorSparseMatrix_double");
#endif
#if HAVE_EIGEN
  BIND_MATRIX(Dune::XT::LA::EigenDenseMatrix<double>, false, eigen_dense_matrix_double, "EigenDenseMatrix_double");
  BIND_MATRIX(Dune::XT::LA::EigenRowMajorSparseMatrix<double>,
              true,
              eigen_row_major_sparse_matrix_double,
              "EigenRowMajorSparseMatrix_double");
#endif
#undef BIND_MATRIX
  Dune::XT::LA::addbind_Matrix_Vector_interaction<Dune::XT::LA::CommonDenseVector<double>>(common_dense_matrix_double);
  Dune::XT::LA::addbind_Matrix_Vector_interaction<Dune::XT::LA::CommonDenseVector<double>>(common_sparse_matrix_double);
#if HAVE_DUNE_ISTL
  Dune::XT::LA::addbind_Matrix_Vector_interaction<Dune::XT::LA::IstlDenseVector<double>>(
      istl_row_major_sparse_matrix_double);
#endif
#if HAVE_EIGEN
  Dune::XT::LA::addbind_Matrix_Vector_interaction<Dune::XT::LA::EigenDenseVector<double>>(eigen_dense_matrix_double);
  Dune::XT::LA::addbind_Matrix_Vector_interaction<Dune::XT::LA::EigenDenseVector<double>>(
      eigen_row_major_sparse_matrix_double);
#endif

  return m.ptr();
} // PYBIND11_PLUGIN(la)

#endif // HAVE_DUNE_PYBINDXI
