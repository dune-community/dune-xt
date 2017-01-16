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
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/configuration.pbh>
#include <dune/xt/common/fvector.pbh>

#include "container/container-interface.pbh"
#include "container/vector-interface.pbh"
#include "container/pattern.pbh"
#include "container/matrix-interface.pbh"
#include "solver.pbh"

#include "container.hh"


namespace py = pybind11;
namespace LA = Dune::XT::LA;
using namespace pybind11::literals;


PYBIND11_PLUGIN(la)
{
  py::module m("la", "dune-xt-la");

  py::module::import("common");

  LA::bind_Backends(m);

  LA::bind_Vector<LA::CommonDenseVector<double>>(m, "CommonDenseVector_double");
#if HAVE_DUNE_ISTL
  LA::bind_Vector<LA::IstlDenseVector<double>>(m, "IstlDenseVector_double");
#endif
#if HAVE_EIGEN
  LA::bind_Vector<LA::EigenDenseVector<double>>(m, "EigenDenseVector_double");
#endif

  LA::bind_SparsityPatternDefault(m);

#define BIND_MATRIX(C, s, c, id) auto c = LA::bind_Matrix<C, s>(m, id);

  BIND_MATRIX(LA::CommonDenseMatrix<double>, false, common_dense_matrix_double, "CommonDenseMatrix_double");
  BIND_MATRIX(LA::CommonSparseMatrix<double>, true, common_sparse_matrix_double, "CommonSparseMatrix_double");
#if HAVE_DUNE_ISTL
  BIND_MATRIX(LA::IstlRowMajorSparseMatrix<double>,
              true,
              istl_row_major_sparse_matrix_double,
              "IstlRowMajorSparseMatrix_double");
#endif
#if HAVE_EIGEN
  BIND_MATRIX(LA::EigenDenseMatrix<double>, false, eigen_dense_matrix_double, "EigenDenseMatrix_double");
  BIND_MATRIX(LA::EigenRowMajorSparseMatrix<double>,
              true,
              eigen_row_major_sparse_matrix_double,
              "EigenRowMajorSparseMatrix_double");
#endif
#undef BIND_MATRIX
  LA::addbind_Matrix_Vector_interaction<LA::CommonDenseVector<double>>(common_dense_matrix_double);
  LA::addbind_Matrix_Vector_interaction<LA::CommonDenseVector<double>>(common_sparse_matrix_double);
#if HAVE_DUNE_ISTL
  LA::addbind_Matrix_Vector_interaction<LA::IstlDenseVector<double>>(istl_row_major_sparse_matrix_double);
#endif
#if HAVE_EIGEN
  LA::addbind_Matrix_Vector_interaction<LA::EigenDenseVector<double>>(eigen_dense_matrix_double);
  LA::addbind_Matrix_Vector_interaction<LA::EigenDenseVector<double>>(eigen_row_major_sparse_matrix_double);
#endif

  LA::bind_Solver<LA::CommonDenseMatrix<double>>(m, "CommonDenseMatrix_double");
//  LA::bind_Solver<LA::CommonSparseMatrix<double>>(m, "CommonSparseMatrix_double");
#if HAVE_DUNE_ISTL
  LA::bind_Solver<LA::IstlRowMajorSparseMatrix<double>>(m, "IstlRowMajorSparseMatrix_double");
#endif
#if HAVE_EIGEN
  LA::bind_Solver<LA::EigenDenseMatrix<double>>(m, "EigenDenseMatrix");
  LA::bind_Solver<LA::EigenRowMajorSparseMatrix<double>>(m, "EigenRowMajorSparseMatrix");
#endif

  m.def("init_logger",
        [](const ssize_t max_info_level,
           const ssize_t max_debug_level,
           const bool enable_warnings,
           const bool enable_colors,
           const std::string& info_color,
           const std::string& debug_color,
           const std::string& warning_color) {
          Dune::XT::Common::TimedLogger().create(
              max_info_level, max_debug_level, enable_warnings, enable_colors, info_color, debug_color, warning_color);
        },
        "max_info_level"_a = -1,
        "max_debug_level"_a = -1,
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  return m.ptr();
} // PYBIND11_PLUGIN(la)

#endif // HAVE_DUNE_PYBINDXI
