// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2018)

#include "config.h"

#include <string>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>

#include <python/dune/xt/la/container/container-interface.hh>
#include <python/dune/xt/la/container/vector-interface.hh>
#include <python/dune/xt/la/container/pattern.hh>
#include <python/dune/xt/la/container/matrix-interface.hh>
#include <python/dune/xt/la/solver.hh>

#include <dune/xt/la/container.hh>


PYBIND11_MODULE(_la, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  namespace LA = Dune::XT::LA;

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");

  LA::bind_Backends(m);

  auto common_dense_vector_double = LA::bind_Vector<LA::CommonDenseVector<double>>(m);
#if HAVE_DUNE_ISTL
  auto istl_dense_vector_double = LA::bind_Vector<LA::IstlDenseVector<double>>(m);
#endif
#if HAVE_EIGEN
  auto eigen_dense_vector_double = LA::bind_Vector<LA::EigenDenseVector<double>>(m);
#endif

  LA::bind_SparsityPatternDefault(m);

#define BIND_MATRIX(C, s, c) auto c = LA::bind_Matrix<C, s>(m);

  BIND_MATRIX(LA::CommonDenseMatrix<double>, false, common_dense_matrix_double);
//  BIND_MATRIX(LA::CommonSparseMatrix<double>, true, common_sparse_matrix_double);
#if HAVE_DUNE_ISTL
  BIND_MATRIX(LA::IstlRowMajorSparseMatrix<double>, true, istl_row_major_sparse_matrix_double);
#endif
#if HAVE_EIGEN
  //  BIND_MATRIX(LA::EigenDenseMatrix<double>, false, eigen_dense_matrix_double);
  BIND_MATRIX(LA::EigenRowMajorSparseMatrix<double>, true, eigen_row_major_sparse_matrix_double);
#endif
#undef BIND_MATRIX
  LA::addbind_Matrix_Vector_interaction(common_dense_matrix_double, common_dense_vector_double);
//  LA::addbind_Matrix_Vector_interaction(common_sparse_matrix_double, common_dense_vector_double);
#if HAVE_DUNE_ISTL
  LA::addbind_Matrix_Vector_interaction(istl_row_major_sparse_matrix_double, istl_dense_vector_double);
#endif
#if HAVE_EIGEN
  //  LA::addbind_Matrix_Vector_interaction(eigen_dense_matrix_double, eigen_dense_vector_double);
  LA::addbind_Matrix_Vector_interaction(eigen_row_major_sparse_matrix_double, eigen_dense_vector_double);
#endif

  LA::bind_Solver<LA::CommonDenseMatrix<double>>(m);
//  LA::bind_Solver<LA::CommonSparseMatrix<double>>(m);
#if HAVE_DUNE_ISTL
  LA::bind_Solver<LA::IstlRowMajorSparseMatrix<double>>(m);
#endif
#if HAVE_EIGEN
  LA::bind_Solver<LA::EigenDenseMatrix<double>>(m);
  LA::bind_Solver<LA::EigenRowMajorSparseMatrix<double>>(m);
#endif

  m.def("_init_mpi",
        [](const std::vector<std::string>& args) {
          int argc = Dune::XT::Common::numeric_cast<int>(args.size());
          char** argv = Dune::XT::Common::vector_to_main_args(args);
          Dune::MPIHelper::instance(argc, argv);
        },
        "args"_a = std::vector<std::string>());

  m.def("_init_logger",
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
        "max_info_level"_a = std::numeric_limits<ssize_t>::max(),
        "max_debug_level"_a = std::numeric_limits<ssize_t>::max(),
        "enable_warnings"_a = true,
        "enable_colors"_a = true,
        "info_color"_a = "blue",
        "debug_color"_a = "darkgray",
        "warning_color"_a = "red");

  m.def("_test_logger",
        [](const bool info, const bool debug, const bool warning) {
          auto logger = Dune::XT::Common::TimedLogger().get("dune.xt.la");
          if (info)
            logger.info() << "info logging works!" << std::endl;
          if (debug)
            logger.debug() << "debug logging works!" << std::endl;
          if (warning)
            logger.warn() << "warning logging works!" << std::endl;
        },
        "info"_a = true,
        "debug"_a = true,
        "warning"_a = true);

} // PYBIND11_PLUGIN(la)
