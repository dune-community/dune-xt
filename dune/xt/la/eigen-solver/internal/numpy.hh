// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_NUMPY_HH
#define DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_NUMPY_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/pybindxi/interpreter.hh>

#include <dune/xt/common/vector.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/exceptions.hh>

#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/fmatrix.hh>

namespace Dune {
namespace XT {
namespace LA {
namespace internal {


bool numpy_eigensolver_available();


/**
 * \todo Extend this for all matrix types M where is_matrix<M>::value or Common::is_matrix<M>::value is true, which
 *       would require pybind11 type casting for these (should happen in <dune/xt/common/matrix/bindings/hh>).
 */
template <class K, int SIZE>
void compute_eigenvalues_and_right_eigenvectors_of_a_fieldmatrix_using_numpy(
    const FieldMatrix<K, SIZE, SIZE>& matrix,
    std::vector<Common::complex_t<K>>& eigenvalues,
    FieldMatrix<Common::complex_t<K>, SIZE, SIZE>& right_eigenvectors)
{
  if (!numpy_eigensolver_available())
    DUNE_THROW(Exceptions::eigen_solver_failed_bc_it_was_not_set_up_correctly,
               "Do not call me if numpy_eigensolver_available() is false!");
  PybindXI::GlobalInterpreter().import_module("warnings").attr("filterwarnings")("error");
  PybindXI::GlobalInterpreter().import_module("numpy").attr("seterr")("raise");
  auto eig = PybindXI::GlobalInterpreter().import_module("numpy.linalg").attr("eig");
  try {
    auto result =
        eig(matrix)
            .template cast<
                std::tuple<FieldVector<Common::complex_t<K>, SIZE>, FieldMatrix<Common::complex_t<K>, SIZE, SIZE>>>();
    eigenvalues = Common::convert_to<std::vector<Common::complex_t<K>>>(std::get<0>(result));
    right_eigenvectors = std::get<1>(result);
  } catch (const pybind11::cast_error&) {
    auto result =
        eig(matrix)
            .template cast<
                std::tuple<FieldVector<Common::real_t<K>, SIZE>, FieldMatrix<Common::real_t<K>, SIZE, SIZE>>>();
    eigenvalues = Common::convert_to<std::vector<Common::complex_t<K>>>(std::get<0>(result));
    right_eigenvectors = Common::convert_to<FieldMatrix<Common::complex_t<K>, SIZE, SIZE>>(std::get<1>(result));
  } catch (const std::runtime_error& ee) {
    DUNE_THROW(Exceptions::eigen_solver_failed,
               "Could not convert result!\n\nThis was the original error: " << ee.what());
  }
} // ... compute_eigenvalues_and_right_eigenvectors_of_a_fieldmatrix_using_numpy(...)


} // namespace internal
} // namespace LA
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_LA_EIGEN_SOLVER_INTERNAL_NUMPY_HH
