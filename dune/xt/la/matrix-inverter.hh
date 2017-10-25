// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_MATRIX_INVERSE_HH
#define DUNE_XT_LA_MATRIX_INVERSE_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/configuration.hh>


namespace Dune {
namespace XT {
namespace LA {


/**
 * \brief A means to obtain available options at compile time.
 * \note  This class needs to bespecialized for each MatrixType, the purpose of this class is merely to document the
 *        expected functionality.
 */
template <class MatrixType>
class InverseMatrixOptions
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

  static std::vector<std::string> types()
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  static Common::Configuration options(const std::string /*type*/ = "")
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }
}; // class InverseMatrixOptions


/**
 * \brief A means to obtain the inverse of a given matrix.
 * \note  This class needs to bespecialized for each MatrixType, the purpose of this class is merely to document the
 *        expected functionality.
 */
template <class MatrixType>
class InverseMatrix
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

  InverseMatrix(const MatrixType& /*matrix*/, const std::string& /*type*/ = "")
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  InverseMatrix(const MatrixType& /*matrix*/, const Common::Configuration /*opts*/)
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  const Common::Configuration& options() const
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  const MatrixType& matrix() const
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  const MatrixType& inverse() const
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  MatrixType inverse()
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }
}; // class InverseMatrix


template <class M>
InverseMatrix<M> make_inverse_matrix(const M& matrix, const std::string& type = "")
{
  return InverseMatrix<M>(matrix, type);
}


template <class M>
InverseMatrix<M> make_inverse_matrix(const M& matrix, const XT::Common::Configuration& options)
{
  return InverseMatrix<M>(matrix, options);
}


} // namespace Dune
} // namespace XT
} // namespace LA

#include "matrix-inverse/eigen.hh"


#endif // DUNE_XT_LA_MATRIX_INVERSE_HH
