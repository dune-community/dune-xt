// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_MATRIX_INVERTER_HH
#define DUNE_XT_LA_MATRIX_INVERTER_HH

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
class MatrixInverterOptions
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
}; // class MatrixInverterOptions


/**
 * \brief A means to obtain the inverse of a given matrix.
 * \note  This class needs to bespecialized for each MatrixType, the purpose of this class is merely to document the
 *        expected functionality.
 */
template <class MatrixType>
class MatrixInverter
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

  MatrixInverter(const MatrixType& /*matrix*/, const std::string& /*type*/ = "")
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
  }

  MatrixInverter(const MatrixType& /*matrix*/, const Common::Configuration /*opts*/)
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
}; // class MatrixInverter


template <class M>
MatrixInverter<M> make_matrix_inverter(const M& matrix, const std::string& type = "")
{
  return MatrixInverter<M>(matrix, type);
}


template <class M>
MatrixInverter<M> make_matrix_inverter(const M& matrix, const XT::Common::Configuration& options)
{
  return MatrixInverter<M>(matrix, options);
}


template <class M>
typename std::enable_if<is_matrix<M>::value, M>::type invert_matrix(const M& matrix,
                                                                    const std::string& inversion_type = "")
{
  return MatrixInverter<M>(matrix, inversion_type).inverse();
}


template <class M>
typename std::enable_if<is_matrix<M>::value, M>::type invert_matrix(const M& matrix,
                                                                    const XT::Common::Configuration& inversion_options)
{
  return MatrixInverter<M>(matrix, inversion_options).inverse();
}


} // namespace Dune
} // namespace XT
} // namespace LA

#include "matrix-inverter/eigen.hh"


#endif // DUNE_XT_LA_MATRIX_INVERTER_HH
