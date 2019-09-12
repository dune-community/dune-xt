// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2017 - 2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_XT_LA_MATRIX_INVERTER_HH
#define DUNE_XT_LA_MATRIX_INVERTER_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/common/configuration.hh>

#include <dune/xt/la/type_traits.hh>

namespace Dune {
namespace XT {
namespace LA {


/**
 * \brief A means to obtain available options at compile time.
 * \note  This class needs to be specialized for each MatrixType, the purpose of this variant is merely to document the
 *        expected functionality.
 */
template <class MatrixType, bool is_matrix = Common::is_matrix<MatrixType>::value>
class MatrixInverterOptions
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

  static std::vector<std::string> types()
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
    return std::vector<std::string>();
  }

  static Common::Configuration options(const std::string /*type*/ = "")
  {
    static_assert(AlwaysFalse<MatrixType>::value,
                  "Please implement for given MatrixType and add the respective include below!");
    return Common::Configuration();
  }
}; // class MatrixInverterOptions


template <class MatrixType>
std::vector<std::string> matrix_inverter_types(const MatrixType& /*matrix*/)
{
  return MatrixInverterOptions<MatrixType>::types();
}


template <class MatrixType>
Common::Configuration matrix_inverter_options(const MatrixType& /*matrix*/, const std::string type = "")
{
  return MatrixInverterOptions<MatrixType>::options(type);
}


/**
 * \brief A means to obtain the inverse of a given matrix.
 * \note  This class needs to bespecialized for each MatrixType, the purpose of this class is merely to document the
 *        expected functionality.
 */
template <class MatrixType, bool is_matrix = Common::is_matrix<MatrixType>::value>
class MatrixInverter
{
  static_assert(AlwaysFalse<MatrixType>::value,
                "Please implement for given MatrixType and add the respective include below!");

public:
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

  const Common::Configuration& options() const;

  const MatrixType& matrix() const;

  const MatrixType& inverse() const;
}; // class MatrixInverter


template <class M>
MatrixInverter<M> make_matrix_inverter(const M& matrix, const std::string& type = "")
{
  return MatrixInverter<M>(matrix, type);
}


template <class M>
MatrixInverter<M> make_matrix_inverter(const M& matrix, const Common::Configuration& options)
{
  return MatrixInverter<M>(matrix, options);
}


template <class M>
typename std::enable_if<is_matrix<M>::value || Common::is_matrix<M>::value, M>::type
invert_matrix(const M& matrix, const std::string& inversion_type = "")
{
  return MatrixInverter<M>(matrix, inversion_type).inverse();
}


template <class M>
typename std::enable_if<is_matrix<M>::value || Common::is_matrix<M>::value, M>::type
invert_matrix(const M& matrix, const Common::Configuration& inversion_options)
{
  return MatrixInverter<M>(matrix, inversion_options).inverse();
}


} // namespace LA
} // namespace XT
} // namespace Dune

#include "matrix-inverter/default.hh"
#include "matrix-inverter/eigen.hh"
#include "matrix-inverter/fmatrix.hh"


#endif // DUNE_XT_LA_MATRIX_INVERTER_HH
