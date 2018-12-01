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

#ifndef DUNE_XT_LA_MATRIX_INVERTER_DEFAULT_HH
#define DUNE_XT_LA_MATRIX_INVERTER_DEFAULT_HH

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverter.hh>

#include "internal/base.hh"

namespace Dune {
namespace XT {
namespace LA {


template <class MatrixType>
class MatrixInverterOptions<MatrixType, true>
{
public:
  static std::vector<std::string> types()
  {
    return {"direct"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_matrix_inverter_type(actual_type, types());
    Common::Configuration opts = internal::default_matrix_inverter_options();
    opts["type"] = actual_type;
    return opts;
  }
}; // class MatrixInverterOptions<MatrixType, true>


template <class MatrixImp>
class MatrixInverter<MatrixImp, true> : public internal::MatrixInverterBase<MatrixImp>
{
  using BaseType = internal::MatrixInverterBase<MatrixImp>;

public:
  using MatrixType = typename BaseType::MatrixType;

  template <class... Args>
  explicit MatrixInverter(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
    const auto type = options_.template get<std::string>("type");
    const XT::Common::Configuration default_opts = MatrixInverterOptions<MatrixType>::options(type);
    if (!options_.get("delay_computation", default_opts.get<bool>("delay_computation")))
      compute();
  }

  void compute() override final
  {
    using M = Common::MatrixAbstraction<MatrixType>;
    const auto type = options_.template get<std::string>("type");
    if (type == "direct") {
      inverse_ = M::make_unique(M::rows(matrix_), M::cols(matrix_));
      auto tmp_matrix = M::make_unique(M::rows(matrix_), M::cols(matrix_));
      *tmp_matrix = matrix_;
      solve_by_qr_decomposition(*tmp_matrix, *inverse_, eye_matrix<MatrixType>(M::rows(matrix_)));
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type
                                << "' is none of MatrixInverterOptions<FieldMatrix<K, ROWS, COLS>>::types(), and "
                                   "internal::MatrixInverterBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << MatrixInverterOptions<MatrixType>::types());

    this->post_checks();
  } // ... compute(...)

protected:
  using BaseType::inverse_;
  using BaseType::matrix_;
  using BaseType::options_;
}; // class MatrixInverter<MatrixType, true>


} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_MATRIX_INVERTER_DEFAULT_HH
