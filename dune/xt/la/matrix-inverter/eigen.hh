// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2018 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_LA_MATRIX_INVERTER_EIGEN_HH
#define DUNE_XT_LA_MATRIX_INVERTER_EIGEN_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverter.hh>

#include "internal/base.hh"
#include "internal/eigen.hh"

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_EIGEN


template <class S>
class MatrixInverterOptions<EigenDenseMatrix<S>, true>
{
public:
  static std::vector<std::string> types()
  {
    return {"direct", "moore_penrose"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_matrix_inverter_type(actual_type, types());
    Common::Configuration opts = internal::default_matrix_inverter_options();
    opts["type"] = actual_type;
    if (actual_type == "moore_penrose")
      opts["eigenvalues_tolerance_factor"] = "1e-15";
    return opts;
  }
}; // class MatrixInverterOptions<EigenDenseMatrix<S>>


template <class S>
class MatrixInverter<EigenDenseMatrix<S>, true> : public internal::MatrixInverterBase<EigenDenseMatrix<S>>
{
  using BaseType = internal::MatrixInverterBase<EigenDenseMatrix<S>>;

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
    const auto type = options_.template get<std::string>("type");
    const XT::Common::Configuration default_opts = MatrixInverterOptions<MatrixType>::options(type);
    if (type == "direct")
      inverse_ = std::make_unique<MatrixType>(matrix_.backend().inverse());
    else if (type == "moore_penrose") {
      const auto eigenvalues_tolerance_factor =
          options_.get("eigenvalues_tolerance_factor", default_opts.get<double>("eigenvalues_tolerance_factor"));
      inverse_ = std::make_unique<MatrixType>(
          internal::compute_moore_penrose_inverse_using_eigen(matrix_.backend(), eigenvalues_tolerance_factor));
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type
                                << "' is none of MatrixInverterOptions<EigenDenseMatrix<S>>::types(), and "
                                   "internal::MatrixInverterBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << MatrixInverterOptions<EigenDenseMatrix<S>>::types());

    this->post_checks();
  } // ... compute(...)

protected:
  using BaseType::inverse_;
  using BaseType::matrix_;
  using BaseType::options_;
}; // class MatrixInverter<EigenDenseMatrix<...>>


#else // HAVE_EIGEN


template <class S>
class MatrixInverterOptions<EigenDenseMatrix<S>, true>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


template <class S>
class MatrixInverter<EigenDenseMatrix<S>, true>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


#endif // HAVE_EIGEN

} // namespace LA
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_LA_MATRIX_INVERTER_EIGEN_HH
