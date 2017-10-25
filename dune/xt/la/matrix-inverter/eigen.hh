// This file is part of the dune-xt-la project:
//   https://github.com/dune-community/dune-xt-la
// Copyright 2009-2017 dune-xt-la developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_XT_LA_MATRIX_INVERSE_EIGEN_HH
#define DUNE_XT_LA_MATRIX_INVERSE_EIGEN_HH

#include <dune/common/typetraits.hh>

#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverse.hh>

#include "internal/base.hh"
#include "internal/eigen.hh"

namespace Dune {
namespace XT {
namespace LA {

#if HAVE_EIGEN


template <class S>
class InverseMatrixOptions<EigenDenseMatrix<S>>
{
public:
  static std::vector<std::string> types()
  {
    return {"direct", "moore_penrose"};
  }

  static Common::Configuration options(const std::string type = "")
  {
    const std::string actual_type = type.empty() ? types()[0] : type;
    internal::ensure_inverse_matrix_type(actual_type, types());
    Common::Configuration opts = internal::default_inverse_matrix_options();
    opts["type"] = actual_type;
    if (actual_type == "moore_penrose")
      opts["eigenvalues_tolerance_factor"] = "1e-15";
    return opts;
  }
}; // class InverseMatrixOptions<EigenDenseMatrix<S>>


template <class S>
class InverseMatrix<EigenDenseMatrix<S>> : public internal::InverseMatrixBase<EigenDenseMatrix<S>>
{
  using BaseType = internal::InverseMatrixBase<EigenDenseMatrix<S>>;

public:
  using MatrixType = typename BaseType::MatrixType;
  using typename BaseType::RealType;

  template <class... Args>
  explicit InverseMatrix(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
    const auto type = options_.template get<std::string>("type");
    const XT::Common::Configuration default_opts = InverseMatrixOptions<MatrixType>::options(type);
    if (!options_.get("delay_computation", default_opts.get<bool>("delay_computation")))
      compute();
  }

  void compute() override final
  {
    const auto type = options_.template get<std::string>("type");
    const XT::Common::Configuration default_opts = InverseMatrixOptions<MatrixType>::options(type);
    if (type == "direct")
      inverse_ = std::make_unique<MatrixType>(matrix_.backend().inverse());
    else if (type == "moore_penrose") {
      const auto eigenvalues_tolerance_factor =
          options_.get("eigenvalues_tolerance_factor", default_opts.get<RealType>("eigenvalues_tolerance_factor"));
      inverse_ = std::make_unique<MatrixType>(
          internal::compute_moore_penrose_inverse_using_eigen(matrix_.backend(), eigenvalues_tolerance_factor));
    } else
      DUNE_THROW(Common::Exceptions::internal_error,
                 "Given type '" << type << "' is none of InverseMatrixOptions<EigenDenseMatrix<S>>::types(), and "
                                           "internal::InverseMatrixBase promised to check this!"
                                << "\n\nThese are the available types:\n\n"
                                << InverseMatrixOptions<EigenDenseMatrix<S>>::types());

    this->post_checks();
  } // ... compute(...)

protected:
  using BaseType::matrix_;
  using BaseType::options_;
  using BaseType::inverse_;
}; // class InverseMatrix<EigenDenseMatrix<...>>


#else // HAVE_EIGEN


template <class S>
class InverseMatrixOptions<EigenDenseMatrix<S>>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


template <class S>
class InverseMatrix<EigenDenseMatrix<S>>
{
  static_assert(AlwaysFalse<S>::value, "You are missing eigen!");
};


#endif // HAVE_EIGEN

} // namespace Dune
} // namespace XT
} // namespace LA


#endif // DUNE_XT_LA_MATRIX_INVERSE_EIGEN_HH
