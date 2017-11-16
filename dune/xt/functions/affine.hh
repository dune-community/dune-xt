// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#ifndef DUNE_XT_FUNCTIONS_AFFINE_HH
#define DUNE_XT_FUNCTIONS_AFFINE_HH

#include "config.h"

#if HAVE_DUNE_XT_LA

#include <memory>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/fmatrix.hh>

#include <dune/xt/functions/constant.hh>

#include <dune/xt/la/container/common.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {
namespace internal {


template <class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class AffineFunctionBase
{
  typedef AffineFunctionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  typedef FieldVector<DomainFieldImp, domainDim> DomainType;
  typedef typename XT::Functions::RangeTypeSelector<RangeFieldImp, rangeDim, rangeDimCols>::type RangeType;
  typedef typename XT::Functions::JacobianRangeTypeSelector<domainDim, RangeFieldImp, rangeDim, rangeDimCols>::type
      JacobianRangeType;
  typedef typename XT::Functions::RangeTypeSelector<RangeFieldImp, rangeDim, 1>::type ColRangeType;
  typedef typename XT::Functions::JacobianRangeTypeSelector<domainDim, RangeFieldImp, rangeDim, 1>::type
      ColJacobianRangeType;
  typedef typename LA::CommonSparseOrDenseMatrixCsr<RangeFieldImp> MatrixType;
  typedef FieldMatrix<RangeFieldImp, rangeDim, domainDim> FieldMatrixType;

  static std::string static_id()
  {
    return "affine";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["A"] = internal::EyeMatrix<RangeFieldImp, rangeDim, domainDim>::value_str();
    config["b"] = internal::EyeMatrix<RangeFieldImp, rangeDim, 1>::value_str();
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  // constructor
  explicit AffineFunctionBase(const std::vector<MatrixType>& A,
                              const RangeType& b = RangeType(0),
                              const std::string name_in = static_id())
    : A_(A)
    , b_(b)
    , name_(name_in)
    , b_zero_(Common::FloatCmp::eq(b_, RangeType(0)))
  {
    assert(A.size() >= rangeDimCols);
  }

  template <class MatrixImp>
  explicit AffineFunctionBase(const FieldVector<MatrixImp, rangeDimCols>& A,
                              const RangeType& b = RangeType(0),
                              const bool prune = true,
                              const std::string name_in = static_id())
    : A_(A.size())
    , b_(b)
    , name_(name_in)
    , b_zero_(Common::FloatCmp::eq(b_, RangeType(0)))
  {
    for (size_t cc = 0; cc < rangeDimCols; ++cc)
      A_[cc] = MatrixType(A[cc], prune);
  }

  explicit AffineFunctionBase(const std::vector<FieldMatrixType>& A,
                              const RangeType& b = RangeType(0),
                              const bool prune = true,
                              const std::string name_in = static_id())
    : A_(A.size())
    , b_(b)
    , name_(name_in)
    , b_zero_(Common::FloatCmp::eq(b_, RangeType(0)))
  {
    assert(A.size() >= rangeDimCols);
    for (size_t cc = 0; cc < rangeDimCols; ++cc)
      A_[cc] = MatrixType(A[cc], prune);
  }

  // constructor for dimRangeCols = 1.
  explicit AffineFunctionBase(const MatrixType& A,
                              const RangeType& b = RangeType(0),
                              const std::string name_in = static_id())
    : AffineFunctionBase(std::vector<MatrixType>(1, A), b, name_in)
  {
    static_assert(rangeDimCols == 1, "Use constructor above for dimRangeCols > 1");
  }

  // constructor for dimRangeCols = 1.
  explicit AffineFunctionBase(const DynamicMatrix<RangeFieldImp>& A,
                              const RangeType& b = RangeType(0),
                              const bool prune = true,
                              const std::string name_in = static_id())
    : AffineFunctionBase(std::vector<MatrixType>(1, MatrixType(A, prune)), b, name_in)
  {
    static_assert(rangeDimCols == 1, "Use constructor above for dimRangeCols > 1");
  }

  AffineFunctionBase(const ThisType& other) = default;

// if HAVE_DUNE_FEM is true, GlobalFunctionInterface is derived from Fem::Function which has a deleted copy assignment
// operator
#if HAVE_DUNE_FEM
  ThisType& operator=(const ThisType& other)
  {
    A_ = other.A_;
    b_ = other.b_;
    b_zero_ = other.b_zero_;
    name_ = other.name_;
    return *this;
  }
#endif

  const std::vector<MatrixType>& A() const
  {
    return A_;
  }

  const RangeType& b() const
  {
    return b_;
  }

protected:
  template <bool is_not_tensor = (rangeDimCols == 1), class anything = void>
  struct helper
  {
    static void evaluate(
        const std::vector<MatrixType>& A, const RangeType& b, const bool b_zero, const DomainType& x, RangeType& ret)
    {
      A[0].mv(x, ret);
      if (!b_zero)
        ret += b;
    }

    static void evaluate_col(const size_t col,
                             const std::vector<MatrixType>& A,
                             const RangeType& b,
                             const bool b_zero,
                             const DomainType& x,
                             RangeType& ret)
    {
      assert(col == 0);
      evaluate(A, b, b_zero, x, ret);
    }

    static void jacobian(const std::vector<MatrixType>& A, JacobianRangeType& ret)
    {
      ret = *(A[0].operator std::unique_ptr<FieldMatrixType>());
    }

    static void jacobian_col(const size_t col, const std::vector<MatrixType>& A, JacobianRangeType& ret)
    {
      assert(col == 0);
      jacobian(A, ret);
    }
  }; // struct helper<true, ...>

  template <class anything>
  struct helper<false, anything>
  {
    static void evaluate(
        const std::vector<MatrixType>& A, const RangeType& b, const bool b_zero, const DomainType& x, RangeType& ret)
    {
      for (size_t cc = 0; cc < rangeDimCols; ++cc) {
        Dune::FieldVector<RangeFieldImp, rangeDim> tmp_col;
        A[cc].mv(x, tmp_col);
        for (size_t rr = 0; rr < rangeDim; ++rr)
          ret[rr][cc] = tmp_col[rr];
      }
      if (!b_zero)
        ret += b;
    }

    static void evaluate_col(const size_t col,
                             const std::vector<MatrixType>& A,
                             const RangeType& b,
                             const bool b_zero,
                             const DomainType& x,
                             ColRangeType& ret)
    {
      A[col].mv(x, ret);
      if (!b_zero)
        for (size_t rr = 0; rr < rangeDim; ++rr)
          ret += b[rr][col];
    }

    static void jacobian(const std::vector<MatrixType>& A, JacobianRangeType& ret)
    {
      for (size_t cc = 0; cc < rangeDimCols; ++cc)
        ret[cc] = *(A[cc].operator std::unique_ptr<FieldMatrixType>());
    }

    static void jacobian_col(const size_t col, const std::vector<MatrixType>& A, ColJacobianRangeType& ret)
    {
      ret = *(A[col].operator std::unique_ptr<FieldMatrixType>());
    }
  }; // struct helper<false, ...>

protected:
  std::vector<MatrixType> A_;
  RangeType b_;
  bool b_zero_;
  std::string name_;
}; // class AffineFunctionBase<...>


} // namespace internal


/**
 * \brief Simple affine function of the form f(x) = A*x + b. For dimRangeCols > 1, there has to be a matrix A_i for
 * every column.
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class AffineFunction
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>,
      public internal::AffineFunctionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      InterfaceType;
  typedef internal::AffineFunctionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
  typedef AffineFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  using typename InterfaceType::DomainType;
  using typename InterfaceType::RangeFieldType;
  using typename InterfaceType::RangeType;
  using typename InterfaceType::JacobianRangeType;
  using InterfaceType::dimDomain;
  using InterfaceType::dimRange;
  using InterfaceType::dimRangeCols;
  using typename BaseType::MatrixType;
  using typename BaseType::FieldMatrixType;

  using BaseType::default_config;

  static std::string static_id()
  {
    return BaseType::static_id() + ".affine";
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
  {
    return 1;
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    std::vector<MatrixType> A_vector(dimRangeCols);
    for (size_t cc = 0; cc < dimRangeCols; ++cc) {
      if (cc == 0 && cfg.has_key("A")) {
        A_vector[0] = cfg.get<MatrixType>("A", dimRange, dimDomain);
      } else {
        A_vector[cc] = cfg.get<MatrixType>("A." + Common::to_string(cc), dimRange, dimDomain);
      }
    }
    return Common::make_unique<ThisType>(
        A_vector, cfg.get<RangeType>("b"), cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  // constructor
  explicit AffineFunction(const std::vector<MatrixType>& A,
                          const RangeType& b = RangeType(0),
                          const std::string name_in = static_id())
    : BaseType(A, b, name_in)
  {
  }

  explicit AffineFunction(const FieldVector<FieldMatrixType, dimRangeCols>& A,
                          const RangeType& b = RangeType(0),
                          const bool prune = true,
                          const std::string name_in = static_id())
    : BaseType(A, b, prune, name_in)
  {
  }

  explicit AffineFunction(const std::vector<FieldMatrixType>& A,
                          const RangeType& b = RangeType(0),
                          const bool prune = true,
                          const std::string name_in = static_id())
    : BaseType(A, b, prune, name_in)
  {
  }

  // constructor for dimRangeCols = 1.
  explicit AffineFunction(const MatrixType& A,
                          const RangeType& b = RangeType(0),
                          const std::string name_in = static_id())
    : BaseType(A, b, name_in)
  {
  }

  using InterfaceType::evaluate;

  virtual void evaluate(const DomainType& x, RangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    BaseType::template helper<>::evaluate(A_, b_, b_zero_, x, ret, mu);
  }

  using InterfaceType::jacobian;

  virtual void
  jacobian(const DomainType& /*x*/, JacobianRangeType& ret, const Common::Parameter& mu = {}) const override final
  {
    BaseType::template helper<>::jacobian(A_, ret, mu);
  }

private:
  using BaseType::A_;
  using BaseType::b_;
  using BaseType::b_zero_;
}; // class AffineFunction<...>


/**
 * \brief Simple affine function of the form f(x, u) = A*u + b. For dimRangeCols > 1, there has to be a matrix A_i for
 * every column.
 */
template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class U_,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class AffineFluxFunction
    : public GlobalFluxFunctionInterface<EntityImp,
                                         DomainFieldImp,
                                         domainDim,
                                         U_,
                                         0,
                                         RangeFieldImp,
                                         rangeDim,
                                         rangeDimCols>,
      public internal::
          AffineFunctionBase<typename U_::RangeFieldType, U_::dimRange, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFluxFunctionInterface<EntityImp,
                                      DomainFieldImp,
                                      domainDim,
                                      U_,
                                      0,
                                      RangeFieldImp,
                                      rangeDim,
                                      rangeDimCols>
      InterfaceType;
  typedef internal::AffineFunctionBase<typename U_::RangeFieldType, U_::dimRange, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef AffineFluxFunction ThisType;

public:
  typedef typename InterfaceType::DomainType DomainType;
  typedef typename InterfaceType::RangeFieldType RangeFieldType;
  typedef typename InterfaceType::RangeType RangeType;
  typedef typename InterfaceType::PartialURangeType PartialURangeType;
  typedef typename InterfaceType::ColRangeType ColRangeType;
  typedef typename InterfaceType::ColPartialURangeType ColPartialURangeType;
  using InterfaceType::dimDomain;
  using InterfaceType::dimRange;
  using InterfaceType::dimRangeCols;
  typedef typename InterfaceType::StateType StateType;
  typedef typename InterfaceType::StateRangeType StateRangeType;
  typedef typename InterfaceType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::MatrixType MatrixType;
  typedef typename BaseType::FieldMatrixType FieldMatrixType;

  using BaseType::default_config;

  static std::string static_id()
  {
    return BaseType::static_id() + ".affineflux";
  }

  virtual size_t order(const Common::Parameter& /*param*/ = {}) const override
  {
    return 1;
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    std::vector<MatrixType> A_vector(dimRangeCols);
    for (size_t cc = 0; cc < dimRangeCols; ++cc) {
      if (cc == 0 && cfg.has_key("A")) {
        A_vector[0] = cfg.get<MatrixType>("A", dimRange, StateType::dimRange);
      } else {
        A_vector[cc] = cfg.get<MatrixType>("A." + Common::to_string(cc), dimRange, StateType::dimRange);
      }
    }
    return Common::make_unique<ThisType>(
        A_vector, cfg.get<RangeType>("b"), cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  // constructor
  explicit AffineFluxFunction(const std::vector<MatrixType>& A,
                              const RangeType& b = RangeType(0),
                              const std::string name_in = static_id())
    : BaseType(A, b, name_in)
  {
  }

  template <class MatrixImp>
  explicit AffineFluxFunction(const Dune::FieldVector<MatrixImp, dimRangeCols>& A,
                              const RangeType& b = RangeType(0),
                              const bool prune = true,
                              const std::string name_in = static_id())
    : BaseType(A, b, prune, name_in)
  {
  }

  explicit AffineFluxFunction(const std::vector<FieldMatrixType>& A,
                              const RangeType& b = RangeType(0),
                              const bool prune = true,
                              const std::string name_in = static_id())
    : BaseType(A, b, prune, name_in)
  {
  }

  // constructor for dimRangeCols = 1.
  explicit AffineFluxFunction(const MatrixType& A,
                              const RangeType& b = RangeType(0),
                              const std::string name_in = static_id())
    : BaseType(A, b, name_in)
  {
  }

  explicit AffineFluxFunction(const DynamicMatrix<RangeFieldType>& A,
                              const RangeType& b = RangeType(0),
                              const bool prune = true,
                              const std::string name_in = static_id())
    : BaseType(A, b, prune, name_in)
  {
  }

  using InterfaceType::evaluate;

  virtual void evaluate(const DomainType& /*x*/,
                        const StateRangeType& u,
                        RangeType& ret,
                        const Common::Parameter& /*mu*/ = {}) const override final
  {
    BaseType::template helper<>::evaluate(A_, b_, b_zero_, u, ret);
  }

  virtual void evaluate_col(const size_t col,
                            const DomainType& /*x*/,
                            const StateRangeType& u,
                            ColRangeType& ret,
                            const Common::Parameter& /*mu*/ = {}) const override final
  {
    BaseType::template helper<>::evaluate_col(col, A_, b_, b_zero_, u, ret);
  }

  using InterfaceType::partial_u;

  virtual void partial_u(const DomainType& /*x*/,
                         const StateRangeType& /*u*/,
                         PartialURangeType& ret,
                         const Common::Parameter& /*mu*/ = {}) const override final
  {
    BaseType::template helper<>::jacobian(A_, ret);
  }

  virtual void partial_u_col(const size_t col,
                             const DomainType& /*x*/,
                             const StateRangeType& /*u*/,
                             ColPartialURangeType& ret,
                             const Common::Parameter& /*mu*/ = {}) const override final
  {
    BaseType::template helper<>::jacobian_col(col, A_, ret);
  }

private:
  using BaseType::A_;
  using BaseType::b_;
  using BaseType::b_zero_;
}; // class AffineFluxFunction<...>


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_XT_LA

#endif // DUNE_XT_FUNCTIONS_AFFINE_HH
