// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2016)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2015)
//   Tobias Leibner  (2014 - 2015)

#ifndef DUNE_XT_FUNCTIONS_AFFINE_HH
#define DUNE_XT_FUNCTIONS_AFFINE_HH

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
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
  typedef AffineFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  typedef typename LA::CommonSparseMatrix<RangeFieldType> MatrixType;
  typedef typename Dune::FieldMatrix<RangeFieldType, dimRange, dimDomain> FieldMatrixType;

  using typename BaseType::LocalfunctionType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".affine";
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
  explicit AffineFunction(const std::vector<MatrixType> A,
                          const RangeType b = RangeType(0),
                          const std::string name_in = static_id())
    : A_(A)
    , b_(b)
    , name_(name_in)
    , b_zero_(Common::FloatCmp::eq(b_, RangeType(0)))
  {
    assert(A.size() >= dimRangeCols);
  }

  explicit AffineFunction(const Dune::FieldVector<FieldMatrixType, dimRangeCols> A,
                          const RangeType b = RangeType(0),
                          const bool prune = true,
                          const std::string name_in = static_id())
    : A_(A.size())
    , b_(b)
    , name_(name_in)
    , b_zero_(Common::FloatCmp::eq(b_, RangeType(0)))
  {
    for (size_t cc = 0; cc < dimRangeCols; ++cc)
      A_[cc] = MatrixType(A[cc], prune);
  }

  // constructor for dimRangeCols = 1.
  explicit AffineFunction(const MatrixType A, const RangeType b = RangeType(0), const std::string name_in = static_id())
    : AffineFunction(std::vector<MatrixType>(1, A), b, name_in)
  {
    static_assert(dimRangeCols == 1, "Use constructor above for dimRangeCols > 1");
  }

  AffineFunction(const ThisType& other) = default;

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

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".affine";
  }

  virtual size_t order() const override final
  {
    return 1;
  }

  using BaseType::evaluate;

  virtual void evaluate(const DomainType& x, RangeType& ret) const override final
  {
    evaluate_helper(x, ret, internal::ChooseVariant<dimRangeCols>());
  }

  using BaseType::jacobian;

  virtual void jacobian(const DomainType& x, JacobianRangeType& ret) const override final
  {
    jacobian_helper(x, ret, internal::ChooseVariant<dimRangeCols>());
  }

  virtual std::string name() const override final
  {
    return name_;
  }

private:
  template <size_t rC>
  void evaluate_helper(const DomainType& x, RangeType& ret, const internal::ChooseVariant<rC>) const
  {
    for (size_t cc = 0; cc < rC; ++cc) {
      Dune::FieldVector<RangeFieldType, dimRange> tmp_col;
      A_[cc].mv(x, tmp_col);
      for (size_t rr = 0; rr < dimRange; ++rr)
        ret[rr][cc] = tmp_col[rr];
    }
    if (!b_zero_)
      ret += b_;
  }

  void evaluate_helper(const DomainType& x, RangeType& ret, const internal::ChooseVariant<1>) const
  {
    A_[0].mv(x, ret);
    if (!b_zero_)
      ret += b_;
  }

  template <size_t rC>
  void jacobian_helper(const DomainType& /*x*/, JacobianRangeType& ret, const internal::ChooseVariant<rC>) const
  {
    for (size_t cc = 0; cc < rC; ++cc)
      ret[cc] = A_[cc].operator FieldMatrixType();
  }

  void jacobian_helper(const DomainType& /*x*/, JacobianRangeType& ret, const internal::ChooseVariant<1>) const
  {
    ret = A_[0].operator FieldMatrixType();
  }

  std::vector<MatrixType> A_;
  RangeType b_;
  bool b_zero_;
  std::string name_;
};

} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // HAVE_DUNE_XT_LA

#endif // DUNE_XT_FUNCTIONS_AFFINE_HH
