// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_STUFF_FUNCTIONS_EXPRESSION_HH
#define DUNE_STUFF_FUNCTIONS_EXPRESSION_HH

#include <vector>
#include <limits>

#include <dune/common/fvector.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>

#include "expression/base.hh"
#include "interfaces.hh"
#include "default.hh"
#include "constant.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


/**
 *  \attention  If you add the create() and default_config() method, do not forget to enable the matrix valued
 *              versions in test/function_expression.cc
 */
template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim,
          size_t rangeDimCols = 1>
class Expression
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
      BaseType;
  typedef Expression<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim * rangeDimCols>
      MathExpressionFunctionType;

  class Localfunction
      : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const size_t dimRange     = BaseType::dimRange;
    static const size_t dimRangeCols = BaseType::dimRangeCols;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const std::shared_ptr<const MathExpressionFunctionType>& function,
                  const size_t ord)
      : BaseType(ent)
      , function_(function)
      , order_(ord)
      , tmp_vector_(0)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const override
    {
      return order_;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      function_->evaluate(this->entity().geometry().global(xx), tmp_vector_);
      for (size_t ii = 0; ii < dimRange; ++ii) {
        auto& retRow = ret[ii];
        for (size_t jj = 0; jj < dimRangeCols; ++jj) {
          retRow[jj] = tmp_vector_[ii * dimRange + jj];
        }
      }
    } // ... evaluate(...)

    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const override
    {
      DUNE_THROW(NotImplemented,
                 "Once we decided on the JacobianRangeType of matrix valued functions we have to implement "
                     << "gradients for this function!");
    }

  private:
    const std::shared_ptr<const MathExpressionFunctionType> function_;
    const size_t order_;
    mutable FieldVector<RangeFieldType, dimRange * dimRangeCols> tmp_vector_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  Expression(const std::string variable, const std::string expression, const size_t ord = 0,
             const std::string nm = static_id())
    : function_(new MathExpressionFunctionType(variable, expression))
    , order_(ord)
    , name_(nm)
  {
  }

  Expression(const std::string variable, const std::vector<std::string> expressions, const size_t ord = 0,
             const std::string nm = static_id())
    : function_(new MathExpressionFunctionType(variable, expressions))
    , order_(ord)
    , name_(nm)
  {
  }

  Expression(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      function_ = other.function_;
      order_    = other.order_;
      name_     = other.name_;
    }
    return *this;
  }

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".expression";
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, function_, order_));
  }

private:
  std::shared_ptr<const MathExpressionFunctionType> function_;
  size_t order_;
  std::string name_;
}; // class Expression


template <class EntityImp, class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class Expression<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef Expression<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> MathExpressionFunctionType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, domainDim> MathExpressionGradientType;

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  static const size_t dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  static const size_t dimRange = BaseType::dimRange;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["variable"]   = "x";
    config["expression"] = "[x[0] sin(x[0]) exp(x[0])]";
    config["order"]      = "3";
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
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // get gradient
    std::vector<std::vector<std::string>> gradient_as_vectors;
    if (cfg.has_key("gradient")) {
      // get gradient as FieldMatrix
      typedef typename Dune::FieldMatrix<std::string, dimRange, dimDomain> JacobianMatrixType;
      const JacobianMatrixType gradient_as_matrix = cfg.get<JacobianMatrixType>("gradient");
      // convert FieldMatrix to std::vector< std::vector < std::string > >
      for (size_t rr = 0; rr < dimRange; ++rr) {
        std::vector<std::string> gradient_expression;
        for (size_t cc = 0; cc < dimDomain; ++cc)
          gradient_expression.emplace_back(gradient_as_matrix[rr][cc]);
        gradient_as_vectors.emplace_back(gradient_expression);
      }
    }
    // create
    return Common::make_unique<ThisType>(cfg.get("variable", default_cfg.get<std::string>("variable")),
                                         cfg.get("expression", default_cfg.get<std::vector<std::string>>("expression")),
                                         cfg.get("order", default_cfg.get<size_t>("order")),
                                         cfg.get("name", default_cfg.get<std::string>("name")),
                                         gradient_as_vectors);
  } // ... create(...)

  Expression(const std::string variable, const std::string expression,
             const size_t ord = default_config().get<size_t>("order"), const std::string nm = static_id(),
             const std::vector<std::vector<std::string>> gradient_expressions = std::vector<std::vector<std::string>>())
    : function_(new MathExpressionFunctionType(variable, expression))
    , order_(ord)
    , name_(nm)
  {
    build_gradients(variable, gradient_expressions);
  }

  Expression(const std::string variable, const std::vector<std::string> expressions,
             const size_t ord = default_config().get<size_t>("order"), const std::string nm = static_id(),
             const std::vector<std::vector<std::string>> gradient_expressions = std::vector<std::vector<std::string>>())
    : function_(new MathExpressionFunctionType(variable, expressions))
    , order_(ord)
    , name_(nm)
  {
    build_gradients(variable, gradient_expressions);
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual size_t order() const override
  {
    return order_;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const override
  {
    function_->evaluate(xx, ret);
#ifndef NDEBUG
#ifndef DUNE_STUFF_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
    bool failure = false;
    std::string type;
    if (std::isnan(ret[0])) {
      failure = true;
      type = "NaN";
    } else if (std::isinf(ret[0])) {
      failure = true;
      type = "inf";
    } else if (std::abs(ret[0]) > (0.9 * std::numeric_limits<double>::max())) {
      failure = true;
      type    = "an unlikely value";
    }
    if (failure)
      DUNE_THROW(Stuff::Exceptions::internal_error,
                 "evaluating this function yielded "
                     << type
                     << "!\n"
                     << "The variable of this function is:     "
                     << function_->variable()
                     << "\n"
                     << "The expression of this functional is: "
                     << function_->expression().at(0)
                     << "\n"
                     << "You tried to evaluate it with:   xx = "
                     << xx
                     << "\n"
                     << "The result was:                       "
                     << ret
                     << "\n\n"
                     << "You can disable this check by defining DUNE_STUFF_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
#endif // DUNE_STUFF_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
#endif // NDEBUG
  } // ... evaluate(...)

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override
  {
    if (gradients_.size() == 0)
      DUNE_THROW(NotImplemented, "This function does not provide any gradients!");
    assert(gradients_.size() == dimRange);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      gradients_[ii]->evaluate(xx, ret[ii]);
    }
  } // ... jacobian(...)

private:
  void build_gradients(const std::string variable, const std::vector<std::vector<std::string>>& gradient_expressions)
  {
    assert(gradient_expressions.size() == 0 || gradient_expressions.size() >= dimRange);
    if (gradient_expressions.size() > 0)
      for (size_t rr = 0; rr < dimRange; ++rr) {
        const auto& gradient_expression = gradient_expressions[rr];
        assert(gradient_expression.size() >= dimDomain);
        gradients_.emplace_back(new MathExpressionGradientType(variable, gradient_expression));
      }
  } // ... build_gradients(...)

  std::shared_ptr<const MathExpressionFunctionType> function_;
  size_t order_;
  std::string name_;
  std::vector<std::shared_ptr<const MathExpressionGradientType>> gradients_;
}; // class Expression< ..., 1 >


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_EXPRESSION_HH
