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

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/parameter/configcontainer.hh>
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
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
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

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const unsigned int dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const unsigned int dimRange     = BaseType::dimRange;
    static const unsigned int dimRangeCols = BaseType::dimRangeCols;
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

    virtual size_t order() const DS_OVERRIDE
    {
      return order_;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
    {
      function_->evaluate(this->entity().geometry().global(xx), tmp_vector_);
      for (size_t ii = 0; ii < dimRange; ++ii) {
        auto& retRow = ret[ii];
        for (size_t jj = 0; jj < dimRangeCols; ++jj) {
          retRow[jj] = tmp_vector_[ii * dimRange + jj];
        }
      }
    } // ... evaluate(...)

    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& /*ret*/) const DS_OVERRIDE
    {
      DUNE_THROW_COLORFULLY(NotImplemented,
                            "If we decided on the JacobianRangeType of matrix valued functions we have to implement "
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

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename LocalfunctionType::RangeType RangeType;

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

  Expression(const ThisType& other)
    : function_(other.function_)
    , order_(other.order_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      function_ = other.function_;
      order_    = other.order_;
      name_     = other.name_;
    }
    return *this;
  }

  virtual std::string type() const DS_OVERRIDE DS_FINAL
  {
    return BaseType::static_id() + ".expression";
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, function_, order_));
  }

private:
  std::shared_ptr<const MathExpressionFunctionType> function_;
  size_t order_;
  std::string name_;
}; // class Expression


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
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

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename LocalfunctionType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename LocalfunctionType::RangeType RangeType;

  typedef typename LocalfunctionType::JacobianRangeType JacobianRangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  static Common::ConfigContainer default_config(const std::string sub_name = "")
  {
    Common::ConfigContainer config;
    config["variable"]   = "x";
    config["expression"] = "[x[0] sin(x[0]) exp(x[0])]";
    config["order"]      = "3";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::ConfigContainer tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::ConfigContainer config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::ConfigContainer cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::ConfigContainer default_cfg = default_config();
    // create
    return Common::make_unique<ThisType>(cfg.get("variable", default_cfg.get<std::string>("variable")),
                                         cfg.get("expression", default_cfg.get<std::vector<std::string>>("expression")),
                                         cfg.get("order", default_cfg.get<size_t>("order")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
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

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  virtual size_t order() const DS_OVERRIDE
  {
    return order_;
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
  {
    function_->evaluate(xx, ret);
  }

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
  {
    if (gradients_.size() == 0)
      DUNE_THROW_COLORFULLY(NotImplemented, "This function does not provide any gradients!");
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

#if 0 // def DUNE_STUFF_FUNCTIONS_TO_LIB
the whole  DUNE_STUFF_FUNCTIONS_TO_LIB stuff threw loads of errors, so I disabled it everywhere (Rene)

#define DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(etype, ddim)                                                     \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGECOLS(etype, ddim, 1)                                                    \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGECOLS(etype, ddim, 2)                                                    \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGECOLS(etype, ddim, 3)

#define DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGECOLS(etype, ddim, rdim)                                           \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 2)                                          \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, 3)

#define DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DOMAINFIELDTYPES(etype, ddim, rdim, rcdim)                                \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_RANGEFIELDTYPES(etype, double, ddim, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_RANGEFIELDTYPES(etype, dftype, ddim, rdim, rcdim)                         \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LAST_EXPANSION(etype, dftype, ddim, double, rdim, rcdim)                             \
  DUNE_STUFF_FUNCTIONS_EXPRESSION_LAST_EXPANSION(etype, dftype, ddim, long double, rdim, rcdim)

#define DUNE_STUFF_FUNCTIONS_EXPRESSION_LAST_EXPANSION(etype, dftype, ddim, rftype, rdim, rcdim)                       \
  extern template class Dune::Stuff::Functions::Expression<etype, dftype, ddim, rftype, rdim, rcdim>;

#if HAVE_DUNE_GRID

DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesSGrid3dEntityType, 3)

DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid1dEntityType, 1)
DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesYaspGrid3dEntityType, 3)

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid2dEntityType, 2)
DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluSimplexGrid3dEntityType, 3)
DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE(DuneStuffFunctionsInterfacesAluCubeGrid3dEntityType, 3)

#endif // HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#endif // HAVE_DUNE_GRID

#undef DUNE_STUFF_FUNCTIONS_EXPRESSION_LAST_EXPANSION
#undef DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_RANGEFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DOMAINFIELDTYPES
#undef DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGECOLS
#undef DUNE_STUFF_FUNCTIONS_EXPRESSION_LIST_DIMRANGE
#endif // DUNE_STUFF_FUNCTIONS_TO_LIB

#endif // DUNE_STUFF_FUNCTIONS_EXPRESSION_HH
