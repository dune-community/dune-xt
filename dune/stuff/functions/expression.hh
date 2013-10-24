#ifndef DUNE_STUFF_FUNCTIONS_EXPRESSION_HH
#define DUNE_STUFF_FUNCTIONS_EXPRESSION_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/tree.hh>

#include "expression/base.hh"
#include "interfaces.hh"
#include "default.hh"

namespace Dune {
namespace Stuff {
namespace Function {


/**
 *  \attention  If you add the create() and defaultSettings() method, do not forget to enable the matrix valued
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

    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& ret) const DS_OVERRIDE
    {
      DUNE_THROW(NotImplemented,
                 "If we decided on the JacobianRangeType of matrix valued functions we havo to implement gradients for "
                     << "this function!");
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
}; // class FunctionExpression


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Expression<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;
  typedef Expression<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> MathExpressionFunctionType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, domainDim> MathExpressionGradientType;

  class Localfunction : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

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

    Localfunction(const EntityType& entity, const std::shared_ptr<const MathExpressionFunctionType>& function,
                  std::vector<std::shared_ptr<const MathExpressionGradientType>> gradients, const size_t ord)
      : BaseType(entity)
      , function_(function)
      , gradients_(gradients)
      , order_(ord)
      , global_point_(DomainFieldType(0))
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
      assert(this->is_a_valid_point(xx));
      global_point_ = this->entity().geometry().global(xx);
      function_->evaluate(global_point_, ret);
      assert(this_value_is_sane(ret));
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
    {
      if (gradients_.size() == 0)
        DUNE_THROW(NotImplemented, "This function does not provide any gradients!");
      else {
        assert(gradients_.size() == dimRange);
        global_point_ = this->entity().geometry().global(xx);
        for (size_t ii = 0; ii < dimRange; ++ii) {
          gradients_[ii]->evaluate(global_point_, ret[ii]);
          assert(this_value_is_sane(ret[ii]));
        }
      }
    } // ... jacobian(...)

  private:
    template <class VectorType>
    bool this_value_is_sane(const VectorType& value) const
    {
      for (size_t rr = 0; rr < value.size(); ++rr)
        if (std::abs(value[rr]) > (0.9 * std::numeric_limits<RangeFieldType>::max())) {
          return false;
          //          std::stringstream ss;
          //          ss << "evaluating this function yielded an unlikely value!\n"
          //             << "The variable of this function is: " << function_->variable() << ",\n";
          //          Stuff::Common::print(function_->expression(),"the expression() of this function is", ss);
          //          Stuff::Common::print(xx, "you tried to evaluate it with xx", ss);
          //          Stuff::Common::print(value, "and the result was", ss);
          //          DUNE_THROW(InvalidStateException, ss.str());
        }
      return true;
    }

    const std::shared_ptr<const MathExpressionFunctionType> function_;
    std::vector<std::shared_ptr<const MathExpressionGradientType>> gradients_;
    const size_t order_;
    mutable DomainType global_point_;
  }; // class Localfunction

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

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["variable"]   = "x";
    description["expression"] = "[x[0]; sin(x[0]); exp(x[0])]";
    description["order"]      = "1";
    description["name"] = static_id();
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const DSC::ExtendedParameterTree settings = defaultSettings())
  {
    // get necessary
    const std::string _variable = settings.get<std::string>("variable", "x");
    std::vector<std::string> _expressions;
    // lets see, if there is a key or vector
    if (settings.hasVector("expression")) {
      const std::vector<std::string> expr = settings.getVector<std::string>("expression", 1);
      for (size_t ii = 0; ii < expr.size(); ++ii)
        _expressions.push_back(expr[ii]);
    } else if (settings.hasKey("expression")) {
      const std::string expr = settings.get<std::string>("expression");
      _expressions.push_back(expr);
    } else
      DUNE_THROW(Dune::IOError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                      << " neither key nor vector 'expression' found in the following settings:\n"
                      << settings.reportString("  "));
    // get optional
    const int order        = settings.get<int>("order", -1);
    const std::string name = settings.get<std::string>("name", "function.expression");
    // create and return
    return new ThisType(_variable, _expressions, order, name);
  } // ... create(...)

  Expression(const std::string variable, const std::string expression, const size_t ord = 0,
             const std::string nm                                             = static_id(),
             const std::vector<std::vector<std::string>> gradient_expressions = std::vector<std::vector<std::string>>())
    : function_(new MathExpressionFunctionType(variable, expression))
    , order_(ord)
    , name_(nm)
  {
    build_gradients(variable, gradient_expressions);
  }

  Expression(const std::string variable, const std::vector<std::string> expressions, const size_t ord = 0,
             const std::string nm                                             = static_id(),
             const std::vector<std::vector<std::string>> gradient_expressions = std::vector<std::vector<std::string>>())
    : function_(new MathExpressionFunctionType(variable, expressions))
    , order_(ord)
    , name_(nm)
  {
    build_gradients(variable, gradient_expressions);
  }

  Expression(const ThisType& other)
    : function_(other.function_)
    , order_(other.order_)
    , name_(other.name_)
    , gradients_(other.gradients_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      function_  = other.function_;
      order_     = other.order_;
      name_      = other.name_;
      gradients_ = other.gradients_;
    }
    return *this;
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, function_, gradients_, order_));
  }

private:
  void build_gradients(const std::string variable, const std::vector<std::vector<std::string>>& gradient_expressions)
  {
    assert(gradient_expressions.size() == 0 || gradient_expressions.size() == dimRange);
    for (const auto& gradient_expression : gradient_expressions) {
      assert(gradient_expression.size() == dimDomain);
      gradients_.emplace_back(new MathExpressionGradientType(variable, gradient_expression));
    }
  } // ... build_gradients(...)

  std::shared_ptr<const MathExpressionFunctionType> function_;
  size_t order_;
  std::string name_;
  std::vector<std::shared_ptr<const MathExpressionGradientType>> gradients_;
}; // class Expression< ..., 1 >


} // namespace Function
} // namespace Stuff
} // namespace Dune

#ifdef DUNE_STUFF_FUNCTIONS_TO_LIB
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
  extern template class Dune::Stuff::Function::Expression<etype, dftype, ddim, rftype, rdim, rcdim>;

#ifdef HAVE_DUNE_GRID

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
