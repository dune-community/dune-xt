#ifndef DUNE_STUFF_FUNCTIONS_EXPRESSION_HH
#define DUNE_STUFF_FUNCTIONS_EXPRESSION_HH

#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/print.hh>

#include "expression/base.hh"
#include "interfaces.hh"
#include "default.hh"

namespace Dune {
namespace Stuff {
namespace Function {


template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows,
          int rangeDimCols = 1>
class Expression : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp,
                                                       rangeDimRows, rangeDimCols>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
      BaseType;
  typedef Expression<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> ThisType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows * rangeDimCols>
      MathExpressionFunctionType;

  class Localfunction
      : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const unsigned int dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const unsigned int dimRangeRows = BaseType::dimRangeRows;
    static const unsigned int dimRangeCols = BaseType::dimRangeCols;
    typedef typename BaseType::RangeType RangeType;

    //    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& ent, const std::shared_ptr<const MathExpressionFunctionType>& function,
                  const size_t ord)
      : entity_(ent)
      , function_(function)
      , order_(ord)
      , tmp_vector_(0)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual const EntityType& entity() const override
    {
      return entity_;
    }

    virtual size_t order() const override
    {
      return order_;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      function_->evaluate(entity_.geometry().global(xx), tmp_vector_);
      for (size_t ii = 0; ii < dimRangeRows; ++ii) {
        auto& retRow = ret[ii];
        for (size_t jj = 0; jj < dimRangeCols; ++jj) {
          retRow[jj] = tmp_vector_[ii * dimRangeRows + jj];
        }
      }
    } // ... evaluate(...)

    //    virtual void jacobian(const DomainType& /*xx*/, JacobianRangeType& ret) const override
    //    {
    //      assert(false);
    //    }

  private:
    const EntityType& entity_;
    const std::shared_ptr<const MathExpressionFunctionType> function_;
    const size_t order_;
    mutable FieldVector<RangeFieldType, dimRangeRows * dimRangeCols> tmp_vector_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = BaseType::rangeDimRows;
  static const unsigned int dimRangeCols = BaseType::rangeDimCols;
  typedef typename LocalfunctionType::RangeType RangeType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  //  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  //  {
  //    Dune::ParameterTree description;
  //    description["variable"] = "x";
  //    description["expression"] = "[x[0]; sin(x[0]); x[0]; x[0]]";
  //    description["order"] = "1";
  //    description["name"] = "function.expression";
  //    if (subName.empty())
  //      return description;
  //    else {
  //      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
  //      extendedDescription.add(description, subName);
  //      return extendedDescription;
  //    }
  //  } // ... defaultSettings(...)

  //  static ThisType* create(const DSC::ExtendedParameterTree settings)
  //  {
  //    // get necessary
  //    const std::string _variable = settings.get< std::string >("variable", "x");
  //    std::vector< std::string > _expressions;
  //    // lets see, if there is a key or vector
  //    if (settings.hasVector("expression")) {
  //      const std::vector< std::string > expr = settings.getVector< std::string >("expression", 1);
  //      for (auto vector: expr){
  //        _expressions.emplace_back(vector);
  //      }
  //    } else if (settings.hasKey("expression")) {
  //      const std::string expr = settings.get< std::string >("expression");
  //      _expressions.push_back(expr);
  //    } else
  //      DUNE_THROW(Dune::IOError,
  //                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
  //                 << " neither key nor vector 'expression' found in the following settings:\n"
  //                 << settings.reportString("  "));
  //    // get optional
  //    const size_t order = settings.get< size_t >("order", 0);
  //    const std::string name = settings.get< std::string >("name", "function.expression");
  //    // create and return
  //    return new ThisType(_variable, _expressions, order, name);
  //  } // ... create(...)

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

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::shared_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return std::shared_ptr<Localfunction>(new Localfunction(entity, function_, order_));
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
    static const unsigned int dimRangeRows = BaseType::dimRangeRows;
    static const unsigned int dimRangeCols = BaseType::dimRangeCols;
    static const unsigned int dimRange     = BaseType::dimRange;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const EntityType& entity, const std::shared_ptr<const MathExpressionFunctionType>& function,
                  std::vector<std::shared_ptr<const MathExpressionGradientType>> gradients, const size_t ord)
      : entity_(entity)
      , function_(function)
      , gradients_(gradients)
      , order_(ord)
      , global_point_(DomainFieldType(0))
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual const EntityType& entity() const
    {
      return entity_;
    }

    virtual size_t order() const override
    {
      return order_;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const override
    {
      assert(this->is_a_valid_point(xx));
      global_point_ = entity_.geometry().global(xx);
      function_->evaluate(global_point_, ret);
      assert(this_value_is_sane(ret));
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override
    {
      if (gradients_.size() == 0)
        DUNE_THROW(NotImplemented, "This function does not provide any gradients!");
      else {
        assert(gradients_.size() == dimRange);
        global_point_ = entity_.geometry().global(xx);
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

    const EntityType& entity_;
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
  static const unsigned int dimRangeRows = BaseType::dimRangeRows;
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

  virtual ThisType* copy() const override
  {
    return new ThisType(*this);
  }

  virtual std::string name() const override
  {
    return name_;
  }

  virtual std::shared_ptr<LocalfunctionType> local_function(const EntityType& entity) const override
  {
    return std::shared_ptr<Localfunction>(new Localfunction(entity, function_, gradients_, order_));
  }

private:
  void build_gradients(const std::string variable, const std::vector<std::vector<std::string>> gradient_expressions)
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

#endif // DUNE_STUFF_FUNCTIONS_EXPRESSION_HH
