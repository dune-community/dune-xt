#ifndef DUNE_STUFF_FUNCTION_EXPRESSION_HH
#define DUNE_STUFF_FUNCTION_EXPRESSION_HH

#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>

#include "expression/base.hh"
#include "interfaces.hh"

namespace Dune {
namespace Stuff {


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols = 1>
class FunctionExpression
    : public FunctionExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows * rangeDimCols>,
      public FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
{
  typedef FunctionExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows * rangeDimCols> BaseType;
  typedef FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> InterfaceType;

public:
  typedef FunctionExpression<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> ThisType;

  typedef typename InterfaceType::DomainFieldType DomainFieldType;
  static const int dimDomain = InterfaceType::dimDomain;
  typedef typename InterfaceType::DomainType DomainType;
  typedef typename InterfaceType::RangeFieldType RangeFieldType;
  static const int dimRange = InterfaceType::dimRange;
  typedef typename InterfaceType::RangeType RangeType;

  static const std::string id()
  {
    return InterfaceType::id() + ".expression";
  }

  FunctionExpression(const std::string _variable, const std::string _expression, const int orderIn = -1,
                     const std::string nameIn = id())
    : BaseType(_variable, _expression)
    , order_(orderIn)
    , name_(nameIn)
  {
  }

  FunctionExpression(const std::string _variable, const std::vector<std::string> _expressions, const int orderIn = -1,
                     const std::string nameIn = id())
    : BaseType(_variable, _expressions)
    , order_(orderIn)
    , name_(nameIn)
  {
  }

public:
  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["variable"]   = "x";
    description["expression"] = "[x[0]; sin(x[0]); x[0]; x[0]]";
    description["order"]      = "1";
    description["name"] = "function.expression";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  } // ... defaultSettings(...)

  static ThisType* create(const DSC::ExtendedParameterTree settings)
  {
    // get necessary
    const std::string _variable = settings.get<std::string>("variable", "x");
    std::vector<std::string> _expressions;
    // lets see, if there is a key or vector
    if (settings.hasVector("expression")) {
      const std::vector<std::string> expr = settings.getVector<std::string>("expression", 1);
      for (auto vector : expr) {
        _expressions.emplace_back(vector);
      }
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

  virtual int order() const
  {
    return order_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  virtual void evaluate(const DomainType& _x, RangeType& _ret) const
  {
    FieldVector<RangeFieldType, rangeDimRows * rangeDimCols> vector(0.0);
    BaseType::evaluate(_x, vector);
    for (int ii = 0; ii < rangeDimRows; ii++) {
      for (int jj = 0; jj < rangeDimCols; jj++) {
        _ret[ii][jj] = vector[ii * rangeDimRows + jj];
      }
    }
  }

  using InterfaceType::localFunction;

private:
  int order_;
  std::string name_;
}; // class FunctionExpression

template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class FunctionExpression<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public FunctionExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>,
      public FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef FunctionExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> InterfaceType;

public:
  typedef FunctionExpression<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

  typedef typename InterfaceType::DomainFieldType DomainFieldType;
  static const int dimDomain = InterfaceType::dimDomain;
  typedef typename InterfaceType::DomainType DomainType;
  typedef typename InterfaceType::RangeFieldType RangeFieldType;
  static const int dimRange = InterfaceType::dimRange;
  typedef typename InterfaceType::RangeType RangeType;

  static std::string static_id()
  {
    return InterfaceType::static_id() + ".expression";
  }

  FunctionExpression(const std::string _variable, const std::string _expression, const int orderIn = -1,
                     const std::string nameIn = static_id())
    : BaseType(_variable, _expression)
    , order_(orderIn)
    , name_(nameIn)
  {
  }

  FunctionExpression(const std::string _variable, const std::vector<std::string> _expressions, const int orderIn = -1,
                     const std::string nameIn = static_id())
    : BaseType(_variable, _expressions)
    , order_(orderIn)
    , name_(nameIn)
  {
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

  virtual int order() const
  {
    return order_;
  }

  virtual std::string name() const
  {
    return name_;
  }

  virtual void evaluate(const DomainType& _x, RangeType& _ret) const
  {
    BaseType::evaluate(_x, _ret);
  }

  using InterfaceType::localFunction;

private:
  int order_;
  std::string name_;
}; // class FunctionExpression< ..., 1 >


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_EXPRESSION_HH
