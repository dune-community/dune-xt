#ifndef DUNE_STUFF_FUNCTION_EXPRESSION_HH
#define DUNE_STUFF_FUNCTION_EXPRESSION_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>

#include "expression/base.hh"
#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace Function {

template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Expression : public ExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
#if HAVE_DUNE_FEM
                   ,
                   public Dune::Fem::Function<Dune::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
                                              Expression<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>>
#endif
                   ,
                   public Interface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim>
{
public:
  typedef Expression<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;
  typedef ExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> BaseType;
  typedef Interface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> InterfaceType;

  using typename InterfaceType::DomainFieldType;
  using InterfaceType::dimDomain;
  typedef typename InterfaceType::DomainType DomainType;
  using typename InterfaceType::RangeFieldType;
  using InterfaceType::dimRange;
  typedef typename InterfaceType::RangeType RangeType;

  static const std::string id()
  {
    return InterfaceType::id() + ".expression";
  }

  Expression(const std::string _variable, const std::string _expression, const int _order = -1,
             const std::string _name = "function.expression")
    : BaseType(_variable, _expression)
    , order_(_order)
    , name_(_name)
  {
  }

  Expression(const std::string _variable, const std::vector<std::string> _expressions, const int _order = -1,
             const std::string _name = "function.expression")
    : BaseType(_variable, _expressions)
    , order_(_order)
    , name_(_name)
  {
  }

  Expression(const ThisType& other)
    : BaseType(other)
    , order_(other.order())
    , name_(other.name())
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      BaseType::operator=(other);
      order_            = other.order();
      name_             = other.name();
    }
    return this;
  } // ThisType& operator=(const ThisType& other)

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["variable"]   = "x";
    description["expression"] = "[x[0]; sin(x[0])]";
    description["order"]      = "1";
    description["name"] = "function.expression";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  }

  static ThisType* createFromDescription(const DSC::ExtendedParameterTree description)
  {
    // get necessary
    const std::string _variable = description.get<std::string>("variable", "x");
    std::vector<std::string> _expressions;
    // lets see, if there is a key or vector
    if (description.hasVector("expression")) {
      const std::vector<std::string> expr = description.getVector<std::string>("expression", 1);
      for (size_t ii = 0; ii < expr.size(); ++ii)
        _expressions.push_back(expr[ii]);
    } else if (description.hasKey("expression")) {
      const std::string expr = description.get<std::string>("expression");
      _expressions.push_back(expr);
    } else
      DUNE_THROW(Dune::IOError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                      << " neither key nor vector 'expression' found in the following description:\n"
                      << description.reportString("  "));
    // get optional
    const int _order        = description.get<int>("order", -1);
    const std::string _name = description.get<std::string>("name", "function.expression");
    // create and return
    return new ThisType(_variable, _expressions, _order, _name);
  } // static ThisType createFromDescription(const Dune::ParameterTree& _description)

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

private:
  int order_;
  std::string name_;
}; // class Expression


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_EXPRESSION_HH
