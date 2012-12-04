#ifndef DUNE_STUFF_FUNCTION_EXPRESSION_HH
#define DUNE_STUFF_FUNCTION_EXPRESSION_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <sstream>
#include <vector>

#ifdef HAVE_EIGEN
#include <Eigen/Core>
#endif // HAVE_EIGEN

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>

#ifdef HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif // HAVE_DUNE_FEM

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>

#include "expression/mathexpr.hh"
#include "interface.hh"

namespace Dune {
namespace Stuff {
namespace Function {

/**
  \brief  Provides a function which evaluates a given mathematical expression at runtime.

          Given a mathematical expression as a string, a domain \f$ K_d^{m \geq 1} \f$ and a range \f$ K_r^{n \geq 1}
          \f$ this function represents the map
          \f{eqnarray}
            f:K_d^m \to K_r^n\\
            x = (x_1, \dots, x_m)' \mapsto (f_1(x), \dots f_n(x))',
          \f}
          where \f$ K_d \f$ is the DomainType and \f$ K_r \f$ is the RangeType, usually a power of \f$ \mathcal{R} \f$.
          The name of the variable as well as the \f$ n \f$ expressions of \f$f_1, \dots, f_n\f$ have to be given in a
          Dune::ParameterTree in the following form:
\code variable: x
expression.0: 2*x[0]
expression.1: sin(x[1])*x[0]\endcode
          There have to exist at least \f$n\f$ expressions; the entries of the variable are indexed by \f$[i]\f$ for
          \f$ 0 \leq i \leq m - 1 \f$.
 **/
template <class DomainFieldImp, int maxDimDomain, class RangeFieldImp, int maxDimRange>
class Expression : public Interface<DomainFieldImp, maxDimDomain, RangeFieldImp, maxDimRange>
{
public:
  typedef DomainFieldImp DomainFieldType;

  typedef RangeFieldImp RangeFieldType;

  typedef Interface<DomainFieldImp, maxDimDomain, RangeFieldImp, maxDimRange> BaseType;

  typedef Expression<DomainFieldImp, maxDimDomain, RangeFieldImp, maxDimRange> ThisType;

  Expression(const std::string _variable, const std::string _expression)
  {
    const std::vector<std::string> expressions(1, _expression);
    setup(_variable, expressions);
  } // Expression(const std::string variable, const std::string expression)

  Expression(const std::string _variable, const std::vector<std::string> _expressions)
  {
    setup(_variable, _expressions);
  } // Expression(const std::string variable, const std::vector< std::string >& expressions)

  Expression(const ThisType& other)
  {
    setup(other.variable(), other.expression());
  } // Expression(const ThisType& other)

  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree)
  {
    const Dune::Stuff::Common::ExtendedParameterTree extendedParamtree(paramTree);
    // get variable
    if (!extendedParamtree.hasKey("variable"))
      DUNE_THROW(Dune::RangeError, "\nError: missing key 'variable'!");
    const std::string variable = extendedParamtree.get("variable", "meaningless_default_value");
    // get expressions
    if (!extendedParamtree.hasKey("expression"))
      DUNE_THROW(Dune::RangeError, "\nError: missing key or vector 'expression'!");
    const std::vector<std::string> expressions =
        extendedParamtree.getVector<std::string>("expression", "meaningless_default_value", 0);
    // create and return
    return ThisType(variable, expressions);
  } // static ThisType createFromParamTree(const Stuff::Common::ExtendedParameterTree& paramTree)

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      cleanup();
      setup(other.variable(), other.expression());
    }
    return this;
  } // ThisType& operator=(const ThisType& other)

  ~Expression()
  {
    cleanup();
  } // ~Expression()

  void report(const std::string name = "stuff.function.expression", std::ostream& stream = std::cout,
              const std::string& prefix = "") const
  {
    const std::string tmp = name + "(" + variable() + ") = ";
    stream << prefix << tmp;
    if (expression().size() == 1)
      stream << expression()[0] << std::endl;
    else {
      stream << "[ " << expression()[0] << ";" << std::endl;
      const std::string whitespace = Dune::Stuff::Common::whitespaceify(tmp + "[ ");
      for (unsigned int i = 1; i < expression().size() - 1; ++i)
        stream << prefix << whitespace << expression()[i] << ";" << std::endl;
      stream << prefix << whitespace << expression()[expression().size() - 1] << " ]" << std::endl;
    }
  } // void report(const std::string, std::ostream&, const std::string&) const

  std::string variable() const
  {
    return variable_;
  }

  const std::vector<std::string> expression() const
  {
    return expressions_;
  }

  unsigned int dimRange() const
  {
    return std::min(int(actualDimRange_), maxDimRange);
  }

  //! needed for Interface
  virtual void evaluate(const Dune::FieldVector<DomainFieldImp, maxDimDomain>& arg,
                        Dune::FieldVector<RangeFieldImp, maxDimRange>& ret) const
  {
    // ensure right dimensions
    assert(arg.size() <= maxDimDomain);
    assert(ret.size() <= dimRange());
    // arg
    for (typename Dune::FieldVector<DomainFieldImp, maxDimDomain>::size_type i = 0; i < arg.size(); ++i) {
      *(arg_[i]) = arg[i];
    }
    // ret
    for (typename Dune::FieldVector<RangeFieldImp, maxDimRange>::size_type i = 0; i < ret.size(); ++i) {
      ret[i] = op_[i]->Val();
    }
  }

  template <class DomainVectorType, class RangeVectorType>
  void evaluate(const Dune::DenseVector<DomainVectorType>& arg, Dune::DenseVector<RangeVectorType>& ret) const
  {
    // ensure right dimensions
    assert(arg.size() <= maxDimDomain);
    assert(ret.size() <= dimRange());
    // arg
    for (typename Dune::DenseVector<DomainVectorType>::size_type i = 0; i < arg.size(); ++i) {
      *(arg_[i]) = arg[i];
    }
    // ret
    for (typename Dune::DenseVector<RangeVectorType>::size_type i = 0; i < ret.size(); ++i) {
      ret[i] = op_[i]->Val();
    }
  }

#ifdef HAVE_EIGEN
  /**
   * \attention ret is resized to size dimRange()!
   */
  void evaluate(const Eigen::VectorXd& arg, Eigen::VectorXd& ret) const
  {
    // ensure right dimensions
    assert(arg.size() <= maxDimDomain);
    ret.resize(dimRange());
    // arg
    for (int i = 0; i < arg.size(); ++i) {
      *(arg_[i]) = arg(i);
    }
    // ret
    for (int i = 0; i < ret.size(); ++i) {
      ret(i) = op_[i]->Val();
    }
  } // void evaluate(const Eigen::VectorXd& arg, Eigen::VectorXd& ret) const
#endif // HAVE_EIGEN

private:
  void setup(const std::string& _variable, const std::vector<std::string>& _expressions)
  {
    assert(maxDimDomain > 0);
    assert(maxDimRange > 0);
    // set expressions
    if (_expressions.size() < 1)
      DUNE_THROW(Dune::InvalidStateException, "\nError: Given 'expressions'-vector is empty!");
    actualDimRange_ = std::min(int(_expressions.size()), maxDimRange);
    expressions_    = _expressions;
    // set variable (i.e. "x")
    variable_ = _variable;
    // fill variables (i.e. "x[0]", "x[1]", ...)
    for (int i = 0; i < maxDimDomain; ++i) {
      std::stringstream variableStream;
      variableStream << variable_ << "[" << i << "]";
      variables_.push_back(variableStream.str());
    }
    // create epressions
    for (unsigned int i = 0; i < maxDimDomain; ++i) {
      arg_[i]      = new DomainFieldType(0.0);
      var_arg_[i]  = new RVar(variables_[i].c_str(), arg_[i]);
      vararray_[i] = var_arg_[i];
    }
    for (unsigned int i = 0; i < dimRange(); ++i) {
      op_[i] = new ROperation(expressions_[i].c_str(), maxDimDomain, vararray_);
    }
  } // void setup(const std::string& variable, const std::vector< std::string >& expressions)

  void cleanup()
  {
    for (unsigned int i = 0; i < dimRange(); ++i) {
      delete op_[i];
    }
    for (unsigned int i = 0; i < maxDimDomain; ++i) {
      delete var_arg_[i];
      delete arg_[i];
    }
  } // void cleanup()

  std::string variable_;
  std::vector<std::string> variables_;
  std::vector<std::string> expressions_;
  unsigned int actualDimRange_;
  mutable DomainFieldType* arg_[maxDimDomain];
  RVar* var_arg_[maxDimDomain];
  RVar* vararray_[maxDimDomain];
  ROperation* op_[maxDimRange];
}; // class Expression

} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_EXPRESSION_HH
