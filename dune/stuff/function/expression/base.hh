#ifndef DUNE_STUFF_FUNCTION_EXPRESSION_BASE_HH
#define DUNE_STUFF_FUNCTION_EXPRESSION_BASE_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#else
#include "config.h"
#endif // ifdef HAVE_CMAKE_CONFIG

#include <sstream>
#include <vector>

//#if HAVE_EIGEN
//  #include <Eigen/Core>
//#endif // HAVE_EIGEN

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>

#include "mathexpr.hh"

namespace Dune {
namespace Stuff {
namespace Function {


// forward
template <class RangeFieldImp>
class Coefficient;


///**
//  \brief  Provides a function which evaluates a given mathematical expression at runtime.

//          Given a mathematical expression as a string, a domain \f$ K_d^{m \geq 1} \f$ and a range \f$ K_r^{n \geq 1}
//          \f$ this function represents the map
//          \f{eqnarray}
//            f:K_d^m \to K_r^n\\
//            x = (x_1, \dots, x_m)' \mapsto (f_1(x), \dots f_n(x))',
//          \f}
//          where \f$ K_d \f$ is the DomainType and \f$ K_r \f$ is the RangeType, usually a power of \f$ \mathcal{R}
//          \f$.
//          The name of the variable as well as the \f$ n \f$ expressions of \f$f_1, \dots, f_n\f$ have to be given in a
//          Dune::ParameterTree in the following form:
//\code variable: x
// expression.0: 2*x[0]
// expression.1: sin(x[1])*x[0]\endcode
//          There have to exist at least \f$n\f$ expressions; the entries of the variable are indexed by \f$[i]\f$ for
//          \f$ 0 \leq i \leq m - 1 \f$.
// **/
/**
 *  \brief base class that wraps makes a function out of the stuff from mathexpr.hh
 *  \attention  Most surely you do not want to use this class directly!
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class ExpressionBase
{
public:
  typedef ExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

  typedef DomainFieldImp DomainFieldType;
  static const int dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const int dimRange = rangeDim;

  ExpressionBase(const std::string _variable, const std::string _expression)
  {
    const std::vector<std::string> expressions(1, _expression);
    setup(_variable, expressions);
  } // NonparametricExpression(const std::string variable, const std::string expression)

  ExpressionBase(const std::string _variable, const std::vector<std::string> _expressions)
  {
    setup(_variable, _expressions);
  } // NonparametricExpression(const std::string variable, const std::vector< std::string >& expressions)

  ExpressionBase(const ThisType& _other)
  {
    setup(_other.variable(), _other.expression());
  } // NonparametricExpression(const ThisType& other)

  ThisType& operator=(const ThisType& _other)
  {
    if (this != &_other) {
      cleanup();
      variable_    = "";
      variables_   = std::vector<std::string>();
      expressions_ = std::vector<std::string>();
      setup(_other.variable(), _other.expression());
    }
    return this;
  } // ThisType& operator=(const ThisType& other)

  ~ExpressionBase()
  {
    cleanup();
  } // ~NonparametricExpression()

  std::string variable() const
  {
    return variable_;
  }

  const std::vector<std::string>& expression() const
  {
    return expressions_;
  }

  void evaluate(const Dune::FieldVector<DomainFieldType, dimDomain>& arg,
                Dune::FieldVector<RangeFieldType, dimRange>& ret) const
  {
    // copy arg
    for (typename Dune::FieldVector<DomainFieldType, dimDomain>::size_type ii = 0; ii < dimDomain; ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::FieldVector<RangeFieldType, dimRange>::size_type ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  /**
   *  \attention  arg will be used up to its size, ret will be resized!
   */
  void evaluate(const Dune::DynamicVector<DomainFieldType>& arg, Dune::DynamicVector<RangeFieldType>& ret) const
  {
    // check for sizes
    assert(arg.size() > 0);
    if (ret.size() != dimRange)
      ret = Dune::DynamicVector<RangeFieldType>(dimRange);
    // copy arg
    for (int ii = 0; ii < std::min(domainDim, int(arg.size())); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::DynamicVector<RangeFieldType>::size_type ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  void evaluate(const Dune::FieldVector<DomainFieldType, dimDomain>& arg,
                Dune::DynamicVector<RangeFieldType>& ret) const
  {
    // check for sizes
    if (ret.size() != dimRange)
      ret = Dune::DynamicVector<RangeFieldType>(dimRange);
    // copy arg
    for (typename Dune::FieldVector<DomainFieldType, dimDomain>::size_type ii = 0; ii < dimDomain; ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::DynamicVector<RangeFieldType>::size_type ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  /**
   *  \attention  arg will be used up to its size
   */
  void evaluate(const Dune::DynamicVector<DomainFieldType>& arg, Dune::FieldVector<RangeFieldType, dimRange>& ret) const
  {
    assert(arg.size() > 0);
    // copy arg
    for (int ii = 0; ii < std::min(domainDim, int(arg.size())); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::FieldVector<RangeFieldType, dimRange>::size_type ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  //#if HAVE_EIGEN
  //  void evaluate(const Eigen::Matrix< Eigen::Dynamic, 1, DomainFieldType >& arg,
  //                Eigen::Matrix< Eigen::Dynamic, 1, RangeFieldType >& ret)
  //  {
  //    // check for sizes
  //    assert(_arg.size() <= dimDomain);
  //    if (_ret.size != dimRange)
  //      ret.resize(dimRange);
  //    // copy arg
  //    for (typename Dune::DynamicVector< DomainFieldType >::size_type ii = 0; ii < dimDomain; ++ii)
  //      *(arg_[ii]) = arg[ii];
  //    // copy ret
  //    for (typename RangeType::size_type ii = 0; ii < ret.size(); ++ii)
  //      ret[ii] = op_[ii]->Val();
  //  }
  //#endif // HAVE_EIGEN

  void report(const std::string _name = "function.nonparametric.expression", std::ostream& stream = std::cout,
              const std::string& _prefix = "") const
  {
    const std::string tmp = _name + "(" + variable() + ") = ";
    stream << _prefix << tmp;
    if (expression().size() == 1)
      stream << expression()[0] << std::endl;
    else {
      stream << "[ " << expression()[0] << ";" << std::endl;
      const std::string whitespace = Dune::Stuff::Common::whitespaceify(tmp + "[ ");
      for (unsigned int i = 1; i < expression().size() - 1; ++i)
        stream << _prefix << whitespace << expression()[i] << ";" << std::endl;
      stream << _prefix << whitespace << expression()[expression().size() - 1] << " ]" << std::endl;
    }
  } // void report(const std::string, std::ostream&, const std::string&) const

private:
  friend class Coefficient<RangeFieldType>;

  void evaluate(const Dune::DynamicVector<DomainFieldType>& arg, RangeFieldType& ret) const
  {
    assert(dimRange == 1 && "I'm only here to be used by Function::Parametric::Coefficient, which has dimrange == 1");
    // copy arg
    for (int ii = 0; ii < std::min(domainDim, int(arg.size())); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    ret = op_[0]->Val();
  }

  void setup(const std::string& variable, const std::vector<std::string>& _expression)
  {
    dune_static_assert((dimDomain > 0), "Really?");
    dune_static_assert((dimRange > 0), "Really?");
    // set expressions
    if (_expression.size() < dimRange)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " '_expression' too short (is "
                      << _expression.size()
                      << ", should be "
                      << dimRange
                      << ")!");
    for (int ii = 0; ii < dimRange; ++ii)
      expressions_.push_back(_expression[ii]);
    // set variable (i.e. "x")
    variable_ = variable;
    // fill variables (i.e. "x[0]", "x[1]", ...)
    for (int ii = 0; ii < dimDomain; ++ii) {
      std::stringstream variableStream;
      variableStream << variable_ << "[" << ii << "]";
      variables_.push_back(variableStream.str());
    }
    // create epressions
    for (int ii = 0; ii < dimDomain; ++ii) {
      arg_[ii]      = new DomainFieldType(0.0);
      var_arg_[ii]  = new RVar(variables_[ii].c_str(), arg_[ii]);
      vararray_[ii] = var_arg_[ii];
    }
    for (int ii = 0; ii < dimRange; ++ii) {
      op_[ii] = new ROperation(expressions_[ii].c_str(), dimDomain, vararray_);
    }
  } // void setup(const std::string& variable, const std::vector< std::string >& expressions)

  void cleanup()
  {
    for (int ii = 0; ii < dimRange; ++ii) {
      delete op_[ii];
    }
    for (int ii = 0; ii < dimDomain; ++ii) {
      delete var_arg_[ii];
      delete arg_[ii];
    }
  } // void cleanup()

  std::string variable_;
  std::vector<std::string> variables_;
  std::vector<std::string> expressions_;
  unsigned int actualDimRange_;
  mutable DomainFieldType* arg_[dimDomain];
  RVar* var_arg_[dimDomain];
  RVar* vararray_[dimDomain];
  ROperation* op_[dimRange];
}; // class ExpressionBase

} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_EXPRESSION_BASE_HH
