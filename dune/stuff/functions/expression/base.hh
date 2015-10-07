// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber

#ifndef DUNE_STUFF_FUNCTION_EXPRESSION_BASE_HH
#define DUNE_STUFF_FUNCTION_EXPRESSION_BASE_HH

#include <sstream>
#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/common/color.hh>

#include "mathexpr.hh"

namespace Dune {
namespace Stuff {
namespace Functions {

/**
 *  \brief base class that makes a function out of the stuff from mathexpr.hh
 *  \attention  Most surely you do not want to use this class directly, but Functions::Expression!
 */
template <class DomainFieldImp, size_t domainDim, class RangeFieldImp, size_t rangeDim>
class MathExpressionBase
{
public:
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

  typedef DomainFieldImp DomainFieldType;
  static const size_t dimDomain = domainDim;

  typedef RangeFieldImp RangeFieldType;
  static const size_t dimRange = rangeDim;

  MathExpressionBase(const std::string _variable, const std::string _expression)
  {
    const std::vector<std::string> expressions(1, _expression);
    setup(_variable, expressions);
  }

  MathExpressionBase(const std::string _variable, const std::vector<std::string> _expressions)
  {
    setup(_variable, _expressions);
  }

  MathExpressionBase(const ThisType& _other)
  {
    setup(_other.variable(), _other.expression());
  }

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
  }

  ~MathExpressionBase()
  {
    cleanup();
  }

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
    std::lock_guard<std::mutex> guard(mutex_);
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
    std::lock_guard<std::mutex> guard(mutex_);
    // check for sizes
    assert(arg.size() > 0);
    if (ret.size() != dimRange)
      ret = Dune::DynamicVector<RangeFieldType>(dimRange);
    // copy arg
    for (size_t ii = 0; ii < std::min(domainDim, arg.size()); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::DynamicVector<RangeFieldType>::size_type ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  void evaluate(const Dune::FieldVector<DomainFieldType, dimDomain>& arg,
                Dune::DynamicVector<RangeFieldType>& ret) const
  {
    std::lock_guard<std::mutex> guard(mutex_);
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
    std::lock_guard<std::mutex> guard(mutex_);
    assert(arg.size() > 0);
    // copy arg
    for (size_t ii = 0; ii < std::min(dimDomain, arg.size()); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  void report(const std::string _name = "dune.stuff.function.mathexpressionbase", std::ostream& stream = std::cout,
              const std::string& _prefix = "") const
  {
    const std::string tmp = _name + "(" + variable() + ") = ";
    stream << _prefix << tmp;
    if (expression().size() == 1)
      stream << expression()[0] << std::endl;
    else {
      stream << "[ " << expression()[0] << std::endl;
      const std::string whitespace = Dune::Stuff::Common::whitespaceify(tmp + "[ ");
      for (size_t i = 1; i < expression().size() - 1; ++i)
        stream << _prefix << whitespace << expression()[i] << std::endl;
      stream << _prefix << whitespace << expression()[expression().size() - 1] << " ]" << std::endl;
    }
  } // void report(const std::string, std::ostream&, const std::string&) const

private:
  void setup(const std::string& _variable, const std::vector<std::string>& _expression)
  {
    static_assert((dimDomain > 0), "Really?");
    static_assert((dimRange > 0), "Really?");
    // set expressions
    if (_expression.size() < dimRange)
      DUNE_THROW(Dune::InvalidStateException,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " '_expression' too short (is "
                      << _expression.size()
                      << ", should be "
                      << dimRange
                      << ")!");
    for (size_t ii = 0; ii < dimRange; ++ii)
      expressions_.push_back(_expression[ii]);
    // set variable (i.e. "x")
    variable_ = _variable;
    // fill variables (i.e. "x[0]", "x[1]", ...)
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      std::stringstream variableStream;
      variableStream << variable_ << "[" << ii << "]";
      variables_.push_back(variableStream.str());
    }
    // create expressions
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      arg_[ii]      = new DomainFieldType(0.0);
      var_arg_[ii]  = new RVar(variables_[ii].c_str(), arg_[ii]);
      vararray_[ii] = var_arg_[ii];
    }
    for (size_t ii = 0; ii < dimRange; ++ii) {
      op_[ii] = new ROperation(expressions_[ii].c_str(), dimDomain, vararray_);
    }
  } // void setup(const std::string& _variable, const std::vector< std::string >& expressions)

  void cleanup()
  {
    for (size_t ii = 0; ii < dimRange; ++ii) {
      delete op_[ii];
    }
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      delete var_arg_[ii];
      delete arg_[ii];
    }
  } // void cleanup()

  std::string variable_;
  std::vector<std::string> variables_;
  std::vector<std::string> expressions_;
  size_t actualDimRange_;
  mutable DomainFieldType* arg_[dimDomain];
  RVar* var_arg_[dimDomain];
  RVar* vararray_[dimDomain];
  ROperation* op_[dimRange];
  mutable std::mutex mutex_;
}; // class MathExpressionBase

} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_EXPRESSION_BASE_HH
