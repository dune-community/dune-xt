// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013, 2015 - 2016)
//   Tobias Leibner  (2014 - 2015, 2017)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_BASE_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_BASE_HH

#include <mutex>
#include <sstream>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/xt/common/color.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/string.hh>

#include "mathexpr.hh"

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE
#define DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE 64
#endif


namespace Dune {
namespace XT {
namespace Functions {


/**
 *  \brief      Base class that makes a function out of the stuff from mathexpr.hh
 *  \attention  Most surely you do not want to use this class directly, but Functions::ExpressionFunction!
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
      variable_ = "";
      variables_ = std::vector<std::string>();
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

  void report(const std::string _name = "function.mathexpressionbase",
              std::ostream& stream = std::cout,
              const std::string& _prefix = "") const
  {
    const std::string tmp = _name + "(" + variable() + ") = ";
    stream << _prefix << tmp;
    if (expression().size() == 1)
      stream << expression()[0] << std::endl;
    else {
      stream << "[ " << expression()[0] << std::endl;
      const std::string whitespace = Common::whitespaceify(tmp + "[ ");
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
                 "\n" << Common::color_string_red("ERROR:") << " '_expression' too short (is " << _expression.size()
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
      arg_[ii] = new DomainFieldType(0.0);
      var_arg_[ii] = new RVar(variables_[ii].c_str(), arg_[ii]);
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


/**
 *  \brief      Base class that makes a function out of the stuff from mathexpr.hh
 *  \attention  Most surely you do not want to use this class directly, but Functions::ParametricExpressionFunction!
 */
template <class D, class R, size_t r, size_t max_d = DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE>
class DynamicMathExpressionBase
{
public:
  typedef DynamicMathExpressionBase<D, R, r, max_d> ThisType;

  typedef D DomainFieldType;
  static const size_t maxDimDomain = max_d;

  typedef R RangeFieldType;
  static const size_t dimRange = r;

  DynamicMathExpressionBase(const std::vector<std::string>& vars, const std::vector<std::string>& exprs)
  {
    setup(vars, exprs);
  }

  DynamicMathExpressionBase(const ThisType& other)
  {
    setup(other.variables(), other.expressions());
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      cleanup();
      original_variables_ = other.original_variables_;
      variables_ = std::vector<std::string>();
      expressions_ = std::vector<std::string>();
      setup(other.original_variables_, other.expressions_);
    }
    return this;
  }

  ~DynamicMathExpressionBase()
  {
    cleanup();
  }

  const std::vector<std::string>& variables() const
  {
    return original_variables_;
  }

  const std::vector<std::string>& expressions() const
  {
    return expressions_;
  }

  void evaluate(const DynamicVector<DomainFieldType>& arg, FieldVector<RangeFieldType, dimRange>& ret) const
  {
    std::lock_guard<std::mutex> guard(mutex_);
    // check for sizes
    if (arg.size() != original_variables_.size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "arg.size(): " << arg.size() << "\n   "
                                << "variables.size(): "
                                << original_variables_.size());
    // copy arg
    for (size_t ii = 0; ii < original_variables_.size(); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (size_t ii = 0; ii < dimRange; ++ii)
      ret[ii] = op_[ii]->Val();
  }

private:
  void setup(const std::vector<std::string>& vars, const std::vector<std::string>& exprs)
  {
    static_assert((maxDimDomain > 0), "");
    static_assert((dimRange > 0), "");
    // set expressions
    if (exprs.size() != dimRange)
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "expressions.size(): " << exprs.size() << "\n   "
                                        << "dimRange: "
                                        << dimRange);
    expressions_ = exprs;
    for (const auto& ex : expressions_)
      if (ex.empty())
        DUNE_THROW(Common::Exceptions::wrong_input_given, "Given expressions must not be empty!");
    if (vars.size() > maxDimDomain)
      DUNE_THROW(
          Common::Exceptions::shapes_do_not_match,
          "This expression function of dynamic size was compiled to work for up to "
              << DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE
              << " variables, but you provided "
              << vars.size()
              << "!\n\n"
              << "Configure dune-xt-functions with a larger DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE!");
    original_variables_ = vars;
    for (const auto& var : original_variables_)
      if (var.empty())
        DUNE_THROW(Common::Exceptions::wrong_input_given, "Given variables must not be empty!");
    variables_ = original_variables_;
    for (size_t ii = variables_.size(); ii < maxDimDomain; ++ii)
      variables_.push_back("this_is_a_long_dummy_name_to_make_sure_it_is_not_used_" + Common::to_string(ii));
    assert(variables_.size() == maxDimDomain);
    // create expressions
    for (size_t ii = 0; ii < maxDimDomain; ++ii) {
      arg_[ii] = new DomainFieldType(0.0);
      var_arg_[ii] = new RVar(variables_[ii].c_str(), arg_[ii]);
      vararray_[ii] = var_arg_[ii];
    }
    for (size_t ii = 0; ii < dimRange; ++ii) {
      op_[ii] = new ROperation(expressions_[ii].c_str(), maxDimDomain, vararray_);
    }
  } // void setup(const std::string& _variable, const std::vector< std::string >& expressions)

  void cleanup()
  {
    for (size_t ii = 0; ii < dimRange; ++ii) {
      delete op_[ii];
    }
    for (size_t ii = 0; ii < maxDimDomain; ++ii) {
      delete var_arg_[ii];
      delete arg_[ii];
    }
  } // void cleanup()

  std::vector<std::string> original_variables_;
  std::vector<std::string> variables_;
  std::vector<std::string> expressions_;
  size_t actualDimRange_;
  mutable DomainFieldType* arg_[maxDimDomain];
  RVar* var_arg_[maxDimDomain];
  RVar* vararray_[maxDimDomain];
  ROperation* op_[dimRange];
  mutable std::mutex mutex_;
}; // class DynamicMathExpressionBase


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_EXPRESSION_BASE_HH
