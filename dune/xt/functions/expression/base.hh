// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2013 - 2017)
//   René Fritze     (2013, 2015 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014 - 2015, 2017, 2019 - 2020)

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
#  define DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE 64
#endif


namespace Dune::XT::Functions {


/**
 *  \brief      Base class that makes a function out of the stuff from mathexpr.hh
 *  \attention  Most surely you do not want to use this class directly, but Functions::ExpressionFunction!
 */
template <class DomainField, size_t domainDim, class RangeField, size_t rangeDim>
class MathExpressionBase
{
  static_assert((domainDim > 0), "Really?");
  static_assert((rangeDim > 0), "Really?");

public:
  using ThisType = MathExpressionBase;

  using DomainFieldType = DomainField;
  static const size_t domain_dim = domainDim;

  using RangeFieldType = RangeField;
  static const size_t range_dim = rangeDim;

  //  MathExpressionBase(const std::string var, const std::string expr)
  //  {
  //    const std::vector<std::string> expressions(1, expr);
  //    setup(var, expressions);
  //  }

  MathExpressionBase(const std::string& var, const Common::FieldVector<std::string, range_dim>& exprs)
    : variable_(var)
    , expressions_(exprs)
  {
    setup();
  }

  MathExpressionBase(const ThisType& other)
  {
    setup(other.variable(), other.expression());
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      cleanup();
      variable_ = std::string("");
      variables_ = Common::FieldVector<std::string, domain_dim>(std::string(""));
      expressions_ = Common::FieldVector<std::string, range_dim>(std::string(""));
      setup();
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

  const Common::FieldVector<std::string, range_dim>& expression() const
  {
    return expressions_;
  }

  void evaluate(const Dune::FieldVector<DomainFieldType, domain_dim>& arg,
                Dune::FieldVector<RangeFieldType, range_dim>& ret) const
  {
    std::lock_guard<std::mutex> guard(mutex_);
    // copy arg
    for (typename Dune::FieldVector<DomainFieldType, domain_dim>::size_type ii = 0; ii < domain_dim; ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::FieldVector<RangeFieldType, range_dim>::size_type ii = 0; ii < range_dim; ++ii)
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
    if (ret.size() != range_dim)
      ret = Dune::DynamicVector<RangeFieldType>(range_dim);
    // copy arg
    for (size_t ii = 0; ii < std::min(domain_dim, arg.size()); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::DynamicVector<RangeFieldType>::size_type ii = 0; ii < range_dim; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  void evaluate(const Dune::FieldVector<DomainFieldType, domain_dim>& arg,
                Dune::DynamicVector<RangeFieldType>& ret) const
  {
    std::lock_guard<std::mutex> guard(mutex_);
    // check for sizes
    if (ret.size() != range_dim)
      ret = Dune::DynamicVector<RangeFieldType>(range_dim);
    // copy arg
    for (typename Dune::FieldVector<DomainFieldType, domain_dim>::size_type ii = 0; ii < domain_dim; ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (typename Dune::DynamicVector<RangeFieldType>::size_type ii = 0; ii < range_dim; ++ii)
      ret[ii] = op_[ii]->Val();
  }

  /**
   *  \attention  arg will be used up to its size
   */
  void evaluate(const Dune::DynamicVector<DomainFieldType>& arg,
                Dune::FieldVector<RangeFieldType, range_dim>& ret) const
  {
    std::lock_guard<std::mutex> guard(mutex_);
    assert(arg.size() > 0);
    // copy arg
    for (size_t ii = 0; ii < std::min(domain_dim, arg.size()); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (size_t ii = 0; ii < range_dim; ++ii)
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
  void setup()
  {
    // fill variables (i.e. "x[0]", "x[1]", ...)
    for (size_t ii = 0; ii < domain_dim; ++ii)
      variables_[ii] = variable_ + "[" + Common::to_string(ii) + "]";
    // create expressions
    for (size_t ii = 0; ii < domain_dim; ++ii) {
      arg_[ii] = new DomainFieldType(0.0);
      var_arg_[ii] = new RVar(variables_[ii].c_str(), arg_[ii]);
      vararray_[ii] = var_arg_[ii];
    }
    for (size_t ii = 0; ii < range_dim; ++ii) {
      op_[ii] = new ROperation(expressions_[ii].c_str(), domain_dim, vararray_);
    }
  } // ... setup(...)

  void cleanup()
  {
    for (size_t ii = 0; ii < range_dim; ++ii) {
      delete op_[ii];
    }
    for (size_t ii = 0; ii < domain_dim; ++ii) {
      delete var_arg_[ii];
      delete arg_[ii];
    }
  } // void cleanup()

  std::string variable_;
  Common::FieldVector<std::string, domain_dim> variables_;
  Common::FieldVector<std::string, range_dim> expressions_;
  size_t actualDimRange_;
  mutable DomainFieldType* arg_[domain_dim];
  RVar* var_arg_[domain_dim];
  RVar* vararray_[domain_dim];
  ROperation* op_[range_dim];
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
  using ThisType = DynamicMathExpressionBase;

  using DomainFieldType = D;
  static const size_t maxDimDomain = max_d;

  using RangeFieldType = R;
  static const size_t range_dim = r;

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
      originalvars_ = other.originalvars_;
      variables_ = std::vector<std::string>();
      expressions_ = std::vector<std::string>();
      setup(other.originalvars_, other.expressions_);
    }
    return this;
  }

  ~DynamicMathExpressionBase()
  {
    cleanup();
  }

  const std::vector<std::string>& variables() const
  {
    return originalvars_;
  }

  const std::vector<std::string>& expressions() const
  {
    return expressions_;
  }

  void evaluate(const DynamicVector<DomainFieldType>& arg, FieldVector<RangeFieldType, range_dim>& ret) const
  {
    std::lock_guard<std::mutex> guard(mutex_);
    // check for sizes
    if (arg.size() != originalvars_.size())
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "arg.size(): " << arg.size() << "\n   "
                                << "variables.size(): " << originalvars_.size());
    // copy arg
    for (size_t ii = 0; ii < originalvars_.size(); ++ii)
      *(arg_[ii]) = arg[ii];
    // copy ret
    for (size_t ii = 0; ii < range_dim; ++ii)
      ret[ii] = op_[ii]->Val();
  }

private:
  void setup(const std::vector<std::string>& vars, const std::vector<std::string>& exprs)
  {
    static_assert((maxDimDomain > 0), "");
    static_assert((range_dim > 0), "");
    // set expressions
    if (exprs.size() != range_dim)
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "expressions.size(): " << exprs.size() << "\n   "
                                        << "range_dim: " << range_dim);
    expressions_ = exprs;
    for (const auto& ex : expressions_)
      if (ex.empty())
        DUNE_THROW(Common::Exceptions::wrong_input_given, "Given expressions must not be empty!");
    if (vars.size() > maxDimDomain)
      DUNE_THROW(Common::Exceptions::shapes_do_not_match,
                 "This expression function of dynamic size was compiled to work for up to "
                     << DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE << " variables, but you provided "
                     << vars.size() << "!\n\n"
                     << "Configure dune-xt with a larger DUNE_XT_FUNCTIONS_EXPRESSION_BASE_MAX_DYNAMIC_SIZE!");
    originalvars_ = vars;
    for (const auto& var : originalvars_)
      if (var.empty())
        DUNE_THROW(Common::Exceptions::wrong_input_given, "Given variables must not be empty!");
    variables_ = originalvars_;
    for (size_t ii = variables_.size(); ii < maxDimDomain; ++ii)
      variables_.push_back("this_is_a_long_dummy_name_to_make_sure_it_is_not_used_" + Common::to_string(ii));
    assert(variables_.size() == maxDimDomain);
    // create expressions
    for (size_t ii = 0; ii < maxDimDomain; ++ii) {
      arg_[ii] = new DomainFieldType(0.0);
      var_arg_[ii] = new RVar(variables_[ii].c_str(), arg_[ii]);
      vararray_[ii] = var_arg_[ii];
    }
    for (size_t ii = 0; ii < range_dim; ++ii) {
      op_[ii] = new ROperation(expressions_[ii].c_str(), maxDimDomain, vararray_);
    }
  } // void setup(const std::string& var, const std::vector< std::string >& expressions)

  void cleanup()
  {
    for (size_t ii = 0; ii < range_dim; ++ii) {
      delete op_[ii];
    }
    for (size_t ii = 0; ii < maxDimDomain; ++ii) {
      delete var_arg_[ii];
      delete arg_[ii];
    }
  } // void cleanup()

  std::vector<std::string> originalvars_;
  std::vector<std::string> variables_;
  std::vector<std::string> expressions_;
  size_t actualDimRange_;
  mutable DomainFieldType* arg_[maxDimDomain];
  RVar* var_arg_[maxDimDomain];
  RVar* vararray_[maxDimDomain];
  ROperation* op_[range_dim];
  mutable std::mutex mutex_;
}; // class DynamicMathExpressionBase


} // namespace Dune::XT::Functions

#endif // DUNE_XT_FUNCTIONS_EXPRESSION_BASE_HH
