// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2018)
//   Kirsten Weber   (2013)
//   René Fritze     (2012 - 2016, 2018 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2014 - 2020)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DEFAULT_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_DEFAULT_HH

#include <limits>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/parallel/threadstorage.hh>

#include "base.hh"
#include <dune/xt/functions/interfaces/function.hh>

namespace Dune {
namespace XT {
namespace Functions {


/**
 * \note This is the matrix valued version, see below for scalar- and vectorvalued case.
 * \note We have a specialization for the other case on purpose, disabling all ctors and using helpers got too
 * complicated!
 */
template <size_t d, size_t r = 1, size_t rC = 1, class RangeField = double>
class ExpressionFunction : public FunctionInterface<d, r, rC, RangeField>
{
  using BaseType = FunctionInterface<d, r, rC, RangeField>;
  using ThisType = ExpressionFunction;

  using typename BaseType::DomainFieldType;
  using MathExpressionFunctionType = MathExpressionBase<DomainFieldType, d, RangeField, r * rC>;
  using MathExpressionGradientType = MathExpressionBase<DomainFieldType, d, RangeField, d>;

public:
  using typename BaseType::RangeReturnType;
  using DerivativeRangeReturnType = typename BaseType::DerivativeRangeReturnType;

  using DomainType = typename BaseType::DomainType;
  using RangeFieldType = typename BaseType::RangeFieldType;
  using SingleDerivativeRangeReturnType = typename BaseType::SingleDerivativeRangeReturnType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  static Common::Configuration defaults()
  {
    Common::Configuration config;
    config["variable"] = "x";
    config["expression"] = "[x[0] sin(x[0]*pi) exp(x[0]); x[0] sin(x[0]*pi) exp(x[0]); x[0] sin(x[0]*pi) exp(x[0])]";
    // seperation for range
    config["gradient.0"] = "[1 0 0; pi*cos(x[0]*pi) 0 0; exp(x[0]) 0 0]";
    config["gradient.1"] = "[1 0 0; pi*cos(x[0]*pi) 0 0; exp(x[0]) 0 0]";
    config["gradient.2"] = "[1 0 0; pi*cos(x[0]*pi) 0 0; exp(x[0]) 0 0]";
    config["order"] = "3";
    config["name"] = static_id();
    return config;
  } // ... defaults(...)

private:
  // This function is required because MathExpressionFunctionBase is only valid for a scalar or vector
  static Common::FieldVector<std::string, r * rC> matrix_to_vector(const Common::FieldMatrix<std::string, r, rC>& mat)
  {
    Common::FieldVector<std::string, r * rC> ret("");
    for (size_t rr = 0; rr < r; ++rr) {
      for (size_t cc = 0; cc < rC; ++cc)
        ret[rr * rC + cc] = mat[rr][cc];
    }
    return ret;
  }

public:
  /**
   * \brief Creates an Expression function
   *
   * \param variable variable of the Expression function, e.g. "x"
   * \param expressions FieldMatrix with size r x rC for the range of the Expression function
   *  range, e.g. [1 sin(x[0]); 2 x[1]] with r = rC = 2
   * \param gradient_expressions FieldVector<FieldMatrix>, vector of the jacobian matrices, where the vector has the
   size of r and the inner Matrix the size rC x d,
   e.g.
   *  [[0 0 ; cos(x[0]) 0] , [0 0 ; 0 1]] would be the gradient_expression corresponding to the expression above
   (if
   *  d = r = rC = 2)
   * Note: it is optional to provide a gradient if you do not want to use jacobian().
   * \param ord order of the Expression function
   * \param nm name of the Expression function
   *
   * For this example the correct constructor call is
   * FunctionType function("x", {{"1", "sin(x[0])"}, {"2", "x[1]"}}, {{{"0", "0"}, {"cos(x[0])", "0"}}, {{"0", "0"},
   {"0", "1"}}}, 3);
   */
  ExpressionFunction(const std::string& variable,
                     const Common::FieldMatrix<std::string, r, rC>& expressions,
                     const Common::FieldVector<Common::FieldMatrix<std::string, rC, d>, r>& gradient_expressions,
                     const size_t ord,
                     const std::string nm = static_id())
    : BaseType()
    , function_(variable, matrix_to_vector(expressions))
    , order_(ord)
    , name_(nm)
  {
    for (size_t cc = 0; cc < r; ++cc) {
      gradients_.emplace_back(std::vector<MathExpressionGradientType>());
      for (size_t rr = 0; rr < rC; ++rr) {
        const auto& gradient_expression = gradient_expressions[cc][rr];
        assert(gradient_expression.size() >= d);
        gradients_[cc].emplace_back(variable, gradient_expression);
      }
    }
  }

  /**
   * @brief ExpressionFunction without given gradient
   *
   */
  ExpressionFunction(const std::string& variable,
                     const Common::FieldMatrix<std::string, r, rC>& expressions,
                     const size_t ord,
                     const std::string nm = static_id())
    : BaseType()
    , function_(variable, matrix_to_vector(expressions))
    , order_(ord)
    , name_(nm)
  {}

  ExpressionFunction(const ThisType& other)
    : BaseType(other)
    , function_(other.function_)
    , order_(other.order_)
    , name_(other.name_)
  {
    for (size_t cc = 0; cc < other.gradients_.size(); ++cc) {
      gradients_.emplace_back(std::vector<MathExpressionGradientType>());
      for (size_t rr = 0; rr < other.gradients_[cc].size(); ++rr)
        gradients_[cc].emplace_back(other.gradients_[cc][rr]);
    }
  }

  ExpressionFunction(ThisType&&) = default;

private:
  ThisType* copy_as_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_function_impl());
  }
  std::string name() const override final
  {
    return name_;
  }

  int order(const Common::Parameter& /*param*/ = {}) const override final
  {
    return static_cast<int>(order_);
  }

  using BaseType::evaluate;

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& /*param*/ = {}) const override final
  {
    RangeReturnType ret(0.);
    Common::FieldVector<RangeFieldType, r * rC> tmp_vector_;
    function_.evaluate(point_in_global_coordinates, tmp_vector_);
    for (size_t rr = 0; rr < r; ++rr) {
      auto& retRow = ret[rr];
      for (size_t cc = 0; cc < rC; ++cc) {
        retRow[cc] = tmp_vector_[rr * rC + cc];
        check_value(point_in_global_coordinates, retRow);
      }
    }
    return ret;
  } // ... evaluate(...)

  using BaseType::jacobian;

  DerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                     const Common::Parameter& /*param*/ = {}) const override final
  {
    if (gradients_.size() != r)
      DUNE_THROW(NotImplemented, "Do not call jacobian() if no gradients are given on construction!");

    DerivativeRangeReturnType ret(Common::FieldMatrix<double, rC, d>(0.));
    for (size_t cc = 0; cc < r; ++cc) {
      assert(gradients_[cc].size() == rC);
      for (size_t rr = 0; rr < rC; ++rr) {
        gradients_[cc][rr].evaluate(point_in_global_coordinates, ret[cc][rr]);
        check_value(point_in_global_coordinates, ret[cc][rr]);
      }
    }
    return ret;
  } // ... jacobian(...)


  template <class V>
  void check_value(const DomainType& DXTC_DEBUG_ONLY(point_in_global_coordinates),
                   const V& DXTC_DEBUG_ONLY(value)) const
  {
#ifndef NDEBUG
#  ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
    size_t range = value.size();
    bool failure = false;
    std::string error_type;
    for (size_t rr = 0; rr < value.size(); ++rr) {
      if (Common::isnan(value[rr])) {
        failure = true;
        error_type = "NaN";
      } else if (Common::isinf(value[rr])) {
        failure = true;
        error_type = "inf";
      } else if (std::abs(value[rr]) > (0.9 * std::numeric_limits<double>::max())) {
        failure = true;
        error_type = "an unlikely value";
      }
      DUNE_THROW_IF(failure,
                    Common::Exceptions::internal_error,
                    "evaluating this function yielded "
                        << error_type << "!\n"
                        << "The variable of this function is:     " << function_.variable() << "\n"
                        << "The expression of this function is: " << function_.expression() // at
                        << "\n"
                        << "You tried to evaluate it with:   point_in_global_coordinates = "
                        << point_in_global_coordinates << "\n"
                        << "The result was:                       " << value[rr] << "\n\n"
                        << "You can disable this check by defining DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
    } // check_value(...)
#  endif // DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
#endif // NDEBUG
  }

private:
  MathExpressionFunctionType function_;
  size_t order_;
  std::string name_;
  std::vector<std::vector<MathExpressionGradientType>> gradients_;
}; // class ExpressionFunction


template <size_t d, size_t r, class RangeField>
class ExpressionFunction<d, r, 1, RangeField> : public FunctionInterface<d, r, 1, RangeField>
{
  using BaseType = FunctionInterface<d, r, 1, RangeField>;
  using ThisType = ExpressionFunction;

  using typename BaseType::DomainFieldType;
  using MathExpressionFunctionType = MathExpressionBase<DomainFieldType, d, RangeField, r>;
  using MathExpressionGradientType = MathExpressionBase<DomainFieldType, d, RangeField, d>;

public:
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::SingleDerivativeRangeReturnType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  static Common::Configuration defaults(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["type"] = static_id();
    config["variable"] = "x";
    config["expression"] = "[x[0] sin(pi*x[0]) exp(x[0])]";
    config["gradient"] = "[1 0 0; pi*cos(pi*x[0]) 0 0; exp(x[0]) 0 0]";
    config["order"] = "3";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... defaults(...)

  /**
   * \brief Creates an Expression function
   *
   * This constructor just expands expressions and gradient_expressions from a Common::FieldVector and
   * Common::FieldMatrix to ExpressionStringVectorType and GradientStringVectorType, respectively.
   *
   * For instance, if d = 2 and r = 1 and expression = {"x[0]*x[1]"}, gradient should be
   * {"x[1]","x[0]"}.
   * For r = 2 we have for example the expression = {"sin(x[0])", "x[1]*x[1]"} with the gradient
   * {{"cos(x[0])", "0"}, {"0", "2*x[1]"}}.
   */

  ExpressionFunction(const std::string& variable,
                     const Common::FieldVector<std::string, r>& expressions,
                     const Common::FieldMatrix<std::string, r, d>& gradient_expressions,
                     const size_t ord,
                     const std::string nm = static_id())
    : BaseType()
    , function_(variable, expressions)
    , order_(ord)
    , name_(nm)
  {
    for (size_t rr = 0; rr < r; ++rr)
      gradients_.emplace_back(new MathExpressionGradientType(variable, gradient_expressions[rr]));
  }

  /**
   * @brief ExpressionFunction without given gradient
   *
   */
  ExpressionFunction(const std::string& variable,
                     const Common::FieldVector<std::string, r>& expressions,
                     const size_t ord,
                     const std::string nm = static_id())
    : BaseType()
    , function_(variable, expressions)
    , order_(ord)
    , name_(nm)
  {}

  ExpressionFunction(const ThisType& other)
    : BaseType()
    , function_(other.function_)
    , order_(other.order_)
    , name_(other.name_)
  {
    for (size_t ii = 0; ii < other.gradients_.size(); ++ii)
      gradients_.emplace_back(new MathExpressionGradientType(*other.gradients_[ii]));
  }

  ExpressionFunction(ThisType&&) = default;

private:
  ThisType* copy_as_function_impl() const override
  {
    return new ThisType(*this);
  }

public:
  std::unique_ptr<ThisType> copy_as_function() const
  {
    return std::unique_ptr<ThisType>(this->copy_as_function_impl());
  }
  std::string name() const override final
  {
    return name_;
  }

  int order(const Common::Parameter& /*param*/ = {}) const override final
  {
    return static_cast<int>(order_);
  }

  using BaseType::evaluate;

  RangeReturnType evaluate(const DomainType& point_in_global_coordinates,
                           const Common::Parameter& /*param*/ = {}) const override final
  {
    RangeReturnType ret(0.);
    function_.evaluate(point_in_global_coordinates, ret);
    check_value(point_in_global_coordinates, ret);
    return ret;
  } // ... evaluate(...)

  using BaseType::jacobian;

  DerivativeRangeReturnType jacobian(const DomainType& point_in_global_coordinates,
                                     const Common::Parameter& /*param*/ = {}) const override final
  {
    if (gradients_.size() != r)
      DUNE_THROW(NotImplemented, "Do not call jacobian() if no gradients are given on construction!");
    DerivativeRangeReturnType ret(0.);
    for (size_t rr = 0; rr < r; ++rr) {
      gradients_[rr]->evaluate(point_in_global_coordinates, ret[rr]);
      check_value(point_in_global_coordinates, ret[rr]);
    }
    return ret;
  }

private:
  template <class V>
  void check_value(const DomainType& DXTC_DEBUG_ONLY(point_in_global_coordinates),
                   const V& DXTC_DEBUG_ONLY(value)) const
  {
#ifndef NDEBUG
#  ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
    size_t range = value.size();
    bool failure = false;
    std::string error_type;
    for (size_t rr = 0; rr < value.size(); ++rr) {
      if (Common::isnan(value[rr])) {
        failure = true;
        error_type = "NaN";
      } else if (Common::isinf(value[rr])) {
        failure = true;
        error_type = "inf";
      } else if (std::abs(value[rr]) > (0.9 * std::numeric_limits<double>::max())) {
        failure = true;
        error_type = "an unlikely value";
      }
      DUNE_THROW_IF(failure,
                    Common::Exceptions::internal_error,
                    "evaluating this function yielded "
                        << error_type << "!\n"
                        << "The variable of this function is:     " << function_.variable() << "\n"
                        << "The expression of this function is: " << function_.expression() << "\n"
                        << "You tried to evaluate it with:   point_in_global_coordinates = "
                        << point_in_global_coordinates << "\n"
                        << "The result was:                       " << value[rr] << "\n\n"
                        << "You can disable this check by defining DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
    } // check_value(...)
#  endif // DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
#endif // NDEBUG
  } // ... check_value(...)

private:
  MathExpressionFunctionType function_;
  size_t order_;
  std::string name_;
  std::vector<std::shared_ptr<MathExpressionGradientType>> gradients_;
}; // class ExpressionFunction


} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_EXPRESSION_DEFAULT_HH
