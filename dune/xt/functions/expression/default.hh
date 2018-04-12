// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2012 - 2016, 2018)
//   Tobias Leibner  (2014 - 2018)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DEFAULT_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_DEFAULT_HH

#include <limits>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/parallel/threadstorage.hh>

#include "base.hh"
#include "../interfaces.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1>
class ExpressionFunction
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> BaseType;
  typedef ExpressionFunction<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols> ThisType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, rangeDim * rangeDimCols>
      MathExpressionFunctionType;
  typedef MathExpressionBase<DomainFieldImp, domainDim, RangeFieldImp, domainDim> MathExpressionGradientType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeFieldType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;
  typedef typename std::vector<std::vector<std::string>> ExpressionStringVectorType;
  typedef typename std::vector<std::vector<std::vector<std::string>>> GradientStringVectorType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".expression";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["variable"] = "x";
    config["expression"] = "[x[0] sin(x[0]) exp(x[0]); x[0] sin(x[0]) exp(x[0]); x[0] sin(x[0]) exp(x[0])]";
    config["order"] = "3";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // get expression
    ExpressionStringVectorType expression_as_vectors;
    // try to get expression as FieldVector (if dimRangeCols == 1) or as FieldMatrix (else)
    try {
      get_helper<>::get_expression(cfg, expression_as_vectors);
    } catch (std::exception) {
      // if dimRangeCols == 1 and we could not get expression as FieldVector, get it as FieldMatrix with one col
      if (dimRangeCols == 1) { // the 2 here is on purpose, anything > 1 will suffice
        get_helper<2>::get_expression(cfg, expression_as_vectors);
      } else
        std::rethrow_exception(std::current_exception());
    }
    // get gradient
    GradientStringVectorType gradient_as_vectors;
    if (cfg.has_key("gradient")) {
      assert(dimRangeCols == 1
             && "Use gradient.0, gradient.1, ... for the gradient of the first, second, ... row, respectively!");
      get_gradient(cfg, gradient_as_vectors, "gradient");
    } else if (cfg.has_key("gradient.0")) {
      get_gradient(cfg, gradient_as_vectors, "gradient.0");
    }
    // create
    return Common::make_unique<ThisType>(cfg.get("variable", default_cfg.get<std::string>("variable")),
                                         expression_as_vectors,
                                         cfg.get("order", default_cfg.get<size_t>("order")),
                                         cfg.get("name", default_cfg.get<std::string>("name")),
                                         gradient_as_vectors);
  } // ... create(...)

  /**
   * \brief Creates an Expression function where every component is identical
   *
   * For example, if dimDomain = dimRange = dimRangeCols = 2 and expression = "x[0]*x[1]", gradient should be
   * ["x[1]" "x[0]"]. Then the resulting function is [x[0]*x[1] x[0]*x[1]; x[0]*x[1] x[0]*x[1]] and the gradient is
   * [[x[1] x[0]; x[1] x[0]] [x[1] x[0] x[1] x[0]].
   */
  ExpressionFunction(const std::string variable,
                     const std::string expression,
                     const size_t ord = 0,
                     const std::string nm = static_id(),
                     const std::vector<std::string> gradient = std::vector<std::string>())
    : order_(ord)
    , name_(nm)
  {
    // create ExpressionStringVectorType with identical expressions
    const std::vector<std::string> expression_row(dimRangeCols, expression);
    const ExpressionStringVectorType expressions(dimRange, expression_row);
    // create associated gradient vector
    GradientStringVectorType gradient_expressions;
    if (gradient.size() > 0) {
      const std::vector<std::vector<std::string>> gradient_row(dimRange, gradient);
      for (size_t cc = 0; cc < dimRangeCols; ++cc) {
        gradient_expressions.emplace_back(gradient_row);
      }
    }
    // build function and gradient
    build_function(variable, expressions);
    build_gradients(variable, gradient_expressions);
  }

  /**
   * \brief Creates an Expression function with dimRangeCols = 1
   *
   * This constructor just expands expressions and gradient_expressions from a std::vector< std::string > and
   * std::vector< std::vector< std::string > to ExpressionStringVectorType and GradientStringVectorType, respectively.
   */
  template <typename std_str = std::string,
            typename = typename std::enable_if<(dimRangeCols == 1) && sizeof(std_str)>::type>
  ExpressionFunction(
      const std_str variable,
      const std::vector<std::string> expressions,
      const size_t ord = default_config().template get<size_t>("order"),
      const std::string nm = static_id(),
      const std::vector<std::vector<std::string>> gradient_expressions = std::vector<std::vector<std::string>>())
    : function_(new MathExpressionFunctionType(variable, expressions))
    , order_(ord)
    , name_(nm)
  {
    static_assert(dimRangeCols == 1, "This constructor does not make sense for dimRangeCols > 1!");
    GradientStringVectorType gradient_expressions_vec;
    if (gradient_expressions.size() > 0) {
      gradient_expressions_vec.emplace_back(gradient_expressions);
    }
    build_gradients(variable, gradient_expressions_vec);
  }

  /**
   * \brief Creates an Expression function
   *
   * \param variable variable of the Expression function, e.g. "x"
   * \param expressions vector< vector< string > >, where the inner vectors are the rows of the Expression functions
   *  range, e.g. [[1 sin(x[0])] [2*x[0] x[1]]] gives the range [1 sin(x[0]); 2 x[1]]
   * \param ord order of the Expression function
   * \param nm name of the Expression function
   * \param gradient_expressions vector< vector< vector< string > > >, vector of the jacobian matrices (written as
   *  vector< vector< string > >, where the inner vectors are the rows) of the columns of the Expression function, e.g.
   *  [[[0 0] [2 0]] [[cos(x[0]) 0] [0 1]]] would be the gradient_expression corresponding to the expression above (if
   *  dimDomain = dimRange = dimRangeCols = 2)
   */
  ExpressionFunction(const std::string variable,
                     const ExpressionStringVectorType expressions,
                     const size_t ord = 0,
                     const std::string nm = static_id(),
                     const GradientStringVectorType gradient_expressions = GradientStringVectorType())
    : order_(ord)
    , name_(nm)
  {
    build_function(variable, expressions);
    build_gradients(variable, gradient_expressions);
  }

#if !DUNE_XT_WITH_PYTHON_BINDINGS
  ExpressionFunction(const ThisType& other) = default;
#endif

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      function_ = other.function_;
      order_ = other.order_;
      name_ = other.name_;
      gradients_ = other.gradients_;
    }
    return *this;
  }

  std::string type() const override final
  {
    return BaseType::static_id() + ".expression";
  }

  std::string name() const override final
  {
    return name_;
  }

  virtual size_t order(const Common::Parameter& /*mu*/ = {}) const override final
  {
    return order_;
  }

private:
  template <bool is_not_tensor = (rangeDimCols == 1), bool anything = true>
  struct eval_helper
  {
    static void
    evaluate(const std::shared_ptr<const MathExpressionFunctionType>& func,
             typename Common::PerThreadValue<FieldVector<RangeFieldType, dimRange * dimRangeCols>>& /*tmp_vector*/,
             const DomainType& xx,
             RangeType& ret)
    {
      func->evaluate(xx, ret);
    }

    static void jacobian(const std::vector<std::vector<std::shared_ptr<const MathExpressionGradientType>>>& gradients,
                         const DomainType& xx,
                         JacobianRangeType& ret)
    {
      assert(gradients[0].size() == dimRange);
      for (size_t rr = 0; rr < dimRange; ++rr) {
        gradients[0][rr]->evaluate(xx, ret[rr]);
      }
    }
  }; // struct evap_helper<true, ...>

  template <bool anything>
  struct eval_helper<false, anything>
  {
    static void
    evaluate(const std::shared_ptr<const MathExpressionFunctionType>& func,
             typename Common::PerThreadValue<FieldVector<RangeFieldType, dimRange * dimRangeCols>>& tmp_vector,
             const DomainType& xx,
             RangeType& ret)
    {
      func->evaluate(xx, *tmp_vector);
      for (size_t rr = 0; rr < dimRange; ++rr) {
        auto& retRow = ret[rr];
        for (size_t cc = 0; cc < dimRangeCols; ++cc)
          retRow[cc] = (*tmp_vector)[rr * dimRangeCols + cc];
      }
    } // ... evaluate(...)

    static void jacobian(const std::vector<std::vector<std::shared_ptr<const MathExpressionGradientType>>>& gradients,
                         const DomainType& xx,
                         JacobianRangeType& ret)
    {
      for (size_t cc = 0; cc < dimRangeCols; ++cc) {
        assert(gradients[cc].size() == dimRange);
        for (size_t rr = 0; rr < dimRange; ++rr) {
          gradients[cc][rr]->evaluate(xx, ret[cc][rr]);
        }
      }
    } // ... jacobian(...)
  }; // struct eval_helper<false, ...>

public:
  using BaseType::evaluate;

  void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    eval_helper<>::evaluate(function_, tmp_vector_, xx, ret);
#ifndef NDEBUG
#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
    bool failure = false;
    std::string error_type;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      *tmp_row_ = ret[rr];
      for (size_t cc = 0; cc < dimRangeCols; ++cc) {
        if (Dune::XT::Common::isnan(tmp_row_->operator[](cc))) {
          failure = true;
          error_type = "NaN";
        } else if (Dune::XT::Common::isinf(tmp_row_->operator[](cc))) {
          failure = true;
          error_type = "inf";
        } else if (std::abs(tmp_row_->operator[](cc)) > (0.9 * std::numeric_limits<double>::max())) {
          failure = true;
          error_type = "an unlikely value";
        }
        if (failure)
          DUNE_THROW(Common::Exceptions::internal_error,
                     "evaluating this function yielded "
                         << error_type
                         << "!\n"
                         << "The variable of this function is:     "
                         << function_->variable()
                         << "\n"
                         << "The expression of this function is: "
                         << function_->expression().at(rr * dimRangeCols + cc)
                         << "\n"
                         << "You tried to evaluate it with:   xx = "
                         << xx
                         << "\n"
                         << "The result was:                       "
                         << tmp_row_->operator[](cc)
                         << "\n\n"
                         << "You can disable this check by defining DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS\n");
      }
    }
#endif // DUNE_XT_FUNCTIONS_EXPRESSION_DISABLE_CHECKS
#endif // NDEBUG
  } // ... evaluate(...)

  using BaseType::jacobian;

  void jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    if (gradients_.size() == 0) {
      DUNE_THROW(NotImplemented, "This function does not provide any gradients!");
    } else {
      assert(gradients_.size() == dimRangeCols);
      eval_helper<>::jacobian(gradients_, xx, ret);
    }
  } // ... jacobian(...)

private:
  // fill the rows of the dimRange x dimRangeCols matrix (aka vector< vector< string > > expression) in a vector of
  // length dimRange*dimRangeCols, e.g. [3 4; 1 2] becomes [3 4 1 2], in order to create function_
  void build_function(const std::string variable, const ExpressionStringVectorType& expressions)
  {
    assert(expressions.size() >= dimRange);
    std::vector<std::string> reordered_expressions;
    for (size_t rr = 0; rr < dimRange; ++rr) {
      assert(expressions[rr].size() >= dimRangeCols);
      for (size_t cc = 0; cc < dimRangeCols; ++cc) {
        reordered_expressions.emplace_back(expressions[rr][cc]);
      }
    }
    function_ = std::make_shared<MathExpressionFunctionType>(variable, reordered_expressions);
  } // ... build_function(...)

  void build_gradients(const std::string variable, const GradientStringVectorType& gradient_expressions)
  {
    assert(gradient_expressions.size() == 0 || gradient_expressions.size() >= dimRangeCols);
    if (gradient_expressions.size() > 0) {
      for (size_t cc = 0; cc < dimRangeCols; ++cc) {
        gradients_.emplace_back(std::vector<std::shared_ptr<const MathExpressionGradientType>>());
        assert(gradient_expressions[cc].size() >= dimRange);
        for (size_t rr = 0; rr < dimRange; ++rr) {
          const auto& gradient_expression = gradient_expressions[cc][rr];
          assert(gradient_expression.size() >= dimDomain);
          gradients_[cc].emplace_back(new MathExpressionGradientType(variable, gradient_expression));
        }
      }
    }
  } // ... build_gradients(...)

  template <size_t cols = rangeDimCols, bool anything = false>
  struct get_helper
  {
    static void get_expression(const Common::Configuration& cfg, ExpressionStringVectorType& expression_as_vectors)
    {
      typedef typename XT::Common::FieldMatrix<std::string, dimRange, dimRangeCols> ExpressionMatrixType;
      const ExpressionMatrixType expression_as_matrix =
          cfg.get<ExpressionMatrixType>("expression", dimRange, dimRangeCols);
      // convert FieldMatrix to ExpressionStringVectorType
      for (size_t rr = 0; rr < dimRange; ++rr) {
        std::vector<std::string> expression_row;
        for (size_t cc = 0; cc < dimRangeCols; ++cc)
          expression_row.emplace_back(expression_as_matrix[rr][cc]);
        expression_as_vectors.emplace_back(expression_row);
      }
    }
  }; // struct get_helper

  template <bool anything>
  struct get_helper<1, anything>
  {
    static void get_expression(const Common::Configuration& cfg, ExpressionStringVectorType& expression_as_vectors)
    {
      typedef typename Dune::FieldVector<std::string, dimRange> ExpressionVectorType;
      const ExpressionVectorType expression_as_vector = cfg.get<ExpressionVectorType>("expression");
      // convert Vector to ExpressionStringVectorType
      for (size_t rr = 0; rr < dimRange; ++rr) {
        std::vector<std::string> expression_row(1, expression_as_vector[rr]);
        expression_as_vectors.emplace_back(expression_row);
      }
    }
  }; // struct get_helper<1, ...>

  static void get_gradient(const Common::Configuration& cfg,
                           GradientStringVectorType& gradient_as_vectors,
                           const std::string first_gradient_key)
  {
    // create vector of gradient keys
    std::vector<std::string> gradient_keys(1, first_gradient_key);
    for (size_t cc = 1; cc < dimRangeCols; ++cc)
      gradient_keys.emplace_back("gradient." + Common::to_string(cc));
    // get gradient as FieldMatrix for every key
    for (std::string key : gradient_keys) {
      ExpressionStringVectorType gradient_as_vectors_component;
      typedef typename XT::Common::FieldMatrix<std::string, dimRange, dimDomain> JacobianMatrixType;
      const JacobianMatrixType gradient_as_matrix = cfg.get<JacobianMatrixType>(key);
      // convert FieldMatrix to ExpressionStringVectorType
      for (size_t rr = 0; rr < dimRange; ++rr) {
        std::vector<std::string> gradient_expression;
        for (size_t ii = 0; ii < dimDomain; ++ii)
          gradient_expression.emplace_back(gradient_as_matrix[rr][ii]);
        gradient_as_vectors_component.emplace_back(gradient_expression);
      }
      gradient_as_vectors.emplace_back(gradient_as_vectors_component);
    }
  } // ... get_gradient(...)

  std::shared_ptr<const MathExpressionFunctionType> function_;
  size_t order_;
  std::string name_;
  mutable typename Common::PerThreadValue<FieldVector<RangeFieldType, dimRange * dimRangeCols>> tmp_vector_;
  mutable typename Common::PerThreadValue<FieldVector<RangeFieldType, dimRangeCols>> tmp_row_;
  std::vector<std::vector<std::shared_ptr<const MathExpressionGradientType>>> gradients_;
}; // class ExpressionFunction


} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_EXPRESSION_DEFAULT_HH
