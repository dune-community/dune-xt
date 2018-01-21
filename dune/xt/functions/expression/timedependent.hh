// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#ifndef DUNE_XT_FUNCTIONS_EXPRESSION_TIMEDEPENDENT_HH
#define DUNE_XT_FUNCTIONS_EXPRESSION_TIMEDEPENDENT_HH

#include <dune/common/deprecated.hh>

#include "default.hh"

namespace Dune {
namespace XT {
namespace Functions {


template <class EntityImp,
          class DomainFieldImp,
          size_t domainDim,
          class RangeFieldImp,
          size_t rangeDim,
          size_t rangeDimCols = 1,
          class TimeFieldImp = double>
class DUNE_DEPRECATED_MSG("Use ParametricExpressionFuntion instead (18.05.2017)!") TimeDependentExpressionFunction
    : TimeDependentFunctionInterface<ExpressionFunction<EntityImp,
                                                        DomainFieldImp,
                                                        domainDim,
                                                        RangeFieldImp,
                                                        rangeDim,
                                                        rangeDimCols>,
                                     TimeFieldImp>
{
  typedef TimeDependentExpressionFunction<EntityImp,
                                          DomainFieldImp,
                                          domainDim,
                                          RangeFieldImp,
                                          rangeDim,
                                          rangeDimCols,
                                          TimeFieldImp>
      ThisType;
  typedef TimeDependentFunctionInterface<ExpressionFunction<EntityImp,
                                                            DomainFieldImp,
                                                            domainDim,
                                                            RangeFieldImp,
                                                            rangeDim,
                                                            rangeDimCols>,
                                         TimeFieldImp>
      BaseType;

public:
  using typename BaseType::TimeIndependentFunctionType;
  using typename BaseType::TimeFieldType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  typedef ExpressionFunction<EntityImp, DomainFieldImp, dimDomain, RangeFieldImp, dimRange, dimRangeCols>
      ExpressionFunctionType;
  typedef typename ExpressionFunctionType::DomainType DomainType;
  typedef typename ExpressionFunctionType::ExpressionStringVectorType ExpressionStringVectorType;
  typedef typename ExpressionFunctionType::GradientStringVectorType GradientStringVectorType;

  static std::string static_id()
  {
    return ExpressionFunctionType::static_id() + ".timedependentexpression";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["variable"] = "x";
    config["expression"] = "[t*x[0] sin(t*x[0]) exp(t+x[0])]";
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
      ExpressionFunctionType::template get_helper<>::get_expression(cfg, expression_as_vectors);
    } catch (std::exception) {
      // if dimRangeCols == 1 and we could not get expression as FieldVector, get it as FieldMatrix with one col
      if (dimRangeCols == 1) { // the 2 here is on purpose, anything > 1 will suffice
        ExpressionFunctionType::template get_helper<2>::get_expression(cfg, expression_as_vectors);
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
   * \brief Creates a TimeDependentExpression function where every component is identical
   */
  TimeDependentExpressionFunction(const std::string variable,
                                  const std::string expression,
                                  const size_t ord = 0,
                                  const std::string nm = static_id(),
                                  const std::vector<std::string> gradient = std::vector<std::string>())
    : variable_(variable)
    , order_(ord)
    , name_(nm)
  {
    // create ExpressionStringVectorType with identical expressions
    const std::vector<std::string> expression_row(dimRangeCols, expression);
    const ExpressionStringVectorType expressions(dimRange, expression_row);
    expressions_ = expressions;
    // create associated gradient vector
    assert(gradient.size() == 0 || gradient.size() >= dimDomain);
    if (gradient.size() > 0) {
      GradientStringVectorType gradient_expressions(dimRangeCols);
      const ExpressionStringVectorType gradient_row(dimRange, gradient);
      for (size_t cc = 0; cc < dimRangeCols; ++cc) {
        gradient_expressions[cc] = gradient_row;
      }
      gradient_expressions_ = gradient_expressions;
    }
  }

  /**
   * \brief Creates a TimeDependentExpression function with dimRangeCols = 1
   */
  TimeDependentExpressionFunction(
      const std::string variable,
      const std::vector<std::string> expressions,
      const size_t ord = default_config().template get<size_t>("order"),
      const std::string nm = static_id(),
      const std::vector<std::vector<std::string>> gradient_expressions = std::vector<std::vector<std::string>>())
    : variable_(variable)
    , order_(ord)
    , name_(nm)
    , expressions_(ExpressionStringVectorType(1, expressions))
  {
    static_assert(dimRangeCols == 1, "This constructor does not make sense for dimRangeCols > 1!");
    GradientStringVectorType gradient_expressions_vec;
    if (gradient_expressions.size() > 0) {
      gradient_expressions_vec.emplace_back(gradient_expressions);
    }
    gradient_expressions_ = gradient_expressions_vec;
  }

  /**
   * \brief Creates a TimeDependentExpression function
   */
  TimeDependentExpressionFunction(const std::string variable,
                                  const ExpressionStringVectorType expressions,
                                  const size_t ord = 0,
                                  const std::string nm = static_id(),
                                  const GradientStringVectorType gradient_expressions = GradientStringVectorType())
    : variable_(variable)
    , order_(ord)
    , name_(nm)
    , expressions_(expressions)
    , gradient_expressions_(gradient_expressions)
  {
  }

  std::string name() const
  {
    return name_;
  }

  size_t order() const
  {
    return order_;
  }

  virtual std::unique_ptr<TimeIndependentFunctionType> evaluate_at_time(const TimeFieldType t) const
  {
    ExpressionStringVectorType expressions_at_time_t;
    for (size_t rr = 0; rr < expressions_.size(); ++rr) {
      expressions_at_time_t.emplace_back(expressions_[rr]);
      for (size_t cc = 0; cc < expressions_[rr].size(); ++cc) {
        replaceAll(expressions_at_time_t[rr][cc], "t", Common::to_string(t));
      }
    }
    GradientStringVectorType gradients_at_time_t;
    for (size_t cc = 0; cc < gradient_expressions_.size(); ++cc) {
      gradients_at_time_t.emplace_back(gradient_expressions_[cc]);
      for (size_t rr = 0; rr < gradient_expressions_[cc].size(); ++rr) {
        for (size_t ii = 0; ii < gradient_expressions_[cc][rr].size(); ++ii)
          replaceAll(gradients_at_time_t[cc][rr][ii], "t", Common::to_string(t));
      }
    }
    return Common::make_unique<ExpressionFunctionType>(
        variable_, expressions_at_time_t, order_, name_, gradients_at_time_t);
  }

private:
  // from https://stackoverflow.com/questions/2896600/how-to-replace-all-occurrences-of-a-character-in-string
  void replaceAll(std::string& source, const std::string from, const std::string to) const
  {
    std::string newString;
    newString.reserve(source.length()); // avoids a few memory allocations

    std::string::size_type lastPos = 0;
    std::string::size_type findPos;

    while (std::string::npos != (findPos = source.find(from, lastPos))) {
      newString.append(source, lastPos, findPos - lastPos);
      newString += to;
      lastPos = findPos + from.length();
    }

    // Care for the rest after last occurrence
    newString += source.substr(lastPos);

    source.swap(newString);
  } // ... replaceAll(...)

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
      typedef typename Dune::FieldMatrix<std::string, dimRange, dimDomain> JacobianMatrixType;
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

  const std::string variable_;
  size_t order_;
  std::string name_;
  ExpressionStringVectorType expressions_;
  GradientStringVectorType gradient_expressions_;
}; // class TimeDependentExpression


} // namespace Functions
} // namespace XT
} // namespace Dune


#endif // DUNE_XT_FUNCTIONS_EXPRESSION_TIMEDEPENDENT_HH
