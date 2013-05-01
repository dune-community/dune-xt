#include "test_common.hh"

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/function/spe10.hh>
#include <dune/stuff/function/constant.hh>
#include <dune/stuff/function/affineparametric/checkerboard.hh>
#include <dune/stuff/function/fixed.hh>


using namespace Dune::Stuff;


typedef testing::Types<FunctionExpression<double, 1, double, 1>, FunctionCheckerboard<double, 1, double, 1>,
                       FunctionCheckerboard<double, 2, double, 1>, FunctionCheckerboard<double, 3, double, 1>,
                       FunctionConstant<double, 1, double, 1>, FunctionConstant<double, 2, double, 1>,
                       FunctionConstant<double, 3, double, 1>
                       //                      , FunctionSpe10Model1< double, 2, double, 1 > // <- this makes only
                       //                      sense, if the data file is present
                       > Functions;

typedef testing::Types<AffineParametricFunctionCheckerboard<double, 1, double>,
                       AffineParametricFunctionCheckerboard<double, 2, double>,
                       AffineParametricFunctionCheckerboard<double, 3, double>> AffineParametricFunctions;

typedef testing::Types<FunctionConstant<double, 1, double, 1>, FunctionConstant<double, 2, double, 1>,
                       FunctionConstant<double, 3, double, 1>> TimedependentFunctions;

template <class FunctionType>
struct FunctionTest : public ::testing::Test
{
  typedef typename FunctionType::DomainFieldType DomainFieldType;
  static const int dimDomain = FunctionType::dimDomain;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeFieldType RangeFieldType;
  static const int dimRangeRows = FunctionType::dimRangeRows;
  static const int dimRangeCols = FunctionType::dimRangeCols;
  typedef typename FunctionType::RangeType RangeType;

  void check() const
  {
    DomainType x(1);
    RangeType ret(0);
    const std::shared_ptr<const FunctionType> function(FunctionType::create(FunctionType::createSampleDescription()));
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    function->evaluate(x, ret);
  }
}; // struct FunctionTest


template <class FunctionType>
struct ParametricFunctionTest : public ::testing::Test
{
  typedef typename FunctionType::DomainFieldType DomainFieldType;
  static const int dimDomain = FunctionType::dimDomain;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeFieldType RangeFieldType;
  static const int dimRangeRows = FunctionType::dimRangeRows;
  static const int dimRangeCols = FunctionType::dimRangeCols;
  typedef typename FunctionType::RangeType RangeType;
  typedef typename FunctionType::ParamFieldType ParamFieldType;
  static const int maxParamDim = FunctionType::maxParamDim;
  typedef typename FunctionType::ParamType ParamType;

  typedef FunctionFixedParameter<DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
      FixedParameterFunctionType;

  void check() const
  {
    DomainType x(1);
    RangeType ret(0);
    const std::shared_ptr<const FunctionType> function(FunctionType::create(FunctionType::createSampleDescription()));
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    const bool parametric = function->parametric();
    assert(parametric);
    const size_t paramSize = function->paramSize();
    assert(paramSize > 0);
    const std::vector<ParamType>& paramRange = function->paramRange();
    assert(paramRange.size() == 2);
    assert(paramRange[0].size() == paramSize);
    assert(paramRange[1].size() == paramSize);
    const std::vector<std::string>& paramExplanation = function->paramExplanation();
    assert(paramExplanation.size() == paramSize);
    function->evaluate(x, paramRange[0], ret);
    const FixedParameterFunctionType fixedFunction(function, paramRange[0]);
    const std::string DUNE_UNUSED(f_name) = fixedFunction.name();
    const int DUNE_UNUSED(f_order) = fixedFunction.order();
    const bool f_parametric = fixedFunction.parametric();
    assert(!f_parametric);
    fixedFunction.evaluate(x, ret);
  }
}; // struct ParametricFunctionTest


template <class FunctionType>
struct AffineParametricFunctionTest : public ::testing::Test
{
  typedef typename FunctionType::DomainFieldType DomainFieldType;
  static const int dimDomain = FunctionType::dimDomain;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeFieldType RangeFieldType;
  static const int dimRangeRows = FunctionType::dimRangeRows;
  static const int dimRangeCols = FunctionType::dimRangeCols;
  typedef typename FunctionType::RangeType RangeType;
  typedef typename FunctionType::ParamFieldType ParamFieldType;
  static const int maxParamDim = FunctionType::maxParamDim;
  typedef typename FunctionType::ParamType ParamType;
  typedef typename FunctionType::ComponentType ComponentType;
  typedef typename FunctionType::CoefficientType CoefficientType;

  void check() const
  {
    DomainType x(1);
    RangeType ret(0);
    const std::shared_ptr<const FunctionType> function(FunctionType::create(FunctionType::createSampleDescription()));
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    const bool parametric = function->parametric();
    assert(parametric);
    const size_t paramSize = function->paramSize();
    assert(paramSize > 0);
    const std::vector<ParamType>& paramRange = function->paramRange();
    assert(paramRange.size() == 2);
    assert(paramRange[0].size() == paramSize);
    assert(paramRange[1].size() == paramSize);
    const std::vector<std::string>& paramExplanation = function->paramExplanation();
    assert(paramExplanation.size() == paramSize);
    function->evaluate(x, paramRange[0], ret);
    const bool affineParametric = function->affineparametric();
    assert(affineParametric);
    const std::vector<std::shared_ptr<const ComponentType>>& components     = function->components();
    const std::vector<std::shared_ptr<const CoefficientType>>& coefficients = function->coefficients();
    assert(components.size() == coefficients.size());
    const bool hasAffinePart = function->hasAffinePart();
    if (hasAffinePart)
      const std::shared_ptr<ComponentType>& DUNE_UNUSED(affinePart) = function->affinePart();
  }
}; // struct AffineParametricFunctionTest


template <class FunctionType>
struct TimedependentFunctionTest : public ::testing::Test
{
  typedef typename FunctionType::DomainFieldType DomainFieldType;
  static const int dimDomain = FunctionType::dimDomain;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeFieldType RangeFieldType;
  static const int dimRangeRows = FunctionType::dimRangeRows;
  static const int dimRangeCols = FunctionType::dimRangeCols;
  typedef typename FunctionType::RangeType RangeType;

  typedef FunctionFixedTime<DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
      FixedTimeFunctionType;

  void check() const
  {
    DomainType x(1);
    RangeType ret(0);
    double t(0);
    const std::shared_ptr<const FunctionType> function(FunctionType::create(FunctionType::createSampleDescription()));
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    function->evaluate(x, t, ret);
    const FixedTimeFunctionType fixedFunction(function, 0);
    const std::string DUNE_UNUSED(f_name) = fixedFunction.name();
    const int DUNE_UNUSED(f_order) = fixedFunction.order();
    fixedFunction.evaluate(x, ret);
  }
}; // struct TimedependentFunctionTest


TYPED_TEST_CASE(FunctionTest, Functions);
TYPED_TEST(FunctionTest, Function)
{
  this->check();
}


TYPED_TEST_CASE(ParametricFunctionTest, AffineParametricFunctions);
TYPED_TEST(ParametricFunctionTest, ParametricFunction)
{
  this->check();
}


TYPED_TEST_CASE(AffineParametricFunctionTest, AffineParametricFunctions);
TYPED_TEST(AffineParametricFunctionTest, AffineParametricFunction)
{
  this->check();
}


TYPED_TEST_CASE(TimedependentFunctionTest, TimedependentFunctions);
TYPED_TEST(TimedependentFunctionTest, TimedependentFunction)
{
  this->check();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
