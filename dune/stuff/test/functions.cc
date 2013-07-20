#include "test_common.hh"

#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/fixed.hh>


using namespace Dune::Stuff;


typedef testing::Types<FunctionExpression<double, 1, double, 1>, FunctionCheckerboard<double, 1, double, 1>,
                       FunctionCheckerboard<double, 2, double, 1>, FunctionCheckerboard<double, 3, double, 1>,
                       FunctionConstant<double, 1, double, 1>, FunctionConstant<double, 2, double, 1>,
                       FunctionConstant<double, 3, double, 1>
                       //                      , FunctionSpe10Model1< double, 2, double, 1 > // <- this makes only
                       //                      sense, if the data file is present
                       > Functions;

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
    const std::shared_ptr<const FunctionType> function(FunctionType::create(FunctionType::defaultSettings()));
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    function->evaluate(x, ret);
  }
}; // struct FunctionTest


TYPED_TEST_CASE(FunctionTest, Functions);
TYPED_TEST(FunctionTest, Function)
{
  this->check();
}


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
    const std::shared_ptr<const FunctionType> function(FunctionType::create(FunctionType::defaultSettings()));
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    function->evaluate(x, t, ret);
    const FixedTimeFunctionType fixedFunction(function, 0);
    const std::string DUNE_UNUSED(f_name) = fixedFunction.name();
    const int DUNE_UNUSED(f_order) = fixedFunction.order();
    fixedFunction.evaluate(x, ret);
  }
}; // struct TimedependentFunctionTest


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
