#include "test_common.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/function/expression.hh>
#include <dune/stuff/function/checkerboard.hh>
#include <dune/stuff/function/spe10.hh>


using namespace Dune::Stuff;


typedef testing::Types<FunctionExpression<double, 1, double, 1>, FunctionCheckerboard<double, 1, double, 1>,
                       FunctionCheckerboard<double, 2, double, 1>, FunctionCheckerboard<double, 3, double, 1>
                       //                      , FunctionSpe10Model1< double, 2, double, 1 > // <- this makes only
                       //                      sense, if the data file is present
                       > Functions;

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
};


TYPED_TEST_CASE(FunctionTest, Functions);
TYPED_TEST(FunctionTest, Function)
{
  this->check();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
