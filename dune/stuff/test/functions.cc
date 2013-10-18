#include "test_common.hh"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/grid/fakeentity.hh>
#include <dune/stuff/functions/interfaces.hh>
//#include <dune/stuff/functions/expression.hh>
//#include <dune/stuff/functions/checkerboard.hh>
//#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/constant.hh>
//#include <dune/stuff/functions/fixed.hh>


using namespace Dune::Stuff;

typedef Grid::FakeEntity<1> DuneStuffFake1dEntityType;

typedef testing::Types<Function::Constant<DuneStuffFake1dEntityType, double, 1, double, 1>> LocalizableFunctionTypes;

template <class LocalizableFunctionType>
struct LocalizableFunctionTest : public ::testing::Test
{
  typedef typename LocalizableFunctionType::EntityType EntityType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
  typedef typename LocalizableFunctionType::DomainType DomainType;
  typedef typename LocalizableFunctionType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = LocalizableFunctionType::dimRange;
  static const unsigned int dimRangeCols = LocalizableFunctionType::dimRangeCols;
  typedef typename LocalizableFunctionType::RangeType RangeType;
  typedef typename LocalizableFunctionType::JacobianRangeType JacobianRangeType;

  void check() const
  {
    const std::unique_ptr<const LocalizableFunctionType> function(
        LocalizableFunctionType::create(LocalizableFunctionType::defaultSettings()));
  }
}; // struct LocalizableFunctionTest


TYPED_TEST_CASE(LocalizableFunctionTest, LocalizableFunctionTypes);
TYPED_TEST(LocalizableFunctionTest, provides_required_methods)
{
  this->check();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
