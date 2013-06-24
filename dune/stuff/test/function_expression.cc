#include "test_common.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/common/tuple.hh>
#include <dune/stuff/function/expression/base.hh>


struct RunExpressionBaseTest
{
  template <class DimDomain, class DimRange>
  static void run()
  {
    using namespace Dune;
    using namespace Dune::Stuff;
    typedef double DomainFieldType;
    static const size_t dimDomain = DimDomain::value;
    typedef double RangeFieldType;
    static const size_t dimRange = DimRange::value;

    typedef FunctionExpressionBase<DomainFieldType, dimDomain, RangeFieldType, dimRange> FunctionType;
    const std::string variable                 = "x";
    const std::vector<std::string> expressions = {"0", "1", "2"};
    const FunctionType function(variable, expressions);
    // test fieldvector/fieldvector
    FieldVector<DomainFieldType, dimDomain> argFieldVector;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      argFieldVector[ii] = RangeFieldType(ii);
    FieldVector<DomainFieldType, dimRange> retFieldVector;
    function.evaluate(argFieldVector, retFieldVector);
    for (size_t ii = 0; ii < dimRange; ++ii)
      if (retFieldVector[ii] < ii || retFieldVector[ii] > ii)
        DUNE_THROW(RangeError,
                   "\nERROR (dimDomain = " << dimDomain << ", dimRange = " << dimRange
                                           << "): wrong result,  retFieldVector["
                                           << ii
                                           << "] is "
                                           << retFieldVector[ii]
                                           << ", should be "
                                           << ii
                                           << "!");
    // test dynvector/dynvector
    DynamicVector<DomainFieldType> argDynVector(dimDomain);
    DynamicVector<DomainFieldType> retDynVector;
    for (size_t dd = 1; dd <= dimDomain; ++dd) {
      argDynVector = DynamicVector<DomainFieldType>(dd);
      retDynVector = DynamicVector<DomainFieldType>();
      for (size_t ii = 0; ii < dd; ++ii)
        argDynVector[ii] = RangeFieldType(ii);
      function.evaluate(argDynVector, retDynVector);
      if (retDynVector.size() != dimRange)
        DUNE_THROW(RangeError,
                   "\nERROR: retDynVector has wrong size, is " << retDynVector.size() << ", should be " << dimRange
                                                               << "!");
      for (size_t ii = 0; ii < dimRange; ++ii)
        if (retDynVector[ii] < ii || retDynVector[ii] > ii)
          DUNE_THROW(RangeError,
                     "\nERROR (dimDomain = " << dimDomain << ", dimRange = " << dimRange
                                             << "): wrong result,  retDynVector["
                                             << ii
                                             << "] is "
                                             << retFieldVector[ii]
                                             << ", should be "
                                             << ii
                                             << "!");
    }
    // test fieldvector/dynvector
    retDynVector = DynamicVector<DomainFieldType>();
    function.evaluate(argFieldVector, retDynVector);
    if (retDynVector.size() != dimRange)
      DUNE_THROW(RangeError,
                 "\nERROR: retDynVector has wrong size, is " << retDynVector.size() << ", should be " << dimRange
                                                             << "!");
    for (size_t ii = 0; ii < dimRange; ++ii)
      if (retDynVector[ii] < ii || retDynVector[ii] > ii)
        DUNE_THROW(RangeError,
                   "\nERROR (dimDomain = " << dimDomain << ", dimRange = " << dimRange
                                           << "): wrong result,  retDynVector["
                                           << ii
                                           << "] is "
                                           << retDynVector[ii]
                                           << ", should be "
                                           << ii
                                           << "!");
    // test dynvector/fieldvector
    for (size_t dd = 1; dd <= dimDomain; ++dd) {
      argDynVector   = DynamicVector<DomainFieldType>(dd);
      retFieldVector = FieldVector<DomainFieldType, dimRange>();
      for (size_t ii = 0; ii < dd; ++ii)
        argDynVector[ii] = RangeFieldType(ii);
      function.evaluate(argDynVector, retFieldVector);
      for (size_t ii = 0; ii < dimRange; ++ii)
        if (retFieldVector[ii] < ii || retFieldVector[ii] > ii)
          DUNE_THROW(RangeError,
                     "\nERROR (dimDomain = " << dimDomain << ", dimRange = " << dimRange
                                             << "): wrong result,  retFieldVector["
                                             << ii
                                             << "] is "
                                             << retFieldVector[ii]
                                             << ", should be "
                                             << ii
                                             << "!");
    }
  }
}; // struct RunExpressionBaseTest


template <class TestFunctor>
struct FunctionTest : public ::testing::Test
{
  typedef boost::mpl::vector<Int<1>, Int<2>, Int<3>> DomainDims;
  typedef DomainDims RangeDims;
  typedef
      typename DSC::TupleProduct::Combine<DomainDims, RangeDims, TestFunctor>::template Generate<> base_generator_type;
  void run()
  {
    base_generator_type::Run();
  }
}; // struct FunctionTest


typedef ::testing::Types<RunExpressionBaseTest> FunctionTestTypes;
TYPED_TEST_CASE(FunctionTest, FunctionTestTypes);
TYPED_TEST(FunctionTest, All)
{
  this->run();
}


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
