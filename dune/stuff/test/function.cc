#include "test_common.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function.hh>

using namespace Dune::Stuff;

typedef testing::Types<Function::Expression<double, 1, double, 1>, Function::Checkerboard<double, 1, double, 1>>
    NonparametricFunctions;

typedef testing::Types<Function::SeparableCheckerboard<double, 1, double, 1>,
                       Function::SeparableDefault<double, 1, double, 1>> SeparableFunctions;


template <class T>
struct NonparametricTest : public ::testing::Test
{
  typedef T FunctionType;
  typedef typename FunctionType::DomainFieldType DomainFieldType;
  static const int dimDomain = FunctionType::dimDomain;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeFieldType RangeFieldType;
  static const int dimRange = FunctionType::dimRange;
  typedef typename FunctionType::RangeType RangeType;
  typedef Function::Interface<DomainFieldType, dimDomain, RangeFieldType, dimRange> InterfaceType;

  void check() const
  {
    const std::unique_ptr<InterfaceType> function(
        Function::create<DomainFieldType, dimDomain, RangeFieldType, dimRange>(
            FunctionType::id(), FunctionType::createSampleDescription()));
    DomainType x(1);
    RangeType ret;
    if (function->parametric())
      DUNE_THROW(Dune::InvalidStateException, "ERROR: nonparametric function returned parametric() == true");
    const std::string DUNE_UNUSED(name) = function->name();
    const int DUNE_UNUSED(order) = function->order();
    function->evaluate(x, ret);
  }
};


template <class T>
struct SeparableTest : public ::testing::Test
{
  typedef T FunctionType;
  typedef typename FunctionType::DomainFieldType DomainFieldType;
  static const int dimDomain = FunctionType::dimDomain;
  typedef typename FunctionType::DomainType DomainType;
  typedef typename FunctionType::RangeFieldType RangeFieldType;
  static const int dimRange = FunctionType::dimRange;
  typedef typename FunctionType::RangeType RangeType;

  typedef Function::Interface<DomainFieldType, dimDomain, RangeFieldType, dimRange> InterfaceType;
  typedef Function::SeparableDefault<DomainFieldType, dimDomain, RangeFieldType, dimRange> DefaultType;

  typedef Common::Parameter::FieldType ParamFieldType;
  static const int maxParamDim = Common::Parameter::maxDim;
  typedef Common::Parameter::Type ParamType;

  void check() const
  {
    const std::unique_ptr<InterfaceType> function(
        FunctionType::createFromDescription(FunctionType::createSampleDescription()));
    if (!function->parametric())
      DUNE_THROW(Dune::InvalidStateException, "ERROR: separable function returned parametric() == false!");
    if (!function->separable())
      DUNE_THROW(Dune::InvalidStateException, "ERROR: separable function returned separable() == false!");
    const std::string name = function->name();
    const int order        = function->order();
    const size_t paramSize = function->paramSize();
    ParamType mu(paramSize);
    const std::vector<ParamType> paramRange = function->paramRange();
    if (paramRange.size() != 2)
      DUNE_THROW(Dune::InvalidStateException, "ERROR: paramRange() has wrong size!");
    if (paramRange[0].size() != paramSize)
      DUNE_THROW(Dune::InvalidStateException, "ERROR: paramRange()[0] has wrong size!");
    if (paramRange[1].size() != paramSize)
      DUNE_THROW(Dune::InvalidStateException, "ERROR: paramRange()[1] has wrong size!");
    for (size_t pp = 0; pp < paramSize; ++pp)
      mu[pp] = paramRange[0][pp] + 0.5 * (paramRange[1][pp] - paramRange[0][pp]);
    DomainType x(1);
    RangeType ret;
    function->evaluate(x, mu, ret);
    const std::vector<std::string>& paramExplanation = function->paramExplanation();
    if (paramExplanation.size() != paramSize)
      DUNE_THROW(Dune::InvalidStateException, "ERROR: paramExplanation() has wrong size!");
    DefaultType DUNE_UNUSED(separableDefault)(
        paramSize, paramRange, function->components(), function->coefficients(), paramExplanation, order, name);
  }
};


TYPED_TEST_CASE(NonparametricTest, NonparametricFunctions);
TYPED_TEST(NonparametricTest, Nonparametric)
{
  this->check();
}


TYPED_TEST_CASE(SeparableTest, SeparableFunctions);
TYPED_TEST(SeparableTest, Separable)
{
  this->check();
}


int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  DSC::Logger().create(DSC::LOG_CONSOLE | DSC::LOG_ERROR);
  return RUN_ALL_TESTS();
}
