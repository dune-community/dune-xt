#include "test_common.hh"

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/stuff/function.hh>


typedef testing::Types<Dune::Stuff::Function::Expression<double, 1, double, 1>,
                       Dune::Stuff::Function::Checkerboard<double, 1, double, 1>> NonparametricFunctions;

typedef testing::Types<Dune::Stuff::Function::SeparableCheckerboard<double, 1, double, 1>,
                       Dune::Stuff::Function::SeparableDefault<double, 1, double, 1>> SeparableFunctions;


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
  typedef Dune::Stuff::Function::Interface<DomainFieldType, dimDomain, RangeFieldType, dimRange> InterfaceType;

  void check() const
  {
    const Dune::shared_ptr<InterfaceType> function =
        Dune::Stuff::Function::create<DomainFieldType, dimDomain, RangeFieldType, dimRange>(
            FunctionType::id(), FunctionType::createSampleDescription());
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

  typedef Dune::Stuff::Function::Interface<DomainFieldType, dimDomain, RangeFieldType, dimRange> InterfaceType;
  typedef Dune::Stuff::Function::SeparableDefault<DomainFieldType, dimDomain, RangeFieldType, dimRange> DefaultType;

  typedef Dune::Stuff::Common::Parameter::FieldType ParamFieldType;
  static const int maxParamDim = Dune::Stuff::Common::Parameter::maxDim;
  typedef Dune::Stuff::Common::Parameter::Type ParamType;

  void check() const
  {
    const Dune::shared_ptr<InterfaceType> function =
        Dune::Stuff::Function::create<DomainFieldType, dimDomain, RangeFieldType, dimRange>(
            FunctionType::id(), FunctionType::createSampleDescription());
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
  return RUN_ALL_TESTS();
}
