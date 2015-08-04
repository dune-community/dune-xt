// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/indicator.hh>


// we need this nasty code generation because the testing::Types< ... > only accepts 50 arguments
// and all combinations of functions and entities and dimensions and fieldtypes would be way too much
#define TEST_STRUCT_GENERATOR(ftype, etype)                                                                            \
  template <class LocalizableFunctionType>                                                                             \
  struct ftype##etype##Test : public ::testing::Test                                                                   \
  {                                                                                                                    \
    typedef typename LocalizableFunctionType::EntityType EntityType;                                                   \
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;                                     \
    typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;                                         \
    static const size_t dimDomain = LocalizableFunctionType::dimDomain;                                                \
    typedef typename LocalizableFunctionType::DomainType DomainType;                                                   \
    typedef typename LocalizableFunctionType::RangeFieldType RangeFieldType;                                           \
    static const size_t dimRange     = LocalizableFunctionType::dimRange;                                              \
    static const size_t dimRangeCols = LocalizableFunctionType::dimRangeCols;                                          \
    typedef typename LocalizableFunctionType::RangeType RangeType;                                                     \
    typedef typename LocalizableFunctionType::JacobianRangeType JacobianRangeType;                                     \
                                                                                                                       \
    void check() const                                                                                                 \
    {                                                                                                                  \
      const std::unique_ptr<const LocalizableFunctionType> function(                                                   \
          LocalizableFunctionType::create(LocalizableFunctionType::default_config()));                                 \
    }                                                                                                                  \
  };
// TEST_STRUCT_GENERATOR


#if HAVE_DUNE_GRID

#include <dune/grid/yaspgrid.hh>

typedef Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>::Codim<0>::Entity DuneYaspGrid1dEntityType;
typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>::Codim<0>::Entity DuneYaspGrid2dEntityType;
typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>::Codim<0>::Entity DuneYaspGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Functions::Indicator<DuneYaspGrid1dEntityType, double, 1, double, 1, 1>,
                       Dune::Stuff::Functions::Indicator<DuneYaspGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Functions::Indicator<DuneYaspGrid3dEntityType, double, 3, double, 1, 1>>
    IndicatorFunctionYaspGridEntityTypes;

TEST_STRUCT_GENERATOR(IndicatorFunction, YaspGridEntity)
TYPED_TEST_CASE(IndicatorFunctionYaspGridEntityTest, IndicatorFunctionYaspGridEntityTypes);
TYPED_TEST(IndicatorFunctionYaspGridEntityTest, provides_required_methods)
{
  this->check();
}

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>::Codim<0>::Entity DuneAluSimplexGrid2dEntityType;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>::Codim<0>::Entity DuneAluSimplexGrid3dEntityType;
typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>::Codim<0>::Entity DuneAluCubeGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Functions::Indicator<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Functions::Indicator<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Functions::Indicator<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 1>>
    IndicatorFunctionAluGridEntityTypes;

TEST_STRUCT_GENERATOR(IndicatorFunction, AluGridEntity)
TYPED_TEST_CASE(IndicatorFunctionAluGridEntityTest, IndicatorFunctionAluGridEntityTypes);
TYPED_TEST(IndicatorFunctionAluGridEntityTest, provides_required_methods)
{
  this->check();
}

#endif // HAVE_ALUGRID_SERIAL || HAVE_ALUGRID_PARALLEL
#endif // HAVE_DUNE_GRID
