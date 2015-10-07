// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/global.hh>


// we need this nasty code generation because the testing::Types< ... > only accepts 50 arguments
// and all combinations of functions and entities and dimensions and fieldtypes would be way too much
#define TEST_STRUCT_GENERATOR(ftype, etype)                                                                            \
  template <class GlobalFunctionType>                                                                                  \
  struct ftype##etype##Test : public ::testing::Test                                                                   \
  {                                                                                                                    \
    typedef typename GlobalFunctionType::EntityType EntityType;                                                        \
    typedef typename GlobalFunctionType::LocalfunctionType LocalfunctionType;                                          \
    typedef typename GlobalFunctionType::DomainFieldType DomainFieldType;                                              \
    static const size_t dimDomain = GlobalFunctionType::dimDomain;                                                     \
    typedef typename GlobalFunctionType::DomainType DomainType;                                                        \
    typedef typename GlobalFunctionType::RangeFieldType RangeFieldType;                                                \
    static const size_t dimRange     = GlobalFunctionType::dimRange;                                                   \
    static const size_t dimRangeCols = GlobalFunctionType::dimRangeCols;                                               \
    typedef typename GlobalFunctionType::RangeType RangeType;                                                          \
    typedef typename GlobalFunctionType::JacobianRangeType JacobianRangeType;                                          \
                                                                                                                       \
    void check() const                                                                                                 \
    {                                                                                                                  \
      const GlobalFunctionType zero([](DomainType /*x*/) { return RangeType(0); }, 0);                                 \
      const GlobalFunctionType one([](DomainType /*x*/) { return RangeType(1); }, 0);                                  \
      const auto diff = zero - one;                                                                                    \
      const auto xx = DomainType(666);                                                                                 \
      EXPECT_EQ(zero.evaluate(xx), RangeType(0));                                                                      \
      EXPECT_EQ(one.evaluate(xx), RangeType(1));                                                                       \
    }                                                                                                                  \
  };
// TEST_STRUCT_GENERATOR


#if HAVE_DUNE_GRID

#include <dune/grid/yaspgrid.hh>

typedef Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>::Codim<0>::Entity DuneYaspGrid1dEntityType;
typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>::Codim<0>::Entity DuneYaspGrid2dEntityType;
typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>::Codim<0>::Entity DuneYaspGrid3dEntityType;

typedef testing::Types<Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 1, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 1, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 1, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 2, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 2, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 2, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 3, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 3, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid1dEntityType, double, 1, double, 3, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid2dEntityType, double, 2, double, 3, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneYaspGrid3dEntityType, double, 3, double, 3, 3>>
    GlobalLambdaFunctionYaspGridEntityTypes;

TEST_STRUCT_GENERATOR(GlobalLambdaFunctionTest, YaspGridEntity)
TYPED_TEST_CASE(GlobalLambdaFunctionTestYaspGridEntityTest, GlobalLambdaFunctionYaspGridEntityTypes);
TYPED_TEST(GlobalLambdaFunctionTestYaspGridEntityTest, provides_required_methods)
{
  this->check();
}

#if HAVE_ALUGRID
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/alugrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming>::Codim<0>::Entity DuneAluSimplexGrid2dEntityType;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>::Codim<0>::Entity DuneAluSimplexGrid3dEntityType;
typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>::Codim<0>::Entity DuneAluCubeGrid3dEntityType;

typedef testing::Types<Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 3>

                       ,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 3>

                       ,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::GlobalLambdaFunction<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 3>>
    GlobalLambdaFunctionAluGridEntityTypes;

TEST_STRUCT_GENERATOR(GlobalLambdaFunctionTest, AluGridEntity)
TYPED_TEST_CASE(GlobalLambdaFunctionTestAluGridEntityTest, GlobalLambdaFunctionAluGridEntityTypes);
TYPED_TEST(GlobalLambdaFunctionTestAluGridEntityTest, provides_required_methods)
{
  this->check();
}

#endif // HAVE_ALUGRID
#endif // HAVE_DUNE_GRID
