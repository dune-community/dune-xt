// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/checkerboard.hh>


// we need this nasty code generation because the testing::Types< ... > only accepts 50 arguments
// and all combinations of functions and entities and dimensions and fieldtypes would be way too much
#define TEST_STRUCT_GENERATOR(ftype, etype)                                                                            \
  template <class LocalizableFunctionType>                                                                             \
  struct ftype##etype##Test : public ::testing::Test                                                                   \
  {                                                                                                                    \
    typedef typename LocalizableFunctionType::EntityType EntityType;                                                   \
    typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;                                     \
    typedef typename LocalizableFunctionType::DomainFieldType DomainFieldType;                                         \
    static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;                                          \
    typedef typename LocalizableFunctionType::DomainType DomainType;                                                   \
    typedef typename LocalizableFunctionType::RangeFieldType RangeFieldType;                                           \
    static const unsigned int dimRange     = LocalizableFunctionType::dimRange;                                        \
    static const unsigned int dimRangeCols = LocalizableFunctionType::dimRangeCols;                                    \
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

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/sgrid.hh>
#include <dune/stuff/common/reenable_warnings.hh>

typedef Dune::SGrid<1, 1>::Codim<0>::Entity DuneSGrid1dEntityType;
typedef Dune::SGrid<2, 2>::Codim<0>::Entity DuneSGrid2dEntityType;
typedef Dune::SGrid<3, 3>::Codim<0>::Entity DuneSGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 3, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 3, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 3, 3>>
    CheckerboardFunctionSGridEntityTypes;

TEST_STRUCT_GENERATOR(CheckerboardFunction, SGridEntity)
TYPED_TEST_CASE(CheckerboardFunctionSGridEntityTest, CheckerboardFunctionSGridEntityTypes);
TYPED_TEST(CheckerboardFunctionSGridEntityTest, provides_required_methods)
{
  this->check();
}

#include <dune/grid/yaspgrid.hh>

typedef Dune::YaspGrid<1>::Codim<0>::Entity DuneYaspGrid1dEntityType;
typedef Dune::YaspGrid<2>::Codim<0>::Entity DuneYaspGrid2dEntityType;
typedef Dune::YaspGrid<3>::Codim<0>::Entity DuneYaspGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 3, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 3, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 3, 3>>
    CheckerboardFunctionYaspGridEntityTypes;

TEST_STRUCT_GENERATOR(CheckerboardFunction, YaspGridEntity)
TYPED_TEST_CASE(CheckerboardFunctionYaspGridEntityTest, CheckerboardFunctionYaspGridEntityTypes);
TYPED_TEST(CheckerboardFunctionYaspGridEntityTest, provides_required_methods)
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

typedef testing::Types<Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 3>

                       ,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 3>

                       ,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Functions::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 3>>
    CheckerboardFunctionAluGridEntityTypes;

TEST_STRUCT_GENERATOR(CheckerboardFunction, AluGridEntity)
TYPED_TEST_CASE(CheckerboardFunctionAluGridEntityTest, CheckerboardFunctionAluGridEntityTypes);
TYPED_TEST(CheckerboardFunctionAluGridEntityTest, provides_required_methods)
{
  this->check();
}

#endif // HAVE_ALUGRID_SERIAL || HAVE_ALUGRID_PARALLEL
#endif // HAVE_DUNE_GRID
