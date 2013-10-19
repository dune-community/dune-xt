#include "test_common.hh"

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
          LocalizableFunctionType::create(LocalizableFunctionType::defaultSettings()));                                \
    }                                                                                                                  \
  };
// TEST_STRUCT_GENERATOR


#ifdef HAVE_DUNE_GRID

#include <dune/grid/sgrid.hh>

typedef typename Dune::SGrid<1, 1>::template Codim<0>::Entity DuneSGrid1dEntityType;
typedef typename Dune::SGrid<2, 2>::template Codim<0>::Entity DuneSGrid2dEntityType;
typedef typename Dune::SGrid<3, 3>::template Codim<0>::Entity DuneSGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid1dEntityType, double, 1, double, 3, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid2dEntityType, double, 2, double, 3, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneSGrid3dEntityType, double, 3, double, 3, 3>>
    CheckerboardFunctionSGridEntityTypes;

TEST_STRUCT_GENERATOR(CheckerboardFunction, SGridEntity)
TYPED_TEST_CASE(CheckerboardFunctionSGridEntityTest, CheckerboardFunctionSGridEntityTypes);
TYPED_TEST(CheckerboardFunctionSGridEntityTest, provides_required_methods)
{
  this->check();
}

#include <dune/grid/yaspgrid.hh>

typedef typename Dune::YaspGrid<1>::template Codim<0>::Entity DuneYaspGrid1dEntityType;
typedef typename Dune::YaspGrid<2>::template Codim<0>::Entity DuneYaspGrid2dEntityType;
typedef typename Dune::YaspGrid<3>::template Codim<0>::Entity DuneYaspGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid1dEntityType, double, 1, double, 3, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid2dEntityType, double, 2, double, 3, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneYaspGrid3dEntityType, double, 3, double, 3, 3>>
    CheckerboardFunctionYaspGridEntityTypes;

TEST_STRUCT_GENERATOR(CheckerboardFunction, YaspGridEntity)
TYPED_TEST_CASE(CheckerboardFunctionYaspGridEntityTest, CheckerboardFunctionYaspGridEntityTypes);
TYPED_TEST(CheckerboardFunctionYaspGridEntityTest, provides_required_methods)
{
  this->check();
}

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H

#undef HAVE_GRIDTYPE
#undef WORLDDIM
#undef GRIDDIM
#undef GRIDTYPE
#undef ENABLE_ALUGRID
#undef HAVE_ALUGRID
#define ENABLE_ALUGRID 1
#define HAVE_ALUGRID 1

#include <dune/grid/alugrid.hh>

typedef typename Dune::ALUSimplexGrid<2, 2>::template Codim<0>::Entity DuneAluSimplexGrid2dEntityType;
typedef typename Dune::ALUSimplexGrid<3, 3>::template Codim<0>::Entity DuneAluSimplexGrid3dEntityType;
typedef typename Dune::ALUCubeGrid<3, 3>::template Codim<0>::Entity DuneAluCubeGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 3>

                       ,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 3>

                       ,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 1, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 2, 3>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 1>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 2>,
                       Dune::Stuff::Function::Checkerboard<DuneAluCubeGrid3dEntityType, double, 3, double, 3, 3>>
    CheckerboardFunctionAluGridEntityTypes;

TEST_STRUCT_GENERATOR(CheckerboardFunction, AluGridEntity)
TYPED_TEST_CASE(CheckerboardFunctionAluGridEntityTest, CheckerboardFunctionAluGridEntityTypes);
TYPED_TEST(CheckerboardFunctionAluGridEntityTest, provides_required_methods)
{
  this->check();
}

#endif // HAVE_ALUGRID_SERIAL || HAVE_ALUGRID_PARALLEL
#endif // HAVE_DUNE_GRID


int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
