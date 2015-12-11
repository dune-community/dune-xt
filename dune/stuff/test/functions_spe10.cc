// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Rene Milk (2015)

#include "main.hxx"

#include <memory>

#include <dune/common/exceptions.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/spe10.hh>
#include <dune/stuff/functions/spe10model2.hh>

// we need this nasty code generation because the testing::Types< ... > only accepts 50 arguments
// and all combinations of functions and entities and dimensions and fieldtypes would be way too much
#define TEST_STRUCT_GENERATOR(ftype, etype)                                                                            \
  template <class Spe10Model2FunctionType>                                                                             \
  struct ftype##etype##Test : public ::testing::Test                                                                   \
  {                                                                                                                    \
    typedef typename Spe10Model2FunctionType::EntityType EntityType;                                                   \
    typedef typename Spe10Model2FunctionType::LocalfunctionType LocalfunctionType;                                     \
    typedef typename Spe10Model2FunctionType::DomainFieldType DomainFieldType;                                         \
    static const size_t dimDomain = Spe10Model2FunctionType::dimDomain;                                                \
    typedef typename Spe10Model2FunctionType::DomainType DomainType;                                                   \
    typedef typename Spe10Model2FunctionType::RangeFieldType RangeFieldType;                                           \
    static const size_t dimRange     = Spe10Model2FunctionType::dimRange;                                              \
    static const size_t dimRangeCols = Spe10Model2FunctionType::dimRangeCols;                                          \
    typedef typename Spe10Model2FunctionType::RangeType RangeType;                                                     \
    typedef typename Spe10Model2FunctionType::JacobianRangeType JacobianRangeType;                                     \
                                                                                                                       \
    void check() const                                                                                                 \
    {                                                                                                                  \
      const Spe10Model2FunctionType zero;                                                                              \
    }                                                                                                                  \
  };
// TEST_STRUCT_GENERATOR

#if HAVE_DUNE_GRID

//# include <dune/grid/yaspgrid.hh>

// typedef Dune::SGrid< 1, 1 >::Codim< 0 >::Entity DuneSGrid1dEntityType;
// typedef Dune::SGrid< 2, 2 >::Codim< 0 >::Entity DuneSGrid2dEntityType;
// typedef Dune::SGrid< 3, 3 >::Codim< 0 >::Entity DuneSGrid3dEntityType;

// typedef testing::Types< Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 1, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 1, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 1, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 2, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 2, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 2, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 3, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 3, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid1dEntityType, double, 1, double, 3, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 1, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 1, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 1, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 2, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 2, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 2, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 3, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 3, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid2dEntityType, double, 2, double, 3, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 1, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 1, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 1, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 2, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 2, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 2, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 3, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 3, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneSGrid3dEntityType, double, 3, double, 3, 3 >
//                        > Spe10Model2FunctionSGridEntityTypes;

// TEST_STRUCT_GENERATOR(Spe10Model2FunctionTest, SGridEntity)
// TYPED_TEST_CASE(Spe10Model2FunctionTestSGridEntityTest, Spe10Model2FunctionSGridEntityTypes);
// TYPED_TEST(Spe10Model2FunctionTestSGridEntityTest, provides_required_methods) {
//  this->check();
//}

#include <dune/grid/yaspgrid.hh>

typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>>::Codim<0>::Entity DuneYaspGrid3dEntityType;

typedef testing::Types<Dune::Stuff::Functions::Spe10::Model2<DuneYaspGrid3dEntityType, double, 3, double, 3, 3>>
    Spe10Model2FunctionYaspGridEntityTypes;

TEST_STRUCT_GENERATOR(Spe10Model2FunctionTest, YaspGridEntity)
TYPED_TEST_CASE(Spe10Model2FunctionTestYaspGridEntityTest, Spe10Model2FunctionYaspGridEntityTypes);
TYPED_TEST(Spe10Model2FunctionTestYaspGridEntityTest, provides_required_methods)
{
  this->check();
}

//# if HAVE_ALUGRID
//#   include <dune/stuff/common/disable_warnings.hh>
//#     include <dune/grid/alugrid.hh>
//#   include <dune/stuff/common/reenable_warnings.hh>

// typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::nonconforming >::Codim< 0 >::Entity DuneAluSimplexGrid2dEntityType;
// typedef Dune::ALUGrid< 3, 3, Dune::simplex, Dune::nonconforming >::Codim< 0 >::Entity DuneAluSimplexGrid3dEntityType;
// typedef Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming >::Codim< 0 >::Entity     DuneAluCubeGrid3dEntityType;

// typedef testing::Types< Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 1,
// 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 2
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 1, 3
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 1
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 2
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 2, 3
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 1
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 2
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid2dEntityType, double, 2, double, 3, 3
//                      >

//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 1
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 2
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 1, 3
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 1
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 2
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 2, 3
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 1
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 2
//                      >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluSimplexGrid3dEntityType, double, 3, double, 3, 3
//                      >

//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 1, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 1, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 1, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 2, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 2, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 2, 3 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 3, 1 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 3, 2 >
//                      , Dune::Stuff::Functions::Spe10::Model2< DuneAluCubeGrid3dEntityType, double, 3, double, 3, 3 >
//                        > Spe10Model2FunctionAluGridEntityTypes;

// TEST_STRUCT_GENERATOR(Spe10Model2FunctionTest, AluGridEntity)
// TYPED_TEST_CASE(Spe10Model2FunctionTestAluGridEntityTest, Spe10Model2FunctionAluGridEntityTypes);
// TYPED_TEST(Spe10Model2FunctionTestAluGridEntityTest, provides_required_methods) {
//  this->check();
//}

//# endif // HAVE_ALUGRID
#endif // HAVE_DUNE_GRID
