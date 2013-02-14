#include "test_common.hh"

#if HAVE_DUNE_GRID

#include <iostream>
#include <fstream>
#include <utility>

#include <boost/filesystem.hpp>

#include <dune/common/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/timer.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/grid/provider.hh>

using namespace Dune;
using namespace Dune::Stuff;

const std::string id = "grid_provider";

ParameterTree sampleParamTree()
{
  ParameterTree paramTree;
  paramTree[id + ".grid"]                           = "stuff.grid.provider.cube";
  paramTree["stuff.grid.provider.cube.lowerLeft"]   = "[0.0; 0.0; 0.0]";
  paramTree["stuff.grid.provider.cube.upperRight"]  = "[1.0; 1.0; 1.0]";
  paramTree["stuff.grid.provider.cube.numElements"] = "[1; 1; 1]";
  return paramTree;
}

static const int dim = 2;

typedef testing::Types<Dune::YaspGrid<dim>
#if HAVE_ALUGRID
                       ,
                       Dune::ALUCubeGrid<dim, dim>, Dune::ALUConformGrid<dim, dim>, Dune::ALUSimplexGrid<dim, dim>
#endif
#if HAVE_ALBERTA
                       ,
                       Dune::AlbertaGrid<dim>
#endif
#if HAVE_UG
                       ,
                       Dune::UGGrid<dim>
#endif
                       ,
                       Dune::SGrid<dim, dim>> GridTypes;

template <class GridType>
struct CubeTest : public testing::Test
{
  typedef Dune::FieldVector<typename GridType::ctype, dim> CoordinateType;

  CubeTest()
  {
  }

  void test_cube(const ParameterTree& paramTree)
  {
    const Dune::Stuff::Grid::Provider::Interface<>* gridProvider =
        Dune::Stuff::Grid::Provider::create(paramTree.get(id + ".grid", "stuff.grid.provider.cube"), paramTree);
    EXPECT_GT(gridProvider->grid()->size(0), 0);
    EXPECT_GT(gridProvider->grid()->size(1), 0);
  }
};

TYPED_TEST_CASE(CubeTest, GridTypes);
TYPED_TEST(CubeTest, All)
{
  ParameterTree paramTree = sampleParamTree();
  this->test_cube(paramTree);
}

#endif // #if HAVE_DUNE_GRID

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
