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

  void test_cube()
  {
    const std::vector<std::string> types = Dune::Stuff::Grid::Provider::types();
    const ParameterTree description = Dune::Stuff::Grid::Provider::createSampleDescription<GridType>(types[0]);
    const Dune::Stuff::Grid::Provider::Interface<GridType>* gridProvider =
        Dune::Stuff::Grid::Provider::create<GridType>(types[0], description);
    EXPECT_GT(gridProvider->grid()->size(0), 0);
    EXPECT_GT(gridProvider->grid()->size(1), 0);
  }
};

TYPED_TEST_CASE(CubeTest, GridTypes);
TYPED_TEST(CubeTest, All)
{
  this->test_cube();
}

#endif // #if HAVE_DUNE_GRID

int main(int argc, char** argv)
{
  test_init(argc, argv);
  return RUN_ALL_TESTS();
}
