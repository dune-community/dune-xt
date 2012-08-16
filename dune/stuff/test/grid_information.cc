#include "test_common.hh"

#include <dune/stuff/grid/information.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune::Stuff;
using namespace std;

template <int i>
struct Int
{
  static const int value = i;
};

typedef testing::Types<Int<1>, Int<2>, Int<3>> GridDims;

template <class T>
struct GridInfoTest : public ::testing::Test
{
  static const int dim = T::value;
  typedef Dune::YaspGrid<dim, dim> GridType;
  Dune::shared_ptr<GridType> gridPtr;
  GridInfoTest()
    : gridPtr(Grid::Provider::UnitCube(/*level=*/5).gridPtr())
  {
  }

  void checkDimension()
  {
  }
};

TYPED_TEST_CASE(GridInfoTest, GridDims);
TYPED_TEST(GridInfoTest, Dimension)
{
}

// TEST_F(GridInfoTest, Dimension) {
//  Grid::Provider::UnitCube cube;
//}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  Dune::MPIHelper::instance(argc, argv);
  return RUN_ALL_TESTS();
}
