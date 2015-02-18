// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "main.hxx"

#if HAVE_DUNE_GRID
#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/grid/provider/interface.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune;
using namespace Stuff;


#define SGRIDS SGrid<1, 1>, SGrid<2, 2>, SGrid<3, 3>, SGrid<4, 4>, SGrid<1, 2>, SGrid<2, 3>, SGrid<3, 4>, SGrid<4, 5>

#define YASPGRIDS YaspGrid<1>, YaspGrid<2>, YaspGrid<3>, YaspGrid<4>

#if HAVE_ALUGRID
#define ALUGRIDS                                                                                                       \
  ALUGrid<2, 2, cube, nonconforming> /*, ALUGrid< 2, 3, cube, nonconforming >     <- provider does not work */         \
      , ALUGrid<3, 3, cube, nonconforming>,                                                                            \
      ALUGrid<2, 2, simplex, nonconforming> /*, ALUGrid< 2, 3, simplex, nonconforming >  <- same here */               \
      , ALUGrid<3, 3, simplex, nonconforming>,                                                                         \
      ALUGrid<2, 2, simplex, conforming> /*, ALUGrid< 2, 3, simplex, conforming >     <- same here */                  \
      , ALUGrid<3, 3, simplex, conforming>
#endif // HAVE_ALUGRID


template <class GridProviderType>
class ConstGridProviderBase : public testing::Test
{
protected:
  typedef typename GridProviderType::GridType GridType;

private:
  typedef Stuff::Grid::ConstProviderInterface<GridType> ConstInterfaceType;
  typedef Stuff::Grid::ProviderInterface<GridType> InterfaceType;

  static_assert(std::is_base_of<InterfaceType, GridProviderType>::value,
                "GridProviderType has to be derived from Stuff::Grid::ConstProviderInterface!");

protected:
  static std::unique_ptr<GridProviderType> create()
  {
    return GridProviderType::create();
  }

public:
  static void is_default_creatable()
  {
    std::unique_ptr<GridProviderType> DUNE_UNUSED(grd_prvdr) = create();
  }

  static void static_interface()
  {
    std::string static_id_s                         = GridProviderType::static_id();
    auto grid_provider_d                            = create();
    std::string static_id_d                         = grid_provider_d->static_id();
    std::unique_ptr<InterfaceType>& grid_provider_i = grid_provider_d;
    std::string static_id_i                         = grid_provider_i->static_id();
  } // ... static_interface(...)

private:
  template <class GPT>
  static void check_const_interface(const GPT& grid_provider)
  {
    using Stuff::Grid::ChooseLayer;
    using Stuff::Grid::ChoosePartView;
    const auto& grid_i = grid_provider.grid();
    auto DUNE_UNUSED(leaf_grid_view_1) = grid_provider.leaf_view();
    auto DUNE_UNUSED(leaf_grid_view_2) = grid_provider.template leaf<ChoosePartView::view>();
    // auto leaf_grid_part = grid_provider.template leaf< ChoosePartView::part >();                       // <does not
    // work
    for (int level = 0; level <= grid_i.maxLevel(); ++level) {
      auto DUNE_UNUSED(level_grid_view_1) = grid_provider.level_view(level);
      auto DUNE_UNUSED(level_grid_view_2) = grid_provider.template level<ChoosePartView::view>(level);
      // auto DUNE_UNUSED(level_grid_part) = grid_provider.template level< ChoosePartView::part >(level); // <does not
      // work
    }
  } // ... check_const_interface()

public:
  static void const_interface()
  {
    const auto grid_provider = create();
#if defined(__INTEL_COMPILER)
    /**
     * ICC is actually correct and standard compliant in rejecting the #else branch here:
     * 7.3.3 The using declaration [namespace.udecl]
     * 14/ If a function declaration in namespace scope or block scope has the same name and
     * the same parameter types as a function introduced by a using-declaration, and the
     * declarations do not declare the same function, the program is ill-formed.
     *
     * This is just a temporary workaround to have tests build on all compilers again
     **/
    check_const_interface<ConstInterfaceType>(*grid_provider);
#else
    check_const_interface<GridProviderType>(*grid_provider);
    check_const_interface<InterfaceType>(*grid_provider);
#endif
  } // ... const_interface()

private:
  template <class GPT>
  static void check_visualize(const GPT& grid_provider)
  {
    auto type_str = Stuff::Common::Typename<GPT>::value();
    grid_provider.visualize();
    grid_provider.visualize(type_str + "_a");
    grid_provider.visualize(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config());
    grid_provider.visualize(Stuff::Grid::BoundaryInfoConfigs::AllDirichlet::default_config(), type_str + "_b");
  } // ... check_visualize(...)

public:
  static void visualize()
  {
    if (GridType::dimension <= 3) {
      const auto grid_provider = create();
      check_visualize<GridProviderType>(*grid_provider);
      check_visualize<InterfaceType>(*grid_provider);
    }
  } // ... visualize(...)

}; // class ConstGridProviderBase


template <class GridProviderType>
class GridProviderBase : public ConstGridProviderBase<GridProviderType>
{
  typedef ConstGridProviderBase<GridProviderType> BaseType;
  typedef typename BaseType::GridType GridType;
  typedef Stuff::Grid::ProviderInterface<GridType> InterfaceType;

private:
  template <class GPT>
  static void check_non_const_interface(GPT& grid_provider)
  {
    using Stuff::Grid::ChooseLayer;
    using Stuff::Grid::ChoosePartView;
    auto& grid_i = grid_provider.grid();
    auto DUNE_UNUSED(leaf_grid_part_1) = grid_provider.template leaf<ChoosePartView::part>();
    auto DUNE_UNUSED(leaf_grid_part_2) = grid_provider.leaf_part();
    for (int level = 0; level <= grid_i.maxLevel(); ++level) {
      auto DUNE_UNUSED(level_grid_part_1) = grid_provider.template level<ChoosePartView::part>(level);
      auto DUNE_UNUSED(level_grid_part_2) = grid_provider.level_part(level);
    }
  } // ... check_non_const_interface()

public:
  static void non_const_interface()
  {
    auto grid_provider = BaseType::create();
    check_non_const_interface<GridProviderType>(*grid_provider);
    check_non_const_interface<InterfaceType>(*grid_provider);
  } // ... non_const_interface()

}; // class GridProviderBase


template <class GridType>
struct CubeGridProvider : public GridProviderBase<Stuff::Grid::Providers::Cube<GridType>>
{
};


typedef testing::Types<SGRIDS, YASPGRIDS
#if HAVE_ALUGRID
                       ,
                       ALUGRIDS
#endif
                       > GridTypes;

TYPED_TEST_CASE(CubeGridProvider, GridTypes);
TYPED_TEST(CubeGridProvider, is_default_creatable)
{
  this->is_default_creatable();
}
TYPED_TEST(CubeGridProvider, fulfills_const_interface)
{
  this->const_interface();
}
TYPED_TEST(CubeGridProvider, is_visualizable)
{
  this->visualize();
}

#if HAVE_DUNE_FEM

TYPED_TEST(CubeGridProvider, fulfills_non_const_interface)
{
  this->non_const_interface();
}

#else // HAVE_DUNE_FEM

TEST(DISABLED_CubeGridProvider, fulfills_non_const_interface)
{
}

#endif


#else // HAVE_DUNE_GRID

TEST(DISABLED_CubeGridProvider, is_default_creatable)
{
}
TEST(DISABLED_CubeGridProvider, fulfills_const_interface)
{
}
TEST(DISABLED_CubeGridProvider, is_visualizable)
{
}

#endif // HAVE_DUNE_GRID
