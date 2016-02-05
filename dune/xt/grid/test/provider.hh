// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2016)
//   Rene Milk       (2012 - 2016)
//   Tobias Leibner  (2014 - 2016)

#ifndef DUNE_XT_TEST_GRID_PROVIDER_HH
#define DUNE_XT_TEST_GRID_PROVIDER_HH

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/provider/interface.hh>

#include <dune/xt/common/test/gtest/gtest.h>

template <class GridProviderType>
class ConstGridProviderBase : public testing::Test
{
protected:
  typedef typename GridProviderType::GridType GridType;

private:
  typedef Dune::XT::Grid::ConstProviderInterface<GridType> ConstInterfaceType;
  typedef Dune::XT::Grid::ProviderInterface<GridType> InterfaceType;

  static_assert(std::is_base_of<InterfaceType, GridProviderType>::value,
                "GridProviderType has to be derived from Dune::XT::Grid::ConstProviderInterface!");

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
    using Dune::XT::Grid::ChooseLayer;
    using Dune::XT::Grid::ChoosePartView;
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
    auto type_str = Dune::XT::Common::Typename<GPT>::value();
    grid_provider.visualize();
    grid_provider.visualize(type_str + "_a");
    grid_provider.visualize(Dune::XT::Grid::alldirichlet_boundaryinfo_default_config());
    grid_provider.visualize(Dune::XT::Grid::alldirichlet_boundaryinfo_default_config(), type_str + "_b");
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
  typedef Dune::XT::Grid::ProviderInterface<GridType> InterfaceType;

private:
  template <class GPT>
  static void check_non_const_interface(GPT& grid_provider)
  {
    using Dune::XT::Grid::ChooseLayer;
    using Dune::XT::Grid::ChoosePartView;
    auto& grid_i = grid_provider.grid();
#if HAVE_DUNE_FEM
    auto DUNE_UNUSED(leaf_grid_part_1) = grid_provider.template leaf<ChoosePartView::part>();
    auto DUNE_UNUSED(leaf_grid_part_2) = grid_provider.leaf_part();
    for (int level = 0; level <= grid_i.maxLevel(); ++level) {
      auto DUNE_UNUSED(level_grid_part_1) = grid_provider.template level<ChoosePartView::part>(level);
      auto DUNE_UNUSED(level_grid_part_2) = grid_provider.level_part(level);
    }
#else // HAVE_DUNE_FEM
    auto DUNE_UNUSED(leaf_grid_view_1) = grid_provider.template leaf<ChoosePartView::view>();
    auto DUNE_UNUSED(leaf_grid_view_2) = grid_provider.leaf_view();
    for (int level = 0; level <= grid_i.maxLevel(); ++level) {
      auto DUNE_UNUSED(level_grid_view_1) = grid_provider.template level<ChoosePartView::view>(level);
      auto DUNE_UNUSED(level_grid_view_2) = grid_provider.level_view(level);
    }
#endif
  } // ... check_non_const_interface()

public:
  static void non_const_interface()
  {
    auto grid_provider = BaseType::create();
    check_non_const_interface<GridProviderType>(*grid_provider);
    check_non_const_interface<InterfaceType>(*grid_provider);
  } // ... non_const_interface()

}; // class GridProviderBase

#endif // DUNE_XT_TEST_GRID_PROVIDER_HH
