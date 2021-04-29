// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012, 2014, 2016 - 2017, 2020)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Tobias Leibner  (2014, 2016, 2018)

/**
 * This file is intended as a starting point for quick testing.
 */

#include <dune/xt/test/main.hxx> // <- has to come first, include config.h!

#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>
#include <dune/xt/grid/dd/glued.hh>

#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/functors/generic.hh>

#include <dune/xt/grid/view/coupling.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>


/**
 * Element vector for testing purposes
 */
template <class E>
class ElementVector
{
public:
  ElementVector(std::vector<E>& vector)
    : vector_{vector}
  {}
  typename std::vector<E>::iterator begin() const
  {
    return vector_.begin();
  };
  typename std::vector<E>::iterator end() const
  {
    return vector_.end();
  };

private:
  std::vector<E>& vector_;
};

/**
 * Element vector and their intersections for testing purposes
 */
template <class E, class I, class MacroGridViewType>
class ElementVectorAndIntersections
{
private:
  size_t subdomain(const E& macro_entity) const
  {
    //      assert_macro_grid_state();
    assert(macro_grid_view_.indexSet().contains(macro_entity));
    return macro_grid_view_.indexSet().index(macro_entity);
  }

public:
  ElementVectorAndIntersections(MacroGridViewType& mgv,
                                std::vector<E>& elements,
                                std::vector<std::set<I, Dune::XT::Grid::internal::CompareType<I>>>& intersections_)
    : macro_grid_view_(mgv)
    , elements_{elements}
    , intersections_{intersections_}
  {}

  typename std::vector<E>::iterator begin() const
  {
    return elements_.begin();
  };
  typename std::vector<E>::iterator end() const
  {
    return elements_.end();
  };
  auto ibegin(E& element) const
  {
    return intersections_[subdomain(element)].begin();
  };
  auto iend(E& element) const
  {
    return intersections_[subdomain(element)].end();
  };

private:
  MacroGridViewType& macro_grid_view_;
  std::vector<E>& elements_;
  std::vector<std::set<I, Dune::XT::Grid::internal::CompareType<I>>>& intersections_;
};

/**
 * Element vector and their coupling intersections for testing purposes
 */
template <class E, class I, class MacroGridViewType>
class ElementVectorAndCouplingIntersections
{
private:
  size_t subdomain(const E& macro_entity) const
  {
    //      assert_macro_grid_state();
    assert(macro_grid_view_.indexSet().contains(macro_entity));
    return macro_grid_view_.indexSet().index(macro_entity);
  }

public:
  ElementVectorAndCouplingIntersections(
      MacroGridViewType& mgv,
      std::vector<E>& elements,
      std::vector<std::set<I, Dune::XT::Grid::internal::CompareType<I>>>& intersections_)
    : macro_grid_view_(mgv)
    , elements_{elements}
    , intersections_{intersections_}
  {}

  typename std::vector<E>::iterator begin() const
  {
    return elements_.begin();
  };
  typename std::vector<E>::iterator end() const
  {
    return elements_.end();
  };
  typename std::set<I, Dune::XT::Grid::internal::CompareType<I>>::iterator ibegin(E& element) const
  {
    return intersections_[subdomain(element)].begin();
  };
  typename std::set<I, Dune::XT::Grid::internal::CompareType<I>>::iterator iend(E& element) const
  {
    return intersections_[subdomain(element)].end();
  };

private:
  MacroGridViewType& macro_grid_view_;
  std::vector<E>& elements_;
  std::vector<std::set<I, Dune::XT::Grid::internal::CompareType<I>>>& intersections_;
};

/**
 * Inside element vector and their coupling intersections for testing purposes
 */
template <class E, class I, class MacroGridViewType>
class InsideElementsAndCouplingIntersections
{
private:
  size_t subdomain(const E& macro_entity) const
  {
    //      assert_macro_grid_state();
    assert(macro_grid_view_.indexSet().contains(macro_entity));
    return macro_grid_view_.indexSet().index(macro_entity);
  }

public:
  InsideElementsAndCouplingIntersections(
      MacroGridViewType& mgv,
      std::vector<E>& elements,
      std::vector<std::set<I, Dune::XT::Grid::internal::CompareType<I>>>& coupling_intersections,
      std::vector<std::vector<E>>& inside_elements)
    : macro_grid_view_(mgv)
    , elements_{elements}
    , coupling_intersections_{coupling_intersections}
    , inside_elements_(inside_elements)
  {}

  typename std::vector<E>::iterator begin() const
  {
    return elements_.begin();
  };
  typename std::vector<E>::iterator end() const
  {
    return elements_.end();
  };
  typename std::set<I, Dune::XT::Grid::internal::CompareType<I>>::iterator ibegin(E& element) const
  {
    return coupling_intersections_[subdomain(element)].begin();
  };
  typename std::set<I, Dune::XT::Grid::internal::CompareType<I>>::iterator iend(E& element) const
  {
    return coupling_intersections_[subdomain(element)].end();
  };
  typename std::vector<E>::iterator ebegin(E& element) const
  {
    return inside_elements_[subdomain(element)].begin();
  };
  typename std::vector<E>::iterator eend(E& element) const
  {
    return inside_elements_[subdomain(element)].end();
  };

private:
  MacroGridViewType& macro_grid_view_;
  std::vector<E>& elements_;
  std::vector<std::set<I, Dune::XT::Grid::internal::CompareType<I>>>& coupling_intersections_;
  std::vector<std::vector<E>>& inside_elements_;
};


GTEST_TEST(empty, main)
{
  using GridType = YASP_2D_EQUIDISTANT_OFFSET;
  using ElementType = typename GridType::template Codim<0>::Entity;
  static constexpr size_t d = GridType::dimension;
  using D = typename GridType::ctype;
  using GridProviderType = Dune::XT::Grid::GridProvider<GridType>;
  using GridViewType = typename GridProviderType::LeafGridViewType;

  auto grid = Dune::XT::Grid::make_cube_grid<GridType>(0., 1., 2);
  auto dd_grid = std::make_unique<Dune::XT::Grid::DD::Glued<GridType, GridType, Dune::XT::Grid::Layers::leaf>>(grid, 2);
  auto local_grid = dd_grid->local_grid(0);

  std::cout << "dimension of the grid:   " << grid.dimDomain << std::endl;
  std::cout << "dimension of local grid: " << local_grid.dimDomain << std::endl;
  std::cout << "number of subdomains:    " << dd_grid->num_subdomains() << std::endl;

  GridViewType mgv = dd_grid->macro_grid_view();

  std::cout << "\n"
            << "walking the marco grid ..." << std::endl;
  std::vector<ElementType> elements_;
  for (auto&& macro_element : Dune::elements(mgv)) {
    elements_.push_back(macro_element);
    const auto subdomain_id = dd_grid->subdomain(macro_element);
    std::cout << "id of the element: " << subdomain_id << std::endl;
  }

  // using iterators from a class !!
  ElementVector<ElementType> vector(elements_);

  using GridGlueType = Dune::XT::Grid::DD::Glued<GridType, GridType, Dune::XT::Grid::Layers::leaf>;
  using GlueType = typename GridGlueType::GlueType;
  using IntersectionType = typename GridViewType::Intersection;

  std::vector<std::set<IntersectionType, Dune::XT::Grid::internal::CompareType<IntersectionType>>> intersections_;

  std::cout << "\n"
            << "walking the marco grid and walk the interesections..." << std::endl;
  for (auto el = vector.begin(); el != vector.end(); ++el) {
    std::cout << "center of element: " << (*el).geometry().center() << "\n";
    std::set<IntersectionType, Dune::XT::Grid::internal::CompareType<IntersectionType>> intersection_set;
    for (auto& macro_intersection : Dune::intersections(mgv, (*el))) {
      auto outside_element = macro_intersection.outside();
      std::cout << "center of outside element " << outside_element.geometry().center() << std::endl;
      intersection_set.insert(macro_intersection);
      std::cout << "intersections centers: {" << std::endl;
      for (auto it : intersection_set)
        std::cout << it.geometry().center() << std::endl;
      std::cout << "}" << std::endl;
    }
    intersections_.push_back(intersection_set);
  }

  ElementVectorAndIntersections<ElementType, IntersectionType, GridViewType> element_and_intersections(
      mgv, elements_, intersections_);

  std::cout << "\n"
            << "walking the marco grid and walk the interesections..." << std::endl;
  for (auto el = element_and_intersections.begin(); el != element_and_intersections.end(); ++el) {
    std::cout << "center of element: " << (*el).geometry().center() << "\n";
    for (auto it = element_and_intersections.ibegin(*el); it != element_and_intersections.iend(*el); ++it) {
      std::cout << "center of intersection: " << (*it).geometry().center() << "\n";
    }
  }
  std::cout << "\n" << std::endl;

  using CouplingIntersectionType = typename GlueType::Intersection;
  using CorrectedCouplingIntersectionType =
      Dune::XT::Grid::internal::CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, IntersectionType>;

  std::vector<std::set<CorrectedCouplingIntersectionType,
                       Dune::XT::Grid::internal::CompareType<CorrectedCouplingIntersectionType>>>
      coupling_intersections_;
  std::set<CorrectedCouplingIntersectionType, Dune::XT::Grid::internal::CompareType<CorrectedCouplingIntersectionType>>
      coupling_intersection_set;

  for (auto el = element_and_intersections.begin(); el != element_and_intersections.end(); ++el) {
    for (auto it = element_and_intersections.ibegin(*el); it != element_and_intersections.iend(*el); ++it) {
      auto inside_element = it->inside();
      auto outside_element = it->outside();
      std::cout << "center of inside element: " << inside_element.geometry().center() << "\n";
      std::cout << "center of outside element: " << outside_element.geometry().center() << "\n";

      if (it->boundary()) {
        std::cout << "skip this element, it is not inside the domain" << std::endl;
        continue;
      }
      // NOTE: shouldn't coupling() throw an assertion if the outside_element is not inside the domain as it is for
      // boundary intersection??? The above if statement is a workaround for that
      auto& coupling = dd_grid->coupling(inside_element, -1, outside_element, -1, true);
      for (auto coupling_intersection_it = coupling.template ibegin<0>();
           coupling_intersection_it != coupling.template iend<0>();
           ++coupling_intersection_it) {
        std::cout << "current center is " << coupling_intersection_it->geometry().center() << std::endl;
        CorrectedCouplingIntersectionType coupling_intersection(*coupling_intersection_it, *it);
        coupling_intersection_set.insert(coupling_intersection);
      }
    }
    coupling_intersections_.push_back(coupling_intersection_set);
    coupling_intersection_set.clear();
  }

  std::cout << (*next(intersections_.at(3).begin(), 1)).geometry().center() << std::endl;
  std::cout << (*next(coupling_intersections_.at(0).begin(), 0)).inside().geometry().center() << std::endl;

  std::vector<std::vector<ElementType>> inside_elements;
  ElementVectorAndCouplingIntersections<ElementType, CorrectedCouplingIntersectionType, GridViewType>
      element_and_coupling_intersections(mgv, elements_, coupling_intersections_);

  for (auto el = element_and_coupling_intersections.begin(); el != element_and_coupling_intersections.end(); ++el) {
    std::cout << "center of element: " << (*el).geometry().center() << "\n";
    std::vector<ElementType> inside_elements_on_element;
    for (auto it = element_and_coupling_intersections.ibegin(*el); it != element_and_coupling_intersections.iend(*el);
         ++it) {
      std::cout << "center of coupling intersection: " << (*it).geometry().center() << "\n";
      auto inside_element = it->inside();
      inside_elements_on_element.push_back(inside_element);
    }
    inside_elements.push_back(inside_elements_on_element);
  }

  InsideElementsAndCouplingIntersections<ElementType, CorrectedCouplingIntersectionType, GridViewType>
      inside_element_and_intersections(mgv, elements_, coupling_intersections_, inside_elements);
  for (auto el = inside_element_and_intersections.begin(); el != inside_element_and_intersections.end(); ++el) {
    std::cout << "center of element: " << (*el).geometry().center() << "\n";
    for (auto it = inside_element_and_intersections.ebegin(*el); it != inside_element_and_intersections.eend(*el);
         ++it) {
      std::cout << "center of inside element: " << (*it).geometry().center() << "\n";
    }
  }

  using CouplingGridViewType = Dune::XT::Grid::CouplingGridView<GridGlueType>;

  // use the gridview from view/coupling.hh
  for (auto&& macro_element : Dune::elements(mgv)) {
    for (auto& macro_intersection : Dune::intersections(mgv, macro_element)) {
      printf("_______\n");
      if (macro_intersection.boundary()) {
        std::cout << "skip this intersection, it is not inside the domain" << std::endl;
        continue;
      }
      auto inside_element = macro_intersection.inside();
      auto outside_element = macro_intersection.outside();
      auto cgv = Dune::XT::Grid::make_coupling_grid_view<GridGlueType, ElementType, IntersectionType>(
          inside_element, outside_element, *dd_grid, macro_intersection);
      for (auto&& inside_element : Dune::elements(cgv)) {
        std::cout << "I arrived at an element with center: " << inside_element.geometry().center() << std::endl;
        for (auto&& intersection : Dune::intersections(cgv, inside_element)) {
          std::cout << "I arrived at an intersection with center: " << intersection.geometry().center() << std::endl;
        }
      }
      //            auto functor = Dune::XT::Grid::GenericElementFunctor<CouplingGridViewType>([] {},
      //                [](const auto& element) {
      //                std::cout << "WALKER: Element : " << element.geometry().center() << std::endl; },
      //                [] {});

      auto functor = Dune::XT::Grid::GenericElementAndIntersectionFunctor<CouplingGridViewType>(
          [] {},
          [](const auto& element) { std::cout << "WALKER: Element : " << element.geometry().center() << std::endl; },
          [](const auto& intersection, const ElementType&, const ElementType&) {
            std::cout << "WALKER:      Intersection : " << intersection.geometry().center() << std::endl;
          },
          [] {});

      // these lines are not required !!! but they test the CouplingGridProvider
      Dune::XT::Grid::CouplingGridProvider<CouplingGridViewType> coupling_provider(cgv);
      const auto cgv_c = coupling_provider.coupling_view(); // with & would be reference, without is copy
      CouplingGridViewType cgv_c2(cgv_c); // <-- calls copy constructor
      auto walker = Dune::XT::Grid::Walker<CouplingGridViewType>(cgv_c); // or cgv itself
      walker.append(functor);
      walker.walk();
    }
  }
}
