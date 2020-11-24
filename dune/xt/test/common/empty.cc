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
 * Inherits all types and methods from the coupling intersection, but uses the macro intersection to provide a correctly
 * oriented normal.
 *
 * \attention Presumes that the coupling intersection lies exactly within the macro intersection!
 */
template <class CouplingIntersectionType, class MacroIntersectionType>
class CouplingIntersectionWithCorrectNormal : public CouplingIntersectionType
{
  using BaseType = CouplingIntersectionType;

public:
  using typename BaseType::GlobalCoordinate;
  using typename BaseType::LocalCoordinate;

  enum
  {
    dimension = MacroIntersectionType::dimension
  };

  CouplingIntersectionWithCorrectNormal(const CouplingIntersectionType& coupling_intersection,
                                        const MacroIntersectionType& macro_intersection)
    : BaseType(coupling_intersection)
    , macro_intersection_(macro_intersection)
  {}

  GlobalCoordinate outerNormal(const LocalCoordinate& local) const
  {
    global_ = this->geometry().global(local);
    local_ = macro_intersection_.geometry().local(global_);
    return macro_intersection_.outerNormal(local_);
  }

  GlobalCoordinate integrationOuterNormal(const LocalCoordinate& local) const
  {
    auto normal = this->unitOuterNormal(local);
    const auto integration_element = BaseType::integrationOuterNormal(local).two_norm();
    normal *= integration_element;
    return normal;
  }

  GlobalCoordinate unitOuterNormal(const LocalCoordinate& local) const
  {
    global_ = this->geometry().global(local);
    local_ = macro_intersection_.geometry().local(global_);
    return macro_intersection_.unitOuterNormal(local_);
  }

  GlobalCoordinate centerUnitOuterNormal() const
  {
    global_ = this->geometry().center();
    local_ = macro_intersection_.geometry().local(global_);
    return macro_intersection_.unitOuterNormal(local_);
  }

private:
  const MacroIntersectionType& macro_intersection_;
  mutable GlobalCoordinate global_;
  mutable LocalCoordinate local_;
}; // class CouplingIntersectionWithCorrectNormal

//template <class Intersection>
//struct CompareType
//{
//  bool operator()(const Intersection& one, const Intersection& other) const
//  {
//    return (std::addressof(one) < std::addressof(other));
//    // does there exist a better way for this????

////    std::cout << "one center " << one.geometry().center() << std::endl;
////    std::cout << "other center " << other.geometry().center() << std::endl;
//////    std::cout << "seed " << (one.inside().seed() == other.inside().seed()) << std::endl;
////    auto one_center = one.geometry().center();
////    auto other_center = other.geometry().center();
//////    auto decision = (Dune::XT::Common::FloatCmp::lt(one_center[0], other_center[0]) && Dune::XT::Common::FloatCmp::ne(one_center, other_center));
//////    if (!decision)
//////        std::cout << "I do not add this" << std::endl;
//////    auto decision = Dune::XT::Common::FloatCmp::ne(one.geometry().center(), other.geometry().center());
//////    std::cout << "my decision is " << decision << std::endl;
////    return decision;
//  }
//};

template<class E>
class Iterator {
public:
    Iterator(std::vector<E> init)
        : vector_(init) {
           }
    bool operator==(const Iterator<E>&) const {
        std::cout << "hello" << std::endl;
        return true;
    };
    bool operator!=(const Iterator<E>&) const;
    Iterator<E>& operator++();
    Iterator<E> operator++(int);
    E & operator*() const;
    E & operator->() const;
private:
    std::vector<E>& vector_;
};

template<class E>
class ElementVector {
public:
    ElementVector(std::vector<E>& vector)
        : vector_{vector} {}
    typename std::vector<E>::iterator begin() const { return vector_.begin(); };
    typename std::vector<E>::iterator end() const { return vector_.end(); };
private:
    std::vector<E>& vector_;
};

template<class E, class I, class MacroGridViewType>
class ElementVectorAndIntersections {
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
                                  std::vector<std::set<I, CompareType<I>>>& intersections_)
        : macro_grid_view_(mgv),
          elements_{elements},
          intersections_{intersections_} {}

    typename std::vector<E>::iterator begin() const { return elements_.begin(); };
    typename std::vector<E>::iterator end() const { return elements_.end(); };
    typename std::set<I, CompareType<I>>::iterator ibegin(E& element) const { return intersections_[subdomain(element)].begin(); };
    typename std::set<I, CompareType<I>>::iterator iend(E& element) const { return intersections_[subdomain(element)].end(); };
private:
    MacroGridViewType& macro_grid_view_;
    std::vector<E>& elements_;
    std::vector<std::set<I, CompareType<I>>>& intersections_;
};

template<class E, class I, class MacroGridViewType>
class ElementVectorAndCouplingIntersections {
private:
    size_t subdomain(const E& macro_entity) const
    {
//      assert_macro_grid_state();
      assert(macro_grid_view_.indexSet().contains(macro_entity));
      return macro_grid_view_.indexSet().index(macro_entity);
    }
public:
    ElementVectorAndCouplingIntersections(MacroGridViewType& mgv,
                                          std::vector<E>& elements,
                                          std::vector<std::set<I, CompareType<I>>>& intersections_)
        : macro_grid_view_(mgv),
          elements_{elements},
          intersections_{intersections_} {}

    typename std::vector<E>::iterator begin() const { return elements_.begin(); };
    typename std::vector<E>::iterator end() const { return elements_.end(); };
    typename std::set<I, CompareType<I>>::iterator ibegin(E& element) const { return intersections_[subdomain(element)].begin(); };
    typename std::set<I, CompareType<I>>::iterator iend(E& element) const { return intersections_[subdomain(element)].end(); };
private:
    MacroGridViewType& macro_grid_view_;
    std::vector<E>& elements_;
    std::vector<std::set<I, CompareType<I>>>& intersections_;
};

template<class E, class I, class MacroGridViewType>
class InsideElementsAndCouplingIntersections {
private:
    size_t subdomain(const E& macro_entity) const
    {
//      assert_macro_grid_state();
      assert(macro_grid_view_.indexSet().contains(macro_entity));
      return macro_grid_view_.indexSet().index(macro_entity);
    }
public:
    InsideElementsAndCouplingIntersections(MacroGridViewType& mgv,
                                           std::vector<E>& elements,
                                           std::vector<std::set<I, CompareType<I>>>& coupling_intersections,
                                           std::vector<std::vector<E>>& inside_elements)
        : macro_grid_view_(mgv),
          elements_{elements},
          coupling_intersections_{coupling_intersections},
          inside_elements_(inside_elements) {}

    typename std::vector<E>::iterator begin() const { return elements_.begin(); };
    typename std::vector<E>::iterator end() const { return elements_.end(); };
    typename std::set<I, CompareType<I>>::iterator ibegin(E& element) const { return coupling_intersections_[subdomain(element)].begin(); };
    typename std::set<I, CompareType<I>>::iterator iend(E& element) const { return coupling_intersections_[subdomain(element)].end(); };
    typename std::vector<E>::iterator ebegin(E& element) const { return inside_elements_[subdomain(element)].begin(); };
    typename std::vector<E>::iterator eend(E& element) const { return inside_elements_[subdomain(element)].end(); };
private:
    MacroGridViewType& macro_grid_view_;
    std::vector<E>& elements_;
    std::vector<std::set<I, CompareType<I>>>& coupling_intersections_;
    std::vector<std::vector<E>>& inside_elements_;
};

template<class GridGlueType>
class TentativeGridView {
public:
    using MacroGridViewType = typename GridGlueType::MacroGridViewType;
    using MacroGridType = typename GridGlueType::MacroGridType;
    using MacroElementType = typename MacroGridType::template Codim<0>::Entity;
    using MacroIntersectionType = typename MacroGridViewType::Intersection;
    using GlueType = typename GridGlueType::GlueType;
    using CouplingIntersectionType = typename GlueType::Intersection;
    using CorrectedCouplingIntersectionType = CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroIntersectionType>;

    TentativeGridView(GridGlueType& dd_grid)  // <- insert this as a pointer?
        : dd_grid_(dd_grid),
          macro_grid_view_(dd_grid.macro_grid_view())
    {
        for (auto& macro_element : Dune::elements(macro_grid_view_)) {
          macro_elements_.push_back(macro_element);
          std::set<MacroIntersectionType, CompareType<MacroIntersectionType>> intersection_set;
          for (auto& macro_intersection : Dune::intersections(macro_grid_view_, (macro_element))) {
              intersection_set.insert(macro_intersection);
          }
          macro_intersections_.push_back(intersection_set);
        }
    }
    typename std::vector<MacroElementType>::iterator mbegin() { return macro_elements_.begin(); };
    typename std::vector<MacroElementType>::iterator mend() { return macro_elements_.end(); };
    typename std::set<MacroIntersectionType, CompareType<MacroIntersectionType>>::iterator mibegin(MacroElementType& element) {
        return macro_intersections_[subdomain(element)].begin();
    };
    typename std::set<MacroIntersectionType, CompareType<MacroIntersectionType>>::iterator miend(MacroElementType& element) {
        return macro_intersections_[subdomain(element)].end();
    };
    // Note: stoped hear because I realized that this is not neccessary. A coupling gridview only needs to be defined wrt a neighbor.
private:
    size_t subdomain(const MacroElementType& macro_entity) const
    {
//      assert_macro_grid_state();
      assert(macro_grid_view_.indexSet().contains(macro_entity));
      return macro_grid_view_.indexSet().index(macro_entity);
    }

private:
    const GridGlueType& dd_grid_;
    const MacroGridViewType& macro_grid_view_;
    std::vector<MacroElementType> macro_elements_;
    std::vector<std::set<MacroIntersectionType, CompareType<MacroIntersectionType>>> macro_intersections_;
};


template<class GridGlueType>
class TentativeCouplingGridView {
public:
    using MacroGridViewType = typename GridGlueType::MacroGridViewType;
    using MacroGridType = typename GridGlueType::MacroGridType;
    using MacroElementType = typename MacroGridType::template Codim<0>::Entity;
    using MacroIntersectionType = typename MacroGridViewType::Intersection;

    using LocalGridProviderType = typename GridGlueType::LocalGridProviderType;
    using LocalGridViewType = typename GridGlueType::MicroGridViewType;
    using LocalElementType = typename GridGlueType::MicroEntityType;

    using GlueType = typename GridGlueType::GlueType;
    using CouplingIntersectionType = typename GlueType::Intersection;
    // TODO: add the macro intersection to use CorrectedCouplingIntersection
    using CorrectedCouplingIntersectionType = CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroIntersectionType>;

    TentativeCouplingGridView(MacroElementType& ss, MacroElementType& nn, GridGlueType& dd_grid)  // <- insert this as a pointer?
        : inside_element_(ss),
          outside_element_(nn),
          dd_grid_(dd_grid),
          macro_grid_view_(dd_grid.macro_grid_view()),
          local_inside_grid_(dd_grid.local_grid(ss))
    {
        auto& coupling = dd_grid.coupling(ss, -1, nn, -1, true);
        for (auto coupling_intersection_it = coupling.template ibegin<0>();
             coupling_intersection_it != coupling.template iend<0>();
             ++coupling_intersection_it) {
            auto inside = coupling_intersection_it->inside();
            inside_elements_ids_.push_back(local_inside_grid_->leaf_view().indexSet().index(inside));
        }
        // some coupling intersection may have the same inside element, remove duplicates
        std::sort(inside_elements_ids_.begin(), inside_elements_ids_.end());
        auto last = std::unique(inside_elements_ids_.begin(), inside_elements_ids_.end());
        inside_elements_ids_.erase(last, inside_elements_ids_.end());

        for (auto id = inside_elements_ids_.begin(); id != inside_elements_ids_.end(); ++id) {
            for (auto el = local_inside_grid_->leaf_view().template begin<0>();
                 el != local_inside_grid_->leaf_view().template end<0>(); ++el) {
                if (local_inside_grid_->leaf_view().indexSet().index(*el) == *id) {
                    // This is the inside element we are searching for.. add it to the vector
                    inside_elements_.push_back(*el);
                    std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>> coupling_intersection_set;
                    // now iteratate over all intersections to find all coupling intersections
                    for (auto coupling_intersection_it = coupling.template ibegin<0>();
                         coupling_intersection_it != coupling.template iend<0>();
                         ++coupling_intersection_it) {
                        auto inside = coupling_intersection_it->inside();
                        auto inside_id = local_inside_grid_->leaf_view().indexSet().index(inside);
                        if (inside_id == *id) {
//                            CorrectedCouplingIntersectionType coupling_intersection(*coupling_intersection_it, macro_intersection_);
                            coupling_intersection_set.insert(*coupling_intersection_it);
                        }
                    }
                    coupling_intersections_.push_back(coupling_intersection_set);
                }
            }
        }
        // introduce a local to global map
        for (auto id = inside_elements_ids_.begin(); id != inside_elements_ids_.end(); ++id) {
            for (auto el = local_inside_grid_->leaf_view().template begin<0>();
                 el != local_inside_grid_->leaf_view().template end<0>(); ++el) {
                if (local_inside_grid_->leaf_view().indexSet().index(*el) == *id) {
                    // This is the inside element we are searching for..
                    local_to_inside_indices_.push_back({*id, local_to_inside_indices_.size()});
                }
            }
        }
    }
    template<int cd>  // fake template
    typename std::vector<LocalElementType>::iterator begin() { return inside_elements_.begin(); };

    template<int cd>  // fake template
    typename std::vector<LocalElementType>::iterator end() { return inside_elements_.end(); };

    typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::iterator ibegin(LocalElementType& inside_element) {
        return coupling_intersections_[local_to_inside_index(inside_element)].begin();
    };
    typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::iterator iend(LocalElementType& inside_element) {
        return coupling_intersections_[local_to_inside_index(inside_element)].end();
    };
private:
    size_t local_to_inside_index(const LocalElementType& local_element) const
    {
        auto id = local_inside_grid_->leaf_view().indexSet().index(local_element);
        for (auto index_pair = local_to_inside_indices_.begin() ; index_pair != local_to_inside_indices_.end(); ++index_pair) {
            if (id == index_pair->first) {
                return index_pair->second;
            }
        }
    }

private:
    const MacroElementType& inside_element_;
    const MacroElementType& outside_element_;
    const GridGlueType& dd_grid_;
    const MacroGridViewType& macro_grid_view_;
    const std::unique_ptr<LocalGridProviderType> local_inside_grid_;
    std::vector<LocalElementType> inside_elements_;
    std::vector<int> inside_elements_ids_;
    std::vector<std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>> coupling_intersections_;
    std::vector<std::pair<size_t, size_t>> local_to_inside_indices_;
};


GTEST_TEST(empty, main) {
    using GridType = YASP_2D_EQUIDISTANT_OFFSET;
    using ElementType = typename GridType::template Codim<0>::Entity;
    static constexpr size_t d = GridType::dimension;
    using D = typename GridType::ctype;
    using GridProviderType = Dune::XT::Grid::GridProvider<GridType>;
    using GridViewType = typename GridProviderType::LeafGridViewType;

    auto grid = Dune::XT::Grid::make_cube_grid<GridType>(0., 1., 2);
    std::cout << grid.dimDomain << std::endl;
    auto dd_grid = std::make_unique<Dune::XT::Grid::DD::Glued<GridType,GridType,Dune::XT::Grid::Layers::leaf>>(grid, 2);
    std::cout << dd_grid->num_subdomains() << std::endl;
    auto local_grid = dd_grid->local_grid(0);
    std::cout << local_grid->dimDomain << std::endl;
    std::cout << dd_grid->num_subdomains() << std::endl;

    GridViewType mgv = dd_grid->macro_grid_view();

    std::vector<ElementType> elements_;
    for (auto&& macro_element : Dune::elements(mgv)) {
      elements_.push_back(macro_element);
      const auto subdomain_id = dd_grid->subdomain(macro_element);
      std::cout << "id of the element " << subdomain_id << std::endl;
    }

    Iterator<ElementType> iterator(elements_);
    std::cout << (iterator == iterator) << std::endl;

    // iterator into a class like in the book !!
    ElementVector<ElementType> vector(elements_);

    using GridGlueType = Dune::XT::Grid::DD::Glued<GridType,GridType,Dune::XT::Grid::Layers::leaf>;
    using GlueType = typename GridGlueType::GlueType;
    using IntersectionType = typename GridViewType::Intersection;

    std::vector<std::set<IntersectionType, CompareType<IntersectionType>>> intersections_;

    for(auto el = vector.begin(); el != vector.end(); ++el) {
        std::cout << "center of element: " << (*el).geometry().center() << "\n";
        std::set<IntersectionType, CompareType<IntersectionType>> intersection_set;
        for (auto& macro_intersection : Dune::intersections(mgv, (*el))) {
            auto outside_element = macro_intersection.outside();
            std::cout << "outside element " << outside_element.geometry().center() << std::endl;
            intersection_set.insert(macro_intersection);
            std::cout << "elements of the set: {" << std::endl;
            for (auto it : intersection_set)
                std::cout << it.geometry().center() << std::endl;
            std::cout << "}" << std::endl;
        }
        intersections_.push_back(intersection_set);
    }

    ElementVectorAndIntersections<ElementType, IntersectionType, GridViewType> element_and_intersections(mgv, elements_, intersections_);

    for(auto el = element_and_intersections.begin(); el != element_and_intersections.end(); ++el) {
        std::cout << "center of element: " << (*el).geometry().center() << "\n";
        for(auto it = element_and_intersections.ibegin(*el); it != element_and_intersections.iend(*el); ++it) {
            std::cout << "center of intersection: " << (*it).geometry().center() << "\n";
        }
    }

    std::cout << "_____ COUPLING INTERSECTIONS _____" << std::endl;
    using CouplingIntersectionType = typename GlueType::Intersection;
    using CorrectedCouplingIntersectionType = CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, IntersectionType>;

    std::vector<std::set<CorrectedCouplingIntersectionType, CompareType<CorrectedCouplingIntersectionType>>> coupling_intersections_;
    std::set<CorrectedCouplingIntersectionType, CompareType<CorrectedCouplingIntersectionType>> coupling_intersection_set;
    // set does not work. inserting does not work with the compare type. TODO: CHECK THIS AGAIN.
//    std::vector<std::vector<CorrectedCouplingIntersectionType>> coupling_intersections_;

    for(auto el = element_and_intersections.begin(); el != element_and_intersections.end(); ++el) {
//        std::vector<CorrectedCouplingIntersectionType> coupling_intersection_set;
        for(auto it = element_and_intersections.ibegin(*el); it != element_and_intersections.iend(*el); ++it) {
//            std::cout << "boundary? " << it->boundary() << std::endl;
            auto inside_element = it->inside();
            auto outside_element = it->outside();
            std::cout << "center of inside element: " << inside_element.geometry().center() << "\n";
            std::cout << "center of outside element: " << outside_element.geometry().center() << "\n";

            if(it->boundary()) {
                std::cout << "skip this element, it is not inside the domain" << std::endl;
                continue;
            }
            // NOTE: shouldn't coupling() throw an assertion if the outside_element is not inside the domain as it is for boundary intersection???
            // The above if statement is a workaround for that
            auto& coupling = dd_grid->coupling(inside_element, -1, outside_element, -1, true);
            for (auto coupling_intersection_it = coupling.template ibegin<0>();
                 coupling_intersection_it != coupling.template iend<0>();
                 ++coupling_intersection_it) {
                std::cout << "current center is " << coupling_intersection_it->geometry().center() << std::endl;
//                auto coupling_intersection = std::make_unique<CorrectedCouplingIntersectionType>(*coupling_intersection_it, *it);
                CorrectedCouplingIntersectionType coupling_intersection(*coupling_intersection_it, *it);
                coupling_intersection_set.insert(coupling_intersection);
//                coupling_intersection_set.push_back(coupling_intersection);
            }
        }
       coupling_intersections_.push_back(coupling_intersection_set);
       coupling_intersection_set.clear();
    }

    std::cout << (*next(intersections_.at(3).begin(),1)).geometry().center() << std::endl;
    std::cout << (*next(coupling_intersections_.at(0).begin(),0)).inside().geometry().center() << std::endl;

    std::vector<std::vector<ElementType>> inside_elements;
    ElementVectorAndCouplingIntersections<ElementType, CorrectedCouplingIntersectionType, GridViewType> element_and_coupling_intersections(mgv, elements_, coupling_intersections_);

    for(auto el = element_and_coupling_intersections.begin(); el != element_and_coupling_intersections.end(); ++el) {
        std::cout << "center of element: " << (*el).geometry().center() << "\n";
        std::vector<ElementType> inside_elements_on_element;
        for(auto it = element_and_coupling_intersections.ibegin(*el); it != element_and_coupling_intersections.iend(*el); ++it) {
            std::cout << "center of coupling intersection: " << (*it).geometry().center() << "\n";
            auto inside_element = it->inside();
            inside_elements_on_element.push_back(inside_element);
//            std::cout << "center of micro inside element: " << inside_element.geometry().center() << std::endl;
        }
        inside_elements.push_back(inside_elements_on_element);
    }

    InsideElementsAndCouplingIntersections<ElementType, CorrectedCouplingIntersectionType, GridViewType> inside_element_and_intersections(mgv, elements_, coupling_intersections_, inside_elements);
    for(auto el = inside_element_and_intersections.begin(); el != inside_element_and_intersections.end(); ++el) {
        std::cout << "center of element: " << (*el).geometry().center() << "\n";
        for(auto it = inside_element_and_intersections.ebegin(*el); it != inside_element_and_intersections.eend(*el); ++it) {
            std::cout << "center of inside element: " << (*it).geometry().center() << "\n";
            }
        }

    TentativeGridView<GridGlueType> tentative_grid_view(*dd_grid);
    for(auto m_el = tentative_grid_view.mbegin(); m_el != tentative_grid_view.mend(); ++m_el) {
        std::cout << "center of macro element: " << (*m_el).geometry().center() << "\n";
        for(auto m_it = tentative_grid_view.mibegin(*m_el); m_it != tentative_grid_view.miend(*m_el); ++m_it) {
            std::cout << "center of macro intersection: " << (*m_it).geometry().center() << "\n";
        }
    }

    // start to work on Coupling grid view for nn and ss
    for (auto&& macro_element : Dune::elements(mgv)) {
        for (auto& macro_intersection : Dune::intersections(mgv, (macro_element))) {
            if (macro_intersection.boundary()) {
                std::cout << "skip this intersection, it is not inside the domain" << std::endl;
                continue;
            }
            auto inside_element = macro_intersection.inside();
            auto outside_element = macro_intersection.outside();
            TentativeCouplingGridView<GridGlueType> tentative_coupling_grid_view(inside_element,
                                                                                 outside_element,
                                                                                 *dd_grid);
            for (auto el = tentative_coupling_grid_view.template begin<0>(); el != tentative_coupling_grid_view.template end<0>(); ++el) {
                std::cout << "center of inside element is: " << el->geometry().center() << std::endl;
                for (auto intersection = tentative_coupling_grid_view.ibegin(*el);
                     intersection != tentative_coupling_grid_view.iend(*el); ++intersection) {
                    std::cout << "center of the corresponding intersection: " << intersection->geometry().center() << std::endl;
                }
            }
        }
    }

    std::cout << " ____ NOW WE TRY TO PUT IT INTO A GRIDVIEW INTERFACE ____" << std::endl;
    using CouplingGridViewType = Dune::XT::Grid::CouplingGridView<GridGlueType>;

    // use the gridview from view/new_coupling.hh
    for (auto&& macro_element : Dune::elements(mgv)) {
        for (auto& macro_intersection : Dune::intersections(mgv, macro_element)) {
            printf("_______\n");
            if (macro_intersection.boundary()) {
                std::cout << "skip this intersection, it is not inside the domain" << std::endl;
                continue;
            }
            auto inside_element = macro_intersection.inside();
            auto outside_element = macro_intersection.outside();
//            auto pgv = Dune::XT::Grid::make_periodic_grid_view<GridViewType>(mgv);
            auto cgv = Dune::XT::Grid::make_coupling_grid_view<GridViewType, GridGlueType, ElementType>(mgv, inside_element, outside_element, *dd_grid);
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

            auto functor = Dune::XT::Grid::GenericElementAndIntersectionFunctor<CouplingGridViewType>([] {},
                [](const auto& element) {
                std::cout << "WALKER: Element : " << element.geometry().center() << std::endl; },
                [](const auto& intersection, const ElementType&, const ElementType&) {
                std::cout << "WALKER:      Intersection : " << intersection.geometry().center() << std::endl; },
                [] {});

            // these lines are not required !!! but they test the CouplingGridProvider
            Dune::XT::Grid::CouplingGridProvider<CouplingGridViewType> coupling_provider(cgv);
            auto cgv_c = coupling_provider.coupling_view();

            auto walker = Dune::XT::Grid::Walker<CouplingGridViewType>(cgv_c); // or cgv itself
            walker.append(functor);
            walker.walk();
        }
    }


}
