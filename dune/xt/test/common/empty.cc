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

template <class Intersection>
struct CompareType
{
  bool operator()(const Intersection& one, const Intersection& other) const
  {
//    std::cout << "one center " << one.geometry().center() << std::endl;
//    std::cout << "other center " << other.geometry().center() << std::endl;
    auto decision = Dune::XT::Common::FloatCmp::ne(one.geometry().center(), other.geometry().center());
//    if (!decision)
//        std::cout << "I do not add this" << std::endl;
//    auto decision = Dune::XT::Common::FloatCmp::lt(one.geometry().center(), other.geometry().center());
//    std::cout << "my decision is " << decision << std::endl;
    return decision;
  }
};

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
//    class Iterator {
//    public:
//        Iterator(E* init=nullptr)
//            : current_element_{init} {}
//        bool operator==(const Iterator&) const {
//            std::cout << "hello" << std::endl;
//            return true;
//        }
//        bool operator!=(const Iterator&) const;
//        Iterator& operator++();
//        Iterator operator++(int);
//        E & operator*() const;
//        E & operator->() const;
//    private:
//        E* current_element_;
//    };
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

//template<class GluedType>
//class GluedView {
//    using GridType = typename GluedType::GridType;
//    using ElementType = typename GridType::template Codim<0>::Entity;
//public:
//    ElementVector(GluedType& dd_grid)
//        : vector_{vector} {
//        std::vector<ElementType>& vector_;
//    }
//    class Iterator {
//    public:
//        Iterator(E* init=nullptr)
//            : current_element_{init} {}
//        bool operator==(const Iterator&) const {
//            std::cout << "hello" << std::endl;
//            return true;
//        }
//        bool operator!=(const Iterator&) const;
//        Iterator& operator++();
//        Iterator operator++(int);
//        E & operator*() const;
//        E & operator->() const;
//    private:
//        E* current_element_;
//    };
//    typename std::vector<E>::iterator begin() const { return vector_.begin(); };
//    typename std::vector<E>::iterator end() const { return vector_.end(); };
//private:
//    std::vector<E>& vector_;
//};


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

    using GlueType = typename Dune::XT::Grid::DD::Glued<GridType,GridType,Dune::XT::Grid::Layers::leaf>::GlueType;
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

    for(auto el = element_and_intersections.begin(); el != element_and_intersections.end(); ++el) {
        for(auto it = element_and_intersections.ibegin(*el); it != element_and_intersections.iend(*el); ++it) {
//            std::cout << "boundary? " << it->boundary() << std::endl;
            auto inside_element = it->inside();
            auto outside_element = it->outside();
            std::cout << "center of inside element: " << inside_element.geometry().center() << "\n";
            std::cout << "center of outside element: " << outside_element.geometry().center() << "\n";
//            auto local_grid = dd_grid->local_grid(outside_element);

            if(it->boundary()) {
                std::cout << "skip this element, it is not inside the domain" << std::endl;
                continue;
            }
            // NOTE: shouldn't coupling() throw an assertion if the outside_element is not inside the domain as it is for boundary intersection???
            // The above if statement is a workaround for that
            auto coupling = dd_grid->coupling(inside_element, -1, outside_element, -1, true);
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
    int a = 0;
}
