// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tim Keil (2020 - 2021)

#ifndef DUNE_XT_GRID_VIEW_NEW_COUPLING_HH
#define DUNE_XT_GRID_VIEW_NEW_COUPLING_HH

#include <bitset>
#include <iterator>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <dune/geometry/typeindex.hh>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/search.hh>
#include <dune/xt/grid/type_traits.hh>


namespace Dune::XT::Grid {
namespace internal {

/**
 * Used for comparing Intersections in sets. Adress uniquely defines interesections,
 * other properties were not successful for comparison.
 */
template <class type>
struct CompareType
{
  bool operator()(const type& one, const type& other) const
  {
    return (std::addressof(one) < std::addressof(other));
  }
};

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
  using Entity = typename MacroIntersectionType::Entity;
  static constexpr auto dimensionworld = MacroIntersectionType::dimensionworld;
  static constexpr auto dimension = MacroIntersectionType::dimension;
  using typename BaseType::GlobalCoordinate;
  using typename BaseType::LocalCoordinate;

  CouplingIntersectionWithCorrectNormal(const CouplingIntersectionType& coupling_intersection,
                                        const MacroIntersectionType macro_intersection)
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

  bool boundary() const
  {
    return this->macro_intersection_.boundary();
  }

  size_t boundarySegmentIndex() const
  {
    return this->macro_intersection_.boundarySegmentIndex();
  }

  //! Compares two intersections for equality. copied from dune/grid/common/intersection
  bool operator==(const CouplingIntersectionWithCorrectNormal& other) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  //! Compares two intersections for inequality. copied from dune/grid/common/intersection
  bool operator!=(const CouplingIntersectionWithCorrectNormal& other) const
  {
    DUNE_THROW(NotImplemented, "");
  }

private:
  const MacroIntersectionType macro_intersection_;
  mutable GlobalCoordinate global_;
  mutable LocalCoordinate local_;
}; // class CouplingIntersectionWithCorrectNormal


// forward for Traits
template <class GridGlueImp>
class CouplingGridViewWrapper;

//! Traits for CouplingGridViewWrapper
template <class GridGlueImp>
class CouplingGridViewWrapperTraits : public GridGlueImp::MacroGridViewType::Traits
{
public:
  using BaseGridViewType = typename GridGlueImp::MacroGridViewType;
  using GridGlueType = GridGlueImp;
  // use types from BaseGridViewType...
  using GridViewImp = CouplingGridViewWrapper<GridGlueType>;
  using Grid = extract_grid_t<BaseGridViewType>;
  using IndexSet = typename BaseGridViewType::IndexSet;
  using CollectiveCommunication = extract_collective_communication_t<BaseGridViewType>;
  using BaseGridViewTraits = typename BaseGridViewType::Traits;
  using LocalElementType = typename GridGlueType::MicroEntityType;

  using GlueType = typename GridGlueType::GlueType;
  using CouplingIntersectionType = typename GlueType::Intersection;
  using MacroInterSectionType = typename BaseGridViewTraits::Intersection;
  using CorrectedCouplingIntersectionType =
      CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroInterSectionType>;
  using Intersection = CorrectedCouplingIntersectionType;
  using IntersectionIterator = typename std::set<CorrectedCouplingIntersectionType,
                                                 CompareType<CorrectedCouplingIntersectionType>>::const_iterator;

  template <int cd>
  struct Codim : public BaseGridViewTraits::template Codim<cd>
  {
    // We need to define these in case BaseGridViewImp is a grid part.
    using Entity = extract_entity_t<BaseGridViewType, cd>;
    using Geometry = extract_geometry_t<BaseGridViewType, cd>;
    using LocalGeometry = extract_local_geometry_t<BaseGridViewType, cd>;

    using Iterator = typename std::vector<LocalElementType>::const_iterator;

    template <PartitionIteratorType pit>
    struct Partition : public BaseGridViewTraits::template Codim<cd>::template Partition<pit>
    {
      using Iterator = typename std::vector<LocalElementType>::const_iterator;
    }; // struct Partition
  }; // ... struct Codim ...

  enum
  {
    conforming = BaseGridViewTraits::conforming
  };
  enum
  {
    dimension = Grid::dimension
  };
  enum
  {
    dimensionworld = Grid::dimensionworld
  };

  using ctype = typename Grid::ctype;

}; // ... class CouplingGridViewWrapperTraits ...


/** \brief Actual Implementation of CouplingGridView
 *  \see CouplingGridView
 */
template <class GridGlueType>
class CouplingGridViewWrapper : public GridGlueType::MacroGridViewType
{
  using BaseType = typename GridGlueType::MacroGridViewType;
  using ThisType = CouplingGridViewWrapper;
  using Traits = CouplingGridViewWrapperTraits<GridGlueType>;

public:
  using GlueType = typename GridGlueType::GlueType;
  using CouplingIntersectionType = typename GlueType::Intersection;
  using IntersectionType = CouplingIntersectionType;

  using Grid = extract_grid_t<BaseType>;
  using IndexSet = typename BaseType::IndexSet;
  using ElementType = extract_entity_t<BaseType, 0>;
  using IndexType = typename Traits::IndexSet::IndexType;

  using MacroGridViewType = typename GridGlueType::MacroGridViewType;
  using MacroGridType = typename GridGlueType::MacroGridType;
  using MacroElementType = typename MacroGridType::template Codim<0>::Entity;
  using MacroIntersectionType = typename MacroGridViewType::Intersection;

  using LocalGridProviderType = typename GridGlueType::LocalGridProviderType;
  using LocalGridType = typename GridGlueType::LocalGridType;
  using LocalGridViewType = typename GridGlueType::MicroGridViewType;
  using LocalElementType = typename GridGlueType::MicroEntityType;

  using IntersectionIterator = typename Traits::IntersectionIterator;
  using CorrectedCouplingIntersectionType =
      CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroIntersectionType>;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {};

public:
  CouplingGridViewWrapper(const MacroElementType& ss,
                          const MacroElementType& nn,
                          GridGlueType& dd_grid,
                          const MacroIntersectionType macro_intersection)
    : BaseType(dd_grid.macro_grid_view())
    , inside_element_(ss)
    , outside_element_(nn)
    , dd_grid_(dd_grid)
    , macro_intersection_(macro_intersection)
    , macro_grid_view_(dd_grid.macro_grid_view())
    , local_inside_grid_(dd_grid.local_grid(ss))

  {
    if (macro_grid_view_.comm().size() > 1)
      DUNE_THROW(Dune::NotImplemented, "CouplingGridView does not work for MPI-distributed grids!");
    this->update();
  } // constructor CouplingGridViewWrapper(...)

  CouplingGridViewWrapper(const ThisType& other) = default;
  CouplingGridViewWrapper(ThisType&& other) = default;
  ThisType& operator=(const ThisType& other) = default;
  ThisType& operator=(ThisType&& other) = default;

  // This is the only member function that actually changes the state of the CouplingGridView. The copy constructor and
  // copy assignment operators only do a shallow copy of the shared_ptrs, so we have to make sure that we reassign the
  // shared_ptrs in this function to avoid changing state of a copied-from CouplingGridView.
  void update()
  {
    auto& coupling = dd_grid_.coupling(inside_element_, -1, outside_element_, -1, true);
    // Note: This iterator iterates over all interesctions in the coupling grid. Not only over all intersections w.r.t.
    //       a single element
    for (auto coupling_intersection_it = coupling.template ibegin<0>();
         coupling_intersection_it != coupling.template iend<0>();
         ++coupling_intersection_it) {
      auto inside = coupling_intersection_it->inside();
      inside_elements_ids_.push_back(local_inside_grid_.leaf_view().indexSet().index(inside));
    }
    // some coupling intersection may have the same inside element, remove duplicates
    std::sort(inside_elements_ids_.begin(), inside_elements_ids_.end());
    auto last = std::unique(inside_elements_ids_.begin(), inside_elements_ids_.end());
    inside_elements_ids_.erase(last, inside_elements_ids_.end());

    for (auto&& id : inside_elements_ids_) {
      for (auto&& el : elements(local_inside_grid_.leaf_view())) {
        if (local_inside_grid_.leaf_view().indexSet().index(el) == id) {
          // This is the inside element we are searching for... add it to the vector
          inside_elements_.push_back(el);
          std::set<CorrectedCouplingIntersectionType, CompareType<CorrectedCouplingIntersectionType>>
              coupling_intersection_set;
          // now iteratate over all intersections to find all coupling intersections
          for (auto coupling_intersection_it = coupling.template ibegin<0>();
               coupling_intersection_it != coupling.template iend<0>();
               ++coupling_intersection_it) {
            auto inside = coupling_intersection_it->inside();
            auto inside_id = local_inside_grid_.leaf_view().indexSet().index(inside);
            if (inside_id == id) {
              CorrectedCouplingIntersectionType coupling_intersection(*coupling_intersection_it, macro_intersection_);
              coupling_intersection_set.insert(coupling_intersection);
            }
          }
          coupling_intersections_.push_back(coupling_intersection_set);
        }
      }
    }
    // introduce a local to global map
    for (auto&& id : inside_elements_ids_) {
      for (auto&& el : elements(local_inside_grid_.leaf_view())) {
        if (local_inside_grid_.leaf_view().indexSet().index(el) == id) {
          // This is the inside element we are searching for..
          local_to_inside_indices_.push_back({id, local_to_inside_indices_.size()});
        }
      }
    }
  }

  size_t local_to_inside_index(const LocalElementType& local_element) const
  {
    auto id = local_inside_grid_.leaf_view().indexSet().index(local_element);
    for (auto&& index_pair : local_to_inside_indices_) {
      if (id == index_pair.first) {
        return index_pair.second;
      }
    }
  }

  int size(int codim) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  int size(const Dune::GeometryType& type) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  IntersectionIterator ibegin(const LocalElementType& inside_element) const
  {
    return coupling_intersections_[local_to_inside_index(inside_element)].begin();
  };

  IntersectionIterator iend(const LocalElementType& inside_element) const
  {
    return coupling_intersections_[local_to_inside_index(inside_element)].end();
  };

  template <int cd>
  typename Codim<cd>::Iterator begin() const
  {
    return inside_elements_.begin();
  }

  template <int cd>
  typename Codim<cd>::Iterator end() const
  {
    return inside_elements_.end();
  }

  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin() const
  {
    return inside_elements_.begin();
  }

  template <int cd, PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end() const
  {
    return inside_elements_.end();
  }

  const IndexSet& indexSet() const
  {
    DUNE_THROW(NotImplemented, "");
  }

private:
  const MacroElementType& inside_element_;
  const MacroElementType& outside_element_;
  GridGlueType& dd_grid_;
  const MacroIntersectionType macro_intersection_;
  const MacroGridViewType& macro_grid_view_;
  const LocalGridProviderType& local_inside_grid_;
  std::vector<LocalElementType> inside_elements_;
  std::vector<int> inside_elements_ids_;
  std::vector<std::set<CorrectedCouplingIntersectionType, CompareType<CorrectedCouplingIntersectionType>>>
      coupling_intersections_;
  std::vector<std::pair<size_t, size_t>> local_to_inside_indices_;
}; // ... class CouplingGridViewWrapper ...


} // namespace internal


template <class GridGlueImp>
class CouplingGridView
  : XT::Common::StorageProvider<internal::CouplingGridViewWrapper<GridGlueImp>>
  , public Dune::GridView<internal::CouplingGridViewWrapperTraits<GridGlueImp>>
{
  static_assert(is_view<typename GridGlueImp::MacroGridViewType>::value);
  using Implementation = internal::CouplingGridViewWrapper<GridGlueImp>;
  using ImplementationStorage = typename XT::Common::StorageProvider<Implementation>;
  using BaseType = Dune::GridView<internal::CouplingGridViewWrapperTraits<GridGlueImp>>;

public:
  using GridGlueType = GridGlueImp;
  using BaseType::dimension;
  using BaseGridViewType = typename GridGlueImp::MacroGridViewType;
  using MacroGridType = typename GridGlueType::MacroGridType;
  using MacroElementType = typename MacroGridType::template Codim<0>::Entity;
  using MacroIntersectionType = typename BaseGridViewType::Traits::Intersection;

  CouplingGridView(const MacroElementType& ss,
                   const MacroElementType& nn,
                   GridGlueType& dd_grid,
                   const MacroIntersectionType macro_intersection)
    : ImplementationStorage(new Implementation(ss, nn, dd_grid, macro_intersection))
    , BaseType(ImplementationStorage::access())
  {}

  CouplingGridView(const CouplingGridView& other)
    : ImplementationStorage(new Implementation(other.access()))
    , BaseType(ImplementationStorage::access())
  {}

  void update()
  {
    ImplementationStorage::access().update();
  }
}; // class CouplingGridView


template <class GT, class E, class IT>
static CouplingGridView<GT> make_coupling_grid_view(const E& ss, const E& nn, GT& dd_grid, const IT& macro_intersection)
{
  return CouplingGridView<GT>(ss, nn, dd_grid, macro_intersection);
}

} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_VIEW_PERIODIC_HH
