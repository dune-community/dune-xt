// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2015 - 2019)
//   Tobias Leibner  (2015 - 2018, 2020)

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

template <class type>
struct CompareType
{
  bool operator()(const type& one, const type& other) const
  {
    return (std::addressof(one) < std::addressof(other));
  }
};

namespace Dune::XT::Grid {
namespace internal {

// forward for Traits
template <class BaseGridViewImp, class GridGlueType>
class CouplingGridViewWrapper;

//! Traits for CouplingGridViewWrapper
template <class BaseGridViewImp, class GridGlueType>
class CouplingGridViewWrapperTraits : public BaseGridViewImp::Traits
{
public:
  using BaseGridViewType = BaseGridViewImp;
  // use types from BaseGridViewType...
  using GridViewImp = CouplingGridViewWrapper<BaseGridViewType, GridGlueType>;
  using Grid = extract_grid_t<BaseGridViewType>;
  using IndexSet = typename BaseGridViewType::IndexSet;
  using CollectiveCommunication = extract_collective_communication_t<BaseGridViewType>;
  using BaseGridViewTraits = typename BaseGridViewType::Traits;
  using LocalElementType = typename GridGlueType::MicroEntityType;

  using GlueType = typename GridGlueType::GlueType;
  using CouplingIntersectionType = typename GlueType::Intersection;
  using Intersection = CouplingIntersectionType;
  using IntersectionIterator = typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::const_iterator;

  template <int cd>
  struct Codim : public BaseGridViewTraits::template Codim<cd>
  {
    // We need to define these in case BaseGridViewImp is a grid part.
    using Entity = extract_entity_t<BaseGridViewImp, cd>;
    using Geometry = extract_geometry_t<BaseGridViewImp, cd>;
    using LocalGeometry = extract_local_geometry_t<BaseGridViewImp, cd>;

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
template <class BaseGridViewImp, class GridGlueType>
class CouplingGridViewWrapper : public BaseGridViewImp
{
  using BaseType = BaseGridViewImp;
  using ThisType = CouplingGridViewWrapper;
  using Traits = CouplingGridViewWrapperTraits<BaseType, GridGlueType>;

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

//  TODO: add the macro intersection to use CorrectedCouplingIntersection
//  using CorrectedCouplingIntersectionType = CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroIntersectionType>;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {};

public:
  CouplingGridViewWrapper(const BaseType& base_grid_view,
                          const MacroElementType& ss,
                          const MacroElementType& nn,
                          GridGlueType& dd_grid)
    : BaseType(base_grid_view),
      inside_element_(ss),
      outside_element_(nn),
      dd_grid_(dd_grid),
      macro_grid_view_(dd_grid.macro_grid_view()),
      local_inside_grid_(dd_grid.local_grid(ss))

  {
    if (base_grid_view.comm().size() > 1)
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

  size_t local_to_inside_index(const LocalElementType& local_element) const
  {
      auto id = local_inside_grid_->leaf_view().indexSet().index(local_element);
      for (auto index_pair = local_to_inside_indices_.begin() ; index_pair != local_to_inside_indices_.end(); ++index_pair) {
          if (id == index_pair->first) {
              return index_pair->second;
          }
      }
  }

  int size(int codim) const
  {
    DUNE_THROW(NotImplemented,"");
  }

  int size(const Dune::GeometryType& type) const
  {
    DUNE_THROW(NotImplemented,"");
  }

  // define IntersectionIterator type
  typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::const_iterator ibegin(const LocalElementType& inside_element) const
  {
      return coupling_intersections_[local_to_inside_index(inside_element)].begin();
  };

  // define IntersectionIterator type
  typename std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>::const_iterator iend(const LocalElementType& inside_element) const
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
    DUNE_THROW(NotImplemented,"");
  }

private:
    const MacroElementType& inside_element_;
    const MacroElementType& outside_element_;
    GridGlueType& dd_grid_;
    const MacroGridViewType& macro_grid_view_;
    const std::shared_ptr<LocalGridProviderType> local_inside_grid_;
    std::vector<LocalElementType> inside_elements_;
    std::vector<int> inside_elements_ids_;
    std::vector<std::set<CouplingIntersectionType, CompareType<CouplingIntersectionType>>> coupling_intersections_;
    std::vector<std::pair<size_t, size_t>> local_to_inside_indices_;
}; // ... class CouplingGridViewWrapper ...


} // namespace internal


template <class BaseGridViewImp, class GridGlueImp>
class CouplingGridView
  : XT::Common::StorageProvider<internal::CouplingGridViewWrapper<BaseGridViewImp, GridGlueImp>>
  , public Dune::GridView<internal::CouplingGridViewWrapperTraits<BaseGridViewImp, GridGlueImp>>
{
  static_assert(is_view<BaseGridViewImp>::value);
  using Implementation = internal::CouplingGridViewWrapper<BaseGridViewImp, GridGlueImp>;
  using ImplementationStorage = typename XT::Common::StorageProvider<Implementation>;
  using BaseType = Dune::GridView<internal::CouplingGridViewWrapperTraits<BaseGridViewImp, GridGlueImp>>;

public:
  using GridGlueType = GridGlueImp;
  using BaseType::dimension;
  using BaseGridViewType = BaseGridViewImp;
  using MacroGridType = typename GridGlueType::MacroGridType;
  using MacroElementType = typename MacroGridType::template Codim<0>::Entity;

  CouplingGridView(const BaseGridViewType& base_grid_view,
                   const MacroElementType& ss,
                   const MacroElementType& nn,
                   GridGlueType& dd_grid)
    : ImplementationStorage(new Implementation(base_grid_view, ss, nn, dd_grid))
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


//template <bool codim_iters_provided, class GL>
//CouplingGridView<GL, codim_iters_provided>
//make_periodic_grid_view(const GL& base_grid_view,
//                        const std::bitset<GL::dimension> periodic_directions = std::bitset<GL::dimension>{}.set())
//{
//  return CouplingGridView<GL, codim_iters_provided>(base_grid_view, periodic_directions);
//}

template <class MG, class GT, class E>
CouplingGridView<MG, GT>
make_coupling_grid_view(const MG& base_grid_view, const E& ss, const E& nn, GT& dd_grid)
{
  return CouplingGridView<MG, GT>(base_grid_view, ss, nn, dd_grid);
}

} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_VIEW_PERIODIC_HH
