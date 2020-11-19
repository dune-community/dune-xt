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

#ifndef DUNE_XT_GRID_VIEW_COUPLING_HH
#define DUNE_XT_GRID_VIEW_COUPLING_HH

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

template<class GridGlueType>
class FirstCouplingGridView {
public:
    using MacroGridViewType = typename GridGlueType::MacroGridViewType;
    using MacroGridType = typename GridGlueType::MacroGridType;
    using MacroElementType = typename MacroGridType::template Codim<0>::Entity;
    using MacroIntersectionType = typename MacroGridViewType::Intersection;

    using LocalGridProviderType = typename GridGlueType::LocalGridProviderType;
    using LocalGridType = typename GridGlueType::LocalGridType;
    using LocalGridViewType = typename GridGlueType::MicroGridViewType;
    using LocalElementType = typename GridGlueType::MicroEntityType;

    using GlueType = typename GridGlueType::GlueType;
    using CouplingIntersectionType = typename GlueType::Intersection;
    // TODO: add the macro intersection to use CorrectedCouplingIntersection
//    using CorrectedCouplingIntersectionType = CouplingIntersectionWithCorrectNormal<CouplingIntersectionType, MacroIntersectionType>;

    enum { //! \brief The dimension of the grid
      dimension = LocalGridType :: dimension
    };

    enum { //! \brief The dimension of the world the grid lives in
      dimensionworld = LocalGridType :: dimensionworld
    };
public:
    FirstCouplingGridView(const MacroElementType& ss, const MacroElementType& nn, GridGlueType& dd_grid)  // <- insert this as a pointer?
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

    const LocalGridViewType &grid () const
    {
        DUNE_THROW(NotImplemented,"");
    }

    /** \brief obtain the index set
     *
     * The lifetime of the returned index set is bound to the lifetime of the
     * grid view. Keep a copy of the grid view to prevent the index set from
     * becoming a dangling reference.
     */
//    const IndexSet &indexSet () const
//    {
//        DUNE_THROW(NotImplemented,"");
//    }

    /** \brief obtain number of entities in a given codimension */
    int size ( int codim ) const
    {
        DUNE_THROW(NotImplemented,"");
    }

    /** \brief obtain number of entities with a given geometry type */
    int size ( const GeometryType &type ) const
    {
        DUNE_THROW(NotImplemented,"");
    }

    /** @brief Return true if the given entity is contained in this grid view
     * @todo Currently we call the implementation on the IndexSet.  This may lead to suboptimal efficiency.
     *
     * \note If the input element e is not an element of the grid, then
     *       the result of contains() is undefined.
     */
    template<class EntityType>
    bool contains (const EntityType& e) const
    {
        DUNE_THROW(NotImplemented,"");
    }


//    /** \brief obtain begin iterator for this view */
//    template< int cd >
//    typename Codim< cd > :: Iterator begin () const
//    {
//      return impl().template begin<cd>();
//    }

//    /** \brief obtain end iterator for this view */
//    template< int cd >
//    typename Codim< cd > :: Iterator end () const
//    {
//      return impl().template end<cd>();
//    }

//    /** \brief obtain begin iterator for this view */
//    template< int cd , PartitionIteratorType pitype >
//    typename Codim< cd > :: template Partition< pitype > :: Iterator
//    begin () const
//    {
//      return impl().template begin<cd,pitype>();
//    }

//    /** \brief obtain end iterator for this view */
//    template< int cd, PartitionIteratorType pitype >
//    typename Codim< cd > :: template Partition< pitype > :: Iterator
//    end () const
//    {
//      return impl().template end<cd,pitype>();
//    }

//    /** \brief obtain begin intersection iterator with respect to this view */
//    IntersectionIterator
//    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
//    {
//      return impl().ibegin(entity);
//    }

//    /** \brief obtain end intersection iterator with respect to this view */
//    IntersectionIterator
//    iend ( const typename Codim< 0 > :: Entity &entity ) const
//    {
//      return impl().iend(entity);
//    }

//    /** \brief obtain collective communication object */
//    const CollectiveCommunication &comm () const
//    {
//        DUNE_THROW(NotImplemented,"");
//    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
        DUNE_THROW(NotImplemented,"");
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const
    {
        DUNE_THROW(NotImplemented,"");
    }

    /** \brief Communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
        DUNE_THROW(NotImplemented,"");
    }
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

} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_VIEW_COUPLING_HH
