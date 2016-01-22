// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)
//   Rene Milk       (2015)
//   Tobias Leibner  (2015)

#ifndef DUNE_XT_GRID_PERIODICVIEW_HH
#define DUNE_XT_GRID_PERIODICVIEW_HH

#include <bitset>
#include <map>
#include <utility>
#include <vector>

#include <dune/grid/common/gridview.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/search.hh>

namespace Dune {
namespace XT {
namespace Grid {

#if HAVE_DUNE_GRID

namespace internal {

// forward
template <class RealGridViewImp, bool use_less_memory>
class PeriodicGridViewTraits;

template <class RealGridViewImp, bool use_less_memory>
class PeriodicIndexSet : public Dune::IndexSet< typename RealGridViewImp::Grid, PeriodicIndexSet<RealGridViewImp, use_less_memory>, typename RealGridViewImp::IndexSet::IndexType, typename RealGridViewImp::IndexSet::Types >
{
    typedef RealGridViewImp RealGridViewType;
    typedef PeriodicIndexSet<RealGridViewType, use_less_memory> ThisType;
    typedef typename RealGridViewType::IndexSet RealIndexSetType;
    typedef typename Dune::IndexSet< typename RealGridViewType::Grid, ThisType, typename RealGridViewType::IndexSet::IndexType, typename RealGridViewType::IndexSet::Types > BaseType;
    typedef PeriodicGridViewTraits< RealGridViewImp, use_less_memory > Traits;
public:
    typedef typename RealIndexSetType::IndexType IndexType;
    typedef typename RealIndexSetType::Types Types;
private:
    typedef std::map< IndexType, IndexType> IndexMapType;
public:
    static const int dimDomain = RealGridViewType::dimension;

    PeriodicIndexSet(const RealIndexSetType& real_index_set,
                     const std::vector< IndexMapType >& index_maps,
                     std::vector< size_t > entity_counts,
                     std::map< Dune::GeometryType, size_t > geometry_type_counts,
                     const std::vector< std::vector< IndexType > >& new_indices_vector)
        : BaseType()
        , real_index_set_(real_index_set)
        , index_maps_(index_maps)
        , entity_counts_(entity_counts)
        , geometry_type_counts_(geometry_type_counts)
        , new_indices_vector_(new_indices_vector)
    {
        assert(index_maps_.size() >= dimDomain + 1);
        assert(entity_counts_.size() >= dimDomain + 1);
        assert(new_indices_vector_.size() >= dimDomain + 1);
    }

    template< int cd >
    IndexType index(const typename Traits::template Codim< cd >::Entity &entity) const
    {
        IndexType real_entity_index = real_index_set_.template index<cd>(entity);
        if (cd == 0) {
            return real_entity_index;
        } else {
            if (!use_less_memory) {
                return new_indices_vector_[cd][real_entity_index];
             } else {
                auto iterator_to_index = index_maps_[cd].find(real_entity_index);
                if (iterator_to_index != index_maps_[cd].end())
                    return (*iterator_to_index).second;
                else
                    return real_entity_index;
            }
        }
    }

    template<class EntityType >
    IndexType index(const EntityType &entity) const
    {
        return index<EntityType::codimension>(entity);
    }

    template<int cd>
    IndexType subIndex(const typename Traits::template Codim< cd >::Entity &entity, int i, unsigned int codim) const
    {
        IndexType real_sub_index = real_index_set_.template subIndex<cd>(entity,i,codim);
        if (codim == 0) {
            return real_sub_index;
         } else {
            if (!use_less_memory) {
                return new_indices_vector_[codim][real_sub_index];
             } else {
                auto iterator_to_index = index_maps_[codim].find(real_sub_index);
                if (iterator_to_index != index_maps_[codim].end())
                    return (*iterator_to_index).second;
                else
                    return real_sub_index;
            }
        }
    }

    template<class EntityType>
    IndexType subIndex(const EntityType &entity, int i, unsigned int codim) const
    {
        return subIndex<EntityType::codimension>(entity,i,codim);
    }

    Types types(int codim) const
    {
        return real_index_set_.types(codim);
    }

    IndexType size(Dune::GeometryType type) const
    {
        return geometry_type_counts_.count(type) ? geometry_type_counts_.at(type) : 0;
    }

    IndexType size(int codim) const
    {
       assert(codim <= dimDomain);
       return entity_counts_[codim];
    }

    template<class EntityType >
    bool contains(const EntityType &entity) const
    {
        return real_index_set_.contains(entity);
    }

private:
    const RealIndexSetType& real_index_set_;
    const std::vector< IndexMapType > index_maps_;
    std::vector< size_t > entity_counts_;
    std::map< Dune::GeometryType, size_t > geometry_type_counts_;
    const std::vector< std::vector< IndexType > > new_indices_vector_;
}; // class PeriodicIndexSet<...>

/** \brief Intersection for PeriodicGridView
 *
 * PeriodicIntersection is derived from the Intersection of the underlying GridView. On the inside of the grid or if
 * periodic_ is false, the PeriodicIntersection will behave exactly like its BaseType. If periodic_ is true, the
 * PeriodicIntersection will return neighbor == true even if it actually is on the boundary. In this case, outside(),
 * geometryInOutside() and indexInOutside() are well-defined and give the information from the periodically adjacent
 * entity.
 *
 * \see PeriodicGridView
 */
template <class RealGridViewImp>
class PeriodicIntersection : public RealGridViewImp::Intersection
{
  typedef RealGridViewImp RealGridViewType;
  typedef PeriodicIntersection<RealGridViewType> ThisType;
  typedef typename RealGridViewType::Intersection BaseType;

public:
  using typename BaseType::LocalGeometry;
  typedef typename BaseType::Entity EntityType;
  typedef typename RealGridViewType::IntersectionIterator RealIntersectionIteratorType;
  static const size_t dimDomain = RealGridViewType::dimension;

  //! \brief Constructor from real intersection
  PeriodicIntersection(const BaseType& real_intersection, const RealGridViewType& real_grid_view,
                       const std::pair<bool, EntityType>& periodic_pair)
    : BaseType(real_intersection)
    , periodic_(periodic_pair.first)
    , outside_(periodic_pair.second)
    , real_grid_view_(new Common::ConstStorageProvider<RealGridViewType>(real_grid_view))
  {
  }

  // methods that differ from BaseType
  bool neighbor() const
  {
    if (periodic_)
      return true;
    else
      return BaseType::neighbor();
  } // bool neighbor() const

  EntityType outside() const
  {
    if (periodic_)
      return outside_;
    else
      return EntityType(BaseType::outside());
  } // ... outside() const

  LocalGeometry geometryInOutside() const
  {
    if (periodic_) {
      return find_intersection_in_outside().geometryInInside();
    } else {
      return BaseType::geometryInOutside();
    }
  } // ... geometryInOutside() const

  int indexInOutside() const
  {
    if (periodic_) {
      return find_intersection_in_outside().indexInInside();
    } else {
      return BaseType::indexInOutside();
    }
  } // int indexInOutside() const

private:
  // tries to find intersection in outside (works only if periodic_ == true)
  BaseType find_intersection_in_outside() const
  {
    const auto coords                                   = this->geometry().center();
    RealIntersectionIteratorType outside_i_it           = real_grid_view_->access().ibegin(outside_);
    const RealIntersectionIteratorType outside_i_it_end = real_grid_view_->access().iend(outside_);
    // walk over outside intersections and find an intersection on the boundary that differs only in one coordinate
    for (; outside_i_it != outside_i_it_end; ++outside_i_it) {
      const BaseType& curr_outside_intersection = *outside_i_it;
      if (curr_outside_intersection.boundary()) {
        const auto curr_outside_intersection_coords = curr_outside_intersection.geometry().center();
        size_t coord_difference_count = 0;
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          if (Dune::XT::Common::FloatCmp::ne(curr_outside_intersection_coords[ii], coords[ii])) {
            ++coord_difference_count;
          }
        }
        if (coord_difference_count == size_t(1)) {
          return *outside_i_it;
        }
      }
    }
    DUNE_THROW(Dune::InvalidStateException, "Could not find outside intersection!");
    return *(real_grid_view_->access().ibegin(outside_));
  } // ... find_intersection_in_outside() const

protected:
  bool periodic_;
  EntityType outside_;
  std::unique_ptr<Common::ConstStorageProvider<RealGridViewType>> real_grid_view_;
}; // ... class PeriodicIntersection ...

/** \brief IntersectionIterator for PeriodicGridView
 *
 * PeriodicIntersectionIterator is derived from the IntersectionIterator of the underlying GridView and behaves exactly
 * like the underlying IntersectionIterator except that it returns a PeriodicIntersection in its operator* and
 * operator-> methods.
 *
 * \see PeriodicGridView
 */
template <class RealGridViewImp>
class PeriodicIntersectionIterator : public RealGridViewImp::IntersectionIterator
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename RealGridViewType::IntersectionIterator BaseType;

public:
  typedef typename BaseType::Intersection RealIntersectionType;
  typedef int IntersectionIndexType;
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef typename RealGridViewType::template Codim<0>::Entity EntityType;
  typedef std::pair<bool, EntityType> PeriodicPairType;
  static const size_t dimDomain = RealGridViewType::dimension;

  PeriodicIntersectionIterator(BaseType real_intersection_iterator, const RealGridViewType& real_grid_view,
                               const EntityType& entity,
                               const std::map<IntersectionIndexType, PeriodicPairType>& intersection_map)
    : BaseType(real_intersection_iterator)
    , real_grid_view_(real_grid_view)
    , entity_(entity)
    , has_boundary_intersections_(entity_.hasBoundaryIntersections())
    , intersection_map_(intersection_map)
    , nonperiodic_pair_(std::make_pair(bool(false), EntityType(entity_)))
    , current_intersection_(create_current_intersection_safely())
  {
  }

  // methods that differ from BaseType
  const Intersection& operator*() const
  {
    current_intersection_ = create_current_intersection();
    return *current_intersection_;
  }

  const Intersection* operator->() const
  {
    current_intersection_ = create_current_intersection();
    return &(*current_intersection_);
  }

private:
  std::unique_ptr<Intersection> create_current_intersection() const
  {
    return Common::make_unique<Intersection>(BaseType::operator*(),
                                             real_grid_view_,
                                             has_boundary_intersections_
                                                 ? intersection_map_.at((BaseType::operator*()).indexInInside())
                                                 : (const PeriodicPairType&)nonperiodic_pair_);
  } // ... create_current_intersection() const

  std::unique_ptr<Intersection> create_current_intersection_safely() const
  {
    const bool is_iend                            = (*this == real_grid_view_.iend(entity_));
    const RealIntersectionType& real_intersection = is_iend ? *real_grid_view_.ibegin(entity_) : BaseType::operator*();
    return Common::make_unique<Intersection>(real_intersection,
                                             real_grid_view_,
                                             has_boundary_intersections_
                                                 ? intersection_map_.at(real_intersection.indexInInside())
                                                 : (const PeriodicPairType&)nonperiodic_pair_);
  } // ... create_current_intersection_safely() const

  const RealGridViewType& real_grid_view_;
  const EntityType& entity_;
  const bool has_boundary_intersections_;
  const std::map<IntersectionIndexType, PeriodicPairType>& intersection_map_;
  PeriodicPairType nonperiodic_pair_;
  mutable std::unique_ptr<Intersection> current_intersection_;
}; // ... class PeriodicIntersectionIterator ...

// forward
template <class RealGridViewImp, bool use_less_memory>
class PeriodicGridViewImp;

//! Traits for PeriodicGridView
template <class RealGridViewImp, bool use_less_memory>
class PeriodicGridViewTraits
{
public:
    typedef RealGridViewImp RealGridViewType;
    // use types from RealGridViewType...
    typedef PeriodicGridViewImp<RealGridViewType, use_less_memory> GridViewImp;
    typedef typename RealGridViewType::Grid Grid;
    typedef PeriodicIndexSet<RealGridViewType, use_less_memory> IndexSet;
    typedef typename RealGridViewType::CollectiveCommunication CollectiveCommunication;
    typedef typename RealGridViewType::Traits RealGridViewTraits;

    template <int cd>
    struct Codim
    {
        class PeriodicIterator : public RealGridViewTraits::template Codim<cd>::Iterator
        {
            typedef typename RealGridViewTraits::template Codim<cd>::Iterator BaseType;
            typedef PeriodicIterator ThisType;
            typedef typename RealGridViewType::IndexSet RealIndexSetType;

        public:
            typedef typename IndexSet::IndexType IndexType;

            PeriodicIterator(BaseType real_iterator,
                             const std::set< IndexType >* entities_to_skip,
                             const RealIndexSetType* real_index_set,
                             const BaseType real_it_end)
              : BaseType(real_iterator)
              , entities_to_skip_(entities_to_skip)
              , real_index_set_(real_index_set)
              , real_it_end_(real_it_end)
            {
            }

            // methods that differ from BaseType
            ThisType& operator++()
            {
                BaseType::operator++();
                while (cd > 0 && *this != real_it_end_ && entities_to_skip_->count(real_index_set_->index(this->operator*())))
                    BaseType::operator++();
                return *this;
            }

            ThisType operator++(int)
            {
               return this->operator++();
            }

        private:
            const std::set< IndexType >* entities_to_skip_;
            const RealIndexSetType* real_index_set_;
            const BaseType real_it_end_;
        };

        typedef PeriodicIterator Iterator;
        typedef typename RealGridViewTraits::template Codim<cd>::EntityPointer EntityPointer;
        typedef typename RealGridViewTraits::template Codim<cd>::Entity Entity;
        typedef typename RealGridViewTraits::template Codim<cd>::Geometry Geometry;
        typedef typename RealGridViewTraits::template Codim<cd>::LocalGeometry LocalGeometry;

        template <PartitionIteratorType pit>
        struct Partition
        {
            class PeriodicIterator : public RealGridViewTraits::template Codim<cd>::template Partition< pit >::Iterator
            {
                typedef typename RealGridViewTraits::template Codim<cd>::template Partition< pit >::Iterator BaseType;
                typedef PeriodicIterator ThisType;
                typedef typename RealGridViewType::IndexSet RealIndexSetType;

            public:
                typedef typename IndexSet::IndexType IndexType;

                PeriodicIterator(BaseType real_iterator,
                                 const std::set< IndexType >* entities_to_skip,
                                 const RealIndexSetType* real_index_set,
                                 const BaseType real_it_end)
                  : BaseType(real_iterator)
                  , entities_to_skip_(entities_to_skip)
                  , real_index_set_(real_index_set)
                  , real_it_end_(real_it_end)
                {
                }

                // methods that differ from BaseType
                ThisType& operator++()
                {
                    BaseType::operator++();
                    while (cd > 0 && *this != real_it_end_ && entities_to_skip_->count(real_index_set_->index(this->operator*())))
                        BaseType::operator++();
                    return *this;
                }

                ThisType operator++(int)
                {
                   return this->operator++();
                }

            private:
                const std::set< IndexType >* entities_to_skip_;
                const RealIndexSetType* real_index_set_;
                const BaseType real_it_end_;
            };
            typedef PeriodicIterator Iterator;
        };
    }; // ... struct Codim ...

  enum
  {
    conforming = RealGridViewTraits::conforming
  };
  enum
  {
    dimension = Grid::dimension
  };
  enum
  {
    dimensionworld = Grid::dimensionworld
  };

  typedef typename Grid::ctype ctype;

  // ...except for the Intersection and IntersectionIterator
  typedef PeriodicIntersection<RealGridViewType> Intersection;
  typedef PeriodicIntersectionIterator<RealGridViewType> IntersectionIterator;
}; // ... class PeriodicGridViewTraits ...

/** \brief Actual Implementation of PeriodicGridView
 *  \see PeriodicGridView
*/
template <class RealGridViewImp, bool use_less_memory >
class PeriodicGridViewImp : public RealGridViewImp
{
  typedef RealGridViewImp BaseType;
  typedef PeriodicGridViewImp< BaseType, use_less_memory > ThisType;
  typedef PeriodicGridViewTraits<BaseType, use_less_memory> Traits;

public:
  typedef typename BaseType::Grid Grid;
  typedef PeriodicIndexSet< BaseType, use_less_memory > IndexSet;
  typedef typename BaseType::template Codim<0>::Entity EntityType;
  typedef PeriodicIntersectionIterator<BaseType> IntersectionIterator;
  typedef typename IntersectionIterator::RealIntersectionType RealIntersectionType;
  typedef typename Traits::IndexSet::IndexType IndexType;
  typedef int IntersectionIndexType;
  typedef typename RealIntersectionType::GlobalCoordinate DomainType;
  typedef PeriodicIntersection<BaseType> Intersection;
  typedef std::map<IntersectionIndexType, std::pair<bool, EntityType>> IntersectionMapType;
  typedef typename std::map<IndexType, IndexType> IndexMapType;
  static const size_t dimDomain = BaseType::dimension;

  template <int cd>
  struct Codim : public Traits::template Codim<cd>
  {
  };

private:
  // function that is used in the constructor to create the index map and the counter for the PeriodicIndexSet
  template< int codim >
  static void create_index_map(const DomainType& lower_left,
                               const DomainType& upper_right,
                               const std::bitset<dimDomain>& periodic_directions,
                               std::vector< IndexMapType >& index_maps,
                               ThisType* this_ptr,
                               const BaseType& real_grid_view,
                               std::vector< size_t >& entity_counts,
                               std::map< Dune::GeometryType, size_t >& type_counts,
                               std::vector< std::set< IndexType >>& entities_to_skip_vector,
                               std::vector< std::vector< IndexType > >& new_indices_vector)
  {
      EntityInlevelSearch<BaseType,codim> entity_search_codim(*this_ptr);
      IndexMapType index_map;
      std::map< IndexType, size_t > entity_to_vector_index_map;
      std::vector< DomainType > periodic_coords_vector;
      const auto& real_index_set = real_grid_view.indexSet();
      if (!use_less_memory)
          new_indices_vector[codim].resize(real_index_set.size(codim));
      auto& new_indices = new_indices_vector[codim];
      size_t num_codim_entities = 0;
      std::set< IndexType > free_indices;
      std::set< IndexType > entities_to_skip;
      for (const auto& entity : entities(*this_ptr, Dune::Codim< codim >())) {
          // increase codim counter
          ++num_codim_entities;
          // increase GeometryType counter
          const auto geometry_type = entity.type();
          if (type_counts.count(geometry_type))
              ++(type_counts.at(geometry_type));
          else
              type_counts.insert(std::make_pair(geometry_type, size_t(1)));
          // check if entity is on a periodic boundary and, if it is, replace index by indices of periodic entity
          auto periodic_coords = entity.geometry().center();
          std::size_t upper_right_coords = 0;
          for (std::size_t ii = 0; ii < dimDomain; ++ii) {
              if (periodic_directions[ii]) {
                  if (XT::Common::FloatCmp::eq(periodic_coords[ii], upper_right[ii])) {
                      ++upper_right_coords;
                      periodic_coords[ii] = lower_left[ii];
                  }
              }
          }
          if (upper_right_coords > 0) { // find periodic adjacent entity
              // dont count this entity as it is identified with the single periodic equivalent entity
              --num_codim_entities;
              --(type_counts.at(geometry_type));
              free_indices.insert(real_index_set.index(entity));
              entities_to_skip.insert(real_index_set.index(entity));
              periodic_coords_vector.push_back(periodic_coords);
              entity_to_vector_index_map.insert(std::make_pair(real_index_set.index(entity), periodic_coords_vector.size() - 1));
          }
      } // walk entities in a given codimension

      entity_counts[codim] = num_codim_entities;
      entities_to_skip_vector[codim] = entities_to_skip;

      // find periodic entities
      const auto periodic_entity_ptrs = entity_search_codim(periodic_coords_vector);

      // assign new indices to the entities
      IndexType current_new_index = 0;

      // to save memory, only entities with an index greater than the number of entities in this codim will get a new index
      if (use_less_memory) {
          for (const auto& index : free_indices)
              if (index >= num_codim_entities)
                  free_indices.erase(index);
      }

      // assign new indices to entities that are not replaced by their periodic equivalent entity
      for (const auto& entity : entities(*this_ptr, Dune::Codim< codim >()))
      {
          const auto old_index = real_index_set.index(entity);
          if (!entity_to_vector_index_map.count(old_index)) {
              if (!use_less_memory) {
                  new_indices[old_index] = current_new_index;
                  ++current_new_index;
              } else {
                  if (old_index >= num_codim_entities) {
                      assert(free_indices.begin() != free_indices.end());
                      auto free_index = *free_indices.begin();
                      index_map.insert(std::make_pair(old_index, free_index));
                      free_indices.erase(free_index);
                  }
              }
          }
      }

      // assign index of periodic equivalent entity to entities that are replaced
      for (const auto& pair : entity_to_vector_index_map) {
          const auto& periodic_entity_ptr = periodic_entity_ptrs.at(pair.second);
          if (periodic_entity_ptr == nullptr)
              DUNE_THROW(Dune::InvalidStateException, "Could not find periodic neighbor entity");
          auto index_of_periodic_entity = real_index_set.index(*periodic_entity_ptr);
          if (!use_less_memory)
              new_indices[pair.first] = new_indices[index_of_periodic_entity];
          else
              index_map.insert(std::make_pair(pair.first, index_map.count(index_of_periodic_entity) ? index_map.at(index_of_periodic_entity) : index_of_periodic_entity));
      }
      index_maps[codim] = index_map;
  }

  // compile time for loop to loop over the codimensions in constructor, see http://stackoverflow.com/a/11081785
  template< int codim, int to >
  struct static_for_loop_for_index_maps
  {
      void operator()(const DomainType& lower_left,
                      const DomainType& upper_right,
                      const std::bitset<dimDomain>& periodic_directions,
                      std::vector< IndexMapType >& index_maps,
                      ThisType* this_ptr,
                      const BaseType& real_grid_view,
                      std::vector< size_t >& entity_counts,
                      std::map< Dune::GeometryType, size_t >& type_counts,
                      std::vector< std::set< IndexType >>& entities_to_skip_vector,
                      std::vector< std::vector< IndexType > >& new_indices_vector)
      {
          create_index_map< codim >(lower_left, upper_right, periodic_directions, index_maps, this_ptr, real_grid_view, entity_counts, type_counts, entities_to_skip_vector, new_indices_vector);
          static_for_loop_for_index_maps< codim + 1, to >()(lower_left, upper_right, periodic_directions, index_maps, this_ptr, real_grid_view, entity_counts, type_counts, entities_to_skip_vector, new_indices_vector);
      }
  };

  // specialization of static for loop to end the loop
  template< int to >
  struct static_for_loop_for_index_maps< to, to>
  {
      void operator()(const DomainType& /*lower_left*/,
                      const DomainType& /*upper_right*/,
                      const std::bitset<dimDomain>& /*periodic_directions*/,
                      std::vector< IndexMapType >& /*index_maps*/,
                      ThisType* /*this_ptr*/,
                      const BaseType& /*real_grid_view*/,
                      std::vector< size_t >& /*entity_counts*/,
                      std::map< Dune::GeometryType, size_t >& /*type_counts*/,
                      std::vector< std::set< IndexType >>& /*entities_to_skip_vector*/,
                      std::vector< std::vector< IndexType > >& /*new_indices_vector*/)
      {
      }
  };

public:
  PeriodicGridViewImp(const BaseType& real_grid_view, const std::bitset<dimDomain> periodic_directions)
      : BaseType(real_grid_view)
      , index_maps_(dimDomain+1)
      , empty_intersection_map_(IntersectionMapType())
      , periodic_directions_(periodic_directions)
      , entities_to_skip_vector_(dimDomain+1)
      , new_indices_vector_(dimDomain+1)
      , real_index_set_(BaseType::indexSet())
  {
      // find lower left and upper right corner of the grid
      auto entity_it         = BaseType::template begin<0>();
      DomainType lower_left  = entity_it->geometry().center();
      DomainType upper_right = lower_left;
      for (const auto& entity : Dune::elements(*this)) {
          if (entity.hasBoundaryIntersections()) {
              const auto i_it_end = BaseType::iend(entity);
              for (auto i_it = BaseType::ibegin(entity); i_it != i_it_end; ++i_it) {
                  const RealIntersectionType& intersection = *i_it;
                  const auto intersection_coords = intersection.geometry().center();
                  for (std::size_t ii = 0; ii < dimDomain; ++ii) {
                      if (intersection_coords[ii] > upper_right[ii])
                          upper_right[ii] = intersection_coords[ii];
                      if (intersection_coords[ii] < lower_left[ii])
                          lower_left[ii] = intersection_coords[ii];
                  }
              }
          }
      }

      /* walk the grid and create a map that maps each entity(index) on the boundary to a map that maps every
      local intersection index to a std::pair< bool, EntityType >, where the bool component of the pair
      indicates whether the intersection is on a periodic boundary and the second component is the periodic neighbor
      entity if the first component is true. If the first component is false, the Entity is not meant to be used. */
      EntityInlevelSearch<BaseType> entity_search(*this);
      DomainType periodic_neighbor_coords;
      IntersectionMapType intersection_neighbor_map;
      std::map< Dune::GeometryType, size_t > type_counts;
      std::vector< DomainType > periodic_neighbor_coords_vector;
      /* we don't want to do the entitysearch on each entity separately, so we collect the coordinates in a vector to do
       * the search after the gridwalk. This map maps the entity indices to a map mapping local intersection indices to
       * the index of the corresponding coordinate in the periodic_neighbor_coords_vector. */
      std::map< IndexType, std::map< IntersectionIndexType, size_t > > entity_to_intersection_to_vector_index_map;
      std::map< IntersectionIndexType, size_t > intersection_to_vector_index_map;
      for (const auto& entity : Dune::elements(*this)) {
          // count entities per geometry type for the PeriodicIndexSet
          const auto geometry_type = entity.type();
          if (type_counts.count(geometry_type))
              ++(type_counts.at(geometry_type));
          else
              type_counts.insert(std::make_pair(geometry_type, size_t(1)));
          if (entity.hasBoundaryIntersections()) {
              intersection_neighbor_map.clear();
              intersection_to_vector_index_map.clear();
              const auto i_it_end = BaseType::iend(entity);
              for (auto i_it = BaseType::ibegin(entity); i_it != i_it_end; ++i_it) {
                  const RealIntersectionType& intersection    = *i_it;
                  const IntersectionIndexType index_in_inside = intersection.indexInInside();
                  bool is_periodic = false;
                  if (intersection.boundary()) {
                      periodic_neighbor_coords   = intersection.geometry().center();
                      size_t num_boundary_coords = 0;
                      for (std::size_t ii = 0; ii < dimDomain; ++ii) {
                          if (periodic_directions_[ii]) {
                              if (XT::Common::FloatCmp::eq(periodic_neighbor_coords[ii], lower_left[ii])) {
                                  is_periodic = true;
                                  periodic_neighbor_coords[ii] =
                                          upper_right[ii] - 1.0 / 100.0 * (entity.geometry().center()[ii] - lower_left[ii]);
                                  ++num_boundary_coords;
                              } else if (XT::Common::FloatCmp::eq(periodic_neighbor_coords[ii], upper_right[ii])) {
                                  is_periodic = true;
                                  periodic_neighbor_coords[ii] =
                                          lower_left[ii] + 1.0 / 100.0 * (upper_right[ii] - entity.geometry().center()[ii]);
                                  ++num_boundary_coords;
                              }
                          }
                      }
                      assert(num_boundary_coords = 1);
                      if (is_periodic) {
                          periodic_neighbor_coords_vector.push_back(periodic_neighbor_coords);
                          intersection_to_vector_index_map.insert(std::make_pair(index_in_inside, periodic_neighbor_coords_vector.size() - 1));
                      } else {
                          intersection_neighbor_map.insert(
                                      std::make_pair(index_in_inside, std::make_pair(is_periodic, EntityType(entity))));
                      }
                  } else {
                      intersection_neighbor_map.insert(
                                  std::make_pair(index_in_inside, std::make_pair(bool(false), EntityType(entity))));
                  }
              }
              entity_to_intersection_to_vector_index_map.insert(std::make_pair(real_index_set_.index(entity), intersection_to_vector_index_map));
              entity_to_intersection_map_map_.insert(
                          std::make_pair(real_index_set_.index(entity), intersection_neighbor_map));
          }
      } // walk Entities

      // search the periodic neighbor entities
      const auto outside_ptrs = entity_search(periodic_neighbor_coords_vector);

      for (auto& pair : entity_to_intersection_to_vector_index_map) {
          for (auto& intersection_to_vector : pair.second) {
              const auto& outside_ptr = outside_ptrs.at(intersection_to_vector.second);
              if (outside_ptr == nullptr)
                  DUNE_THROW(Dune::InvalidStateException, "Could not find periodic neighbor entity");
              entity_to_intersection_map_map_.at(pair.first).insert(
                          std::make_pair(intersection_to_vector.first, std::make_pair(bool(true), *outside_ptr)));
          }
      }

      /* walk the grid for each codimension from 1 to dimDomain and create a map mapping indices from entitys of that
       * codimension on a periodic boundary to the index of the corresponding periodic equivalent entity that has the
       * most coordinates in common with the lower left corner of the grid */
      std::vector< size_t > entity_counts(dimDomain+1);
      entity_counts[0] = real_grid_view.indexSet().size(0);
      static_for_loop_for_index_maps<1, dimDomain + 1>()(lower_left, upper_right, periodic_directions_, index_maps_, this, real_grid_view, entity_counts, type_counts, entities_to_skip_vector_, new_indices_vector_);
      // create index_set
      index_set_ = std::make_shared< IndexSet >(real_index_set_, index_maps_, entity_counts, type_counts, new_indices_vector_);
  } // constructor PeriodicGridViewImp(...)

  int size (int codim) const
  {
      return index_set_->size(codim);
  }

  int size (const Dune::GeometryType &type) const
  {
      return index_set_->size(type);
  }

  template<int cd>
  typename Codim< cd >::Iterator begin() const
  {
      return typename Codim< cd >::Iterator(BaseType::template begin<cd>(), &(entities_to_skip_vector_[cd]), &real_index_set_, BaseType::template end<cd>());
  }

  template<int cd>
  typename Codim< cd >::Iterator end() const
  {
      return typename Codim< cd >::Iterator(BaseType::template end<cd>(), &(entities_to_skip_vector_[cd]), &real_index_set_, BaseType::template end<cd>());
  }

  template<int cd, PartitionIteratorType pitype>
  typename Codim< cd >::template Partition< pitype >::Iterator 	begin () const
  {
      return typename Codim<cd>::template Partition<pitype>::Iterator(BaseType::template begin<cd,pitype>(), &(entities_to_skip_vector_[cd]), &real_index_set_, BaseType::template end<cd,pitype>());
  }


  template<int cd, PartitionIteratorType pitype>
  typename Codim< cd >::template Partition< pitype >::Iterator end () const
  {
      return typename Codim<cd>::template Partition<pitype>::Iterator(BaseType::template begin<cd,pitype>(), &(entities_to_skip_vector_[cd]), &real_index_set_, BaseType::template end<cd,pitype>());
  }

  const IndexSet& indexSet() const
  {
      return *index_set_;
  }

  IntersectionIterator ibegin(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(BaseType::ibegin(entity),
                                *this,
                                entity,
                                entity.hasBoundaryIntersections()
                                    ? entity_to_intersection_map_map_.at(this->indexSet().index(entity))
                                    : (const IntersectionMapType&)empty_intersection_map_);

  } // ... ibegin(...)

  IntersectionIterator iend(const typename Codim<0>::Entity& entity) const
  {
    return IntersectionIterator(BaseType::iend(entity),
                                *this,
                                entity,
                                entity.hasBoundaryIntersections()
                                    ? entity_to_intersection_map_map_.at(this->indexSet().index(entity))
                                    : (const IntersectionMapType&)empty_intersection_map_);
  } // ... iend(...)

private:
  std::map<IndexType, IntersectionMapType> entity_to_intersection_map_map_;
  std::vector< IndexMapType > index_maps_;
  const IntersectionMapType empty_intersection_map_;
  const std::bitset<dimDomain> periodic_directions_;
  std::shared_ptr< IndexSet > index_set_;
  std::vector<std::set< IndexType >> entities_to_skip_vector_;
  std::vector< std::vector< IndexType > > new_indices_vector_;
  const typename BaseType::IndexSet& real_index_set_;
}; // ... class PeriodicGridViewImp ...

} // namespace internal

/** \brief GridView that takes an arbitrary Dune::GridView and adds periodic boundaries
 *
 * PeriodicGridView is templated by and derived from an arbitrary Dune::GridView. All methods are forwarded to the
 * underlying GridView except for the ibegin and iend methods. These methods return a PeriodicIntersectionIterator
 * which again behaves like the underlying IntersectionIterator except that it returns a PeriodicIntersection in its
 * operator*. The PeriodicIntersection again behaves like an Intersection of the underlying GridView, but may return
 * neighbor() == true and an outside() entity even if it is on the boundary. The outside() entity is the entity
 * adjacent to the intersection if it is identified with the intersection on the other side of the grid.
 * In the constructor, PeriodicGridViewImp will build a map mapping boundary entity indices to a map mapping local
 * intersection indices to a std::pair containing the information whether this intersection shall be periodic and the
 * outside entity. This may take quite long as finding the outside entity requires a grid walk for each periodic
 * intersection.
 * By default, all coordinate directions will be made periodic. By supplying a std::bitset< dimension > you can decide
 * for each direction whether it should be periodic (1 means periodic, 0 means 'behave like underlying GridView in that
 * direction').

   \note
      -  Currently, PeriodicGridView will only work with GridViews on the unit hypercube
      -  Only cube and regular simplex grids have been tested so far. Other grids may not work properly. This is due to
      the heuristics for finding the periodic neighbor entity: Given an intersection on the boundary that shall be
      periodic, the coordinates intersection.geometry().center() are moved to the other side of the grid and then
      supplied to Dune::XT::Grid::EntityInLevelSearch. As the coordinates are on the boundary of the wanted entity,
      this search will fail for some grids. Thus, the coordinates are moved a little to the inside of the grid before
      searching for the entity. The moved coordinates will be inside the wanted entity for cube and usual simplex grids
      but this is not guaranteed for arbitrary grids.
 */
template <class RealGridViewImp, bool use_less_memory = false>
class PeriodicGridView : XT::Common::ConstStorageProvider<internal::PeriodicGridViewImp<RealGridViewImp, use_less_memory>>,
                         public Dune::GridView<internal::PeriodicGridViewTraits<RealGridViewImp, use_less_memory>>
{
  typedef RealGridViewImp RealGridViewType;
  typedef typename Dune::GridView<internal::PeriodicGridViewTraits<RealGridViewType, use_less_memory>> BaseType;
  typedef
      typename XT::Common::ConstStorageProvider<internal::PeriodicGridViewImp<RealGridViewImp, use_less_memory>> ConstStorProv;
  typedef typename RealGridViewType::template Codim<0>::Geometry::GlobalCoordinate DomainType;

public:
  static const size_t dimension = RealGridViewType::dimension;

  PeriodicGridView(const RealGridViewType& real_grid_view,
                   const std::bitset<dimension> periodic_directions = std::bitset<dimension>().set())
    : ConstStorProv(new internal::PeriodicGridViewImp<RealGridViewType, use_less_memory>(real_grid_view, periodic_directions))
    , BaseType(ConstStorProv::access())
  {
  }

  PeriodicGridView(const PeriodicGridView& other)
    : ConstStorProv(new internal::PeriodicGridViewImp<RealGridViewType, use_less_memory>(other.access()))
    , BaseType(ConstStorProv::access())
  {
  }
}; // class PeriodicGridView

#else // HAVE_DUNE_GRID

template <class RealGridViewImp>
class PeriodicGridView
{
  static_assert(AlwaysFalse<RealGridViewImp>::value, "You are missing dune-grid!");
};

#endif // HAVE_DUNE_GRID


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PERIODICVIEW_HH
