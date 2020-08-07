// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2017, 2020)

#ifndef DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_INTERNAL_HH
#define DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_INTERNAL_HH

// This file is copied from the EXADUNE dune-grid module,
// https://gitlab.dune-project.org/exadune/dune-grid, commit 96d51e79acad503bf430e7c3cfde6cbd0515f621,
// with adaptions to remove the dune/grid/utility/entityrange.hh include and to be usable with grid views and grid
// parts alike.

// clang-format off
#include <dune/xt/common/disable_warnings.hh>

#include <iterator>
#include <vector>

#if HAVE_TBB
#include <tbb/tbb_stddef.h>
#endif

#include <dune/xt/grid/type_traits.hh>

namespace Dune::XT::Grid {
    //! Partioning base on remembering iterator ranges
    template<class GridView, int codim, PartitionIteratorType pit = All_Partition>
    class RangedPartitioning {
      typedef XT::Grid::extract_iterator_t<GridView, codim, pit> Iterator;
    public:
      //! type of partitions
      class Partition;

      //! type used to count partitions
      using Size = typename std::iterator_traits<typename std::vector<Iterator>::iterator>::difference_type;

      //! construct
      RangedPartitioning(const GridView &gv, Size partitions) {
        const Iterator end = gv.template end<0, pit>();
        // GridView's size() would return the wrong count for non AllPartition Views/Iterators
        const auto totalsize = std::distance(gv.template begin<0, pit>(), end);
        Iterator it = gv.template begin<0, pit>();
        entry_points_.reserve(partitions + 1);
        stride_ = totalsize / partitions;
        overflow_ = totalsize % partitions;
        for (Size p = 0; p < partitions; ++p) {
          entry_points_.push_back(it);
          std::advance(it, stride_ + (p < overflow_));
        }
        entry_points_.push_back(end);
      }

      //! return maximum number of partitions
      Size partitions() const {
        return entry_points_.size() - 1;
      }

      //! whole partitioning as a partition object
      /**
       * The partition can be split by tbb
       */
      Partition everything() const {
        return Partition(*this, 0, partitions());
      }

      //! return a particular partition
      Partition partition(Size pId) const {
        return Partition(*this, pId, pId + 1);
      }
      //! return a range of partitions
      /**
       * The returned partition object can be split by tbb
       */
      Partition partitions(Size first, Size last) const {
        return Partition(*this, first, last);
      }

    private:
      std::vector<Iterator> entry_points_;
      Size stride_;
      Size overflow_;
    };

    //! Partition type for ranged partitioning
    /**
     * \implements EntityRangeInterface
     *
     * This implements the interface for entity ranges.  It differs from
     * IteratorEntityRange in that it is TBB-splittable and thus needs to hold a
     * reference to the partitioning.
     */
    template<class GridView, int codim, PartitionIteratorType pit>
    class RangedPartitioning<GridView, codim, pit>::Partition {
    public:
      //! type of iterator
      typedef XT::Grid::extract_iterator_t<GridView, codim, pit> Iterator;
      //! type of entity
      using Entity = XT::Grid::extract_entity_t<GridView>;
      //! type used to count entites
      using Size = typename std::iterator_traits<Iterator>::difference_type;

      //! Construct partition from scratch
      Partition(const RangedPartitioning &partitioning,
                RangedPartitioning::Size firstPartition,
                RangedPartitioning::Size lastPartition) :
          partitioning_(partitioning), firstPartition_(firstPartition),
          lastPartition_(lastPartition) {}

      //! Create a begin iterator
      const Iterator &begin() const {
        return partitioning_.entry_points_[firstPartition_];
      }

      //! Create an end iterator
      const Iterator &end() const {
        return partitioning_.entry_points_[lastPartition_];
      }

      //! Number of Elements visited by an iterator
      Size size() const {
        return (lastPartition_ - firstPartition_) * partitioning_.stride_
               + (firstPartition_ < partitioning_.overflow_
                  ? partitioning_.overflow_ - firstPartition_
                  : 0);
      }

#if HAVE_TBB
      //! Splitting Constructor
      /**
       * Construct second half of set, update \c other to represent first half.
       *
       * This is done by looking up the entry point in the partitioning that
       * lies roughly in the middle of the current \c begin() and \c end().
       */
      Partition(Partition &other, tbb::split) :
        partitioning_(other.partitioning_)
      {
        firstPartition_ = (other.firstPartition_ + other.lastPartition_ + 1) / 2;
        lastPartition_ = other.lastPartition_;
        other.lastPartition_ = firstPartition_;
      }
      //! check whether set is empty (for TBB)
      /**
       * This is equivalent to \c size()==0.
       */
      bool empty() const
      {
        return size() == 0;
      }
      //! check whether set can be split (for TBB)
      bool is_divisible() const
      {
        return lastPartition_ - firstPartition_ > 1;
      }
#endif // HAVE_TBB

    private:
      const RangedPartitioning &partitioning_;
      RangedPartitioning::Size firstPartition_;
      RangedPartitioning::Size lastPartition_;
    };

  } // namespace Dune

#include <dune/xt/common/reenable_warnings.hh>
// clang-format on

#endif // DUNE_XT_GRID_PARALLEL_PARTITIONING_RANGED_INTERNAL_HH
