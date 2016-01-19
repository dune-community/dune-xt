// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2015)
//   Rene Milk       (2012 - 2013, 2015)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_INTERSECTION_HH
#define DUNE_XT_GRID_INTERSECTION_HH

#include <dune/common/deprecated.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#if HAVE_DUNE_GRID
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/gridview.hh>
#endif

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_utils.hh>

namespace Dune {
namespace XT {
namespace Grid {

#if HAVE_DUNE_GRID

template <class GridPartOrViewType>
class Intersection
{
  template <class GridViewType, bool is_view>
  struct Choose
  {
    typedef typename GridViewType::Intersection Type;
  };

  template <class GridPartType>
  struct Choose<GridPartType, false>
  {
    typedef typename GridPartType::IntersectionType Type;
  };

  static const bool this_is_a_grid_view =
      std::is_base_of<GridView<typename GridPartOrViewType::Traits>, GridPartOrViewType>::value;

public:
  typedef typename Choose<GridPartOrViewType, this_is_a_grid_view>::Type Type;
}; // class Intersection

#endif // HAVE_DUNE_GRID

/**
  \brief      prints some basic information about a Dune::Intersection, namely the number of its corners and the
              coordinates of those corners.
  \tparam     IntersectionType
              Dune::Intersection compatible
  \param[in]  intersection
              Dune::Intersection, whose information should be printed
  \param[out] stream
              std::ostream, into which the information is printed
  **/
template <class IntersectionType>
void printIntersection(const IntersectionType& intersection, std::ostream& out = std::cout,
                       const std::string prefix = "")
{
  out << prefix << Common::Typename<IntersectionType>::value() << std::endl;
  const auto& geometry = intersection.geometry();
  for (auto ii : Common::value_range(geometry.corners()))
    out << prefix << "  corner " + Common::to_string(ii) << " = " << geometry.corner(ii)
        << " (local: " << geometry.local(geometry.corner(ii)) << ")\n";
} // ... printIntersection(...)

#if HAVE_DUNE_GRID
/**
 * \brief Checks if intersection contains the given global_point.
 *
 *        Returns true, if global_point lies on the line between the corners of intersection.
 */
template <class G, class I, class D>
typename std::enable_if<Dune::Intersection<G, I>::dimension == 2, bool>::type
contains(const Dune::Intersection<G, I>& intersection, const Dune::FieldVector<D, 2>& global_point,
         const D& tolerance = Common::FloatCmp::DefaultEpsilon<D>::value())
{
  const auto& geometry = intersection.geometry();
  // get the global coordinates of the intersections corners
  assert(geometry.corners() == 2);
  const auto corner_0 = geometry.corner(0);
  const auto corner_1 = geometry.corner(1);
  // A line is given by $y = a*x + b$. Computing a and b fails for certain intersections (for instance those
  // parallel to the y axis). So in order to check if the point is on the line between the corners we consider the
  // vectors pointing from the point to each corner. If those are not orthogonal to the intersections normal, the point
  // cannot lie on the line between the two corners.
  const auto normal = intersection.centerUnitOuterNormal();
  for (auto vector : {corner_0 - global_point, corner_1 - global_point})
    if (vector * normal > tolerance)
      return false;
  // Now that we know the point is on the line, check if it is outside the bounding box of the two corners.
  if (Common::FloatCmp::lt(global_point[0], std::min(corner_0[0], corner_1[0]), tolerance)
      || Common::FloatCmp::gt(global_point[0], std::max(corner_0[0], corner_1[0]), tolerance))
    return false;
  if (Common::FloatCmp::lt(global_point[1], std::min(corner_0[1], corner_1[1]), tolerance)
      || Common::FloatCmp::gt(global_point[1], std::max(corner_0[1], corner_1[1]), tolerance))
    return false;
  // At this point we cannot reject the assumption that the point lies on the line between the two corners.
  return true;
} // ... contains(...)
#endif

/** Check whether a spatial point lies on an intersection.
*
* @param[in] intersection The intersection
* @param[in] globalPoint A Dune::FieldVector with the global coordinates of the point
* @return Returns true if the point lies on the intersection, false otherwise.
*/
template <class IntersectionType, class FieldType, int dim>
bool DUNE_DEPRECATED_MSG("This method does not produce correct results, use contains() instead!")
    intersectionContains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, dim>& globalPoint)
{
  // map global coordinates to local coordinates of the intersection
  const auto& intersectionGeometry = intersection.geometry();
  const auto& localPoint           = intersectionGeometry.local(globalPoint);

// get codim 1 reference element
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY, 2, 3)
  const auto& refElement = ReferenceElements<FieldType, dim - 1>::general(intersectionGeometry.type());
#else
  const auto& refElement = GenericReferenceElements<FieldType, dim - 1>::general(intersectionGeometry.type());
#endif
  // check whether reference element contains the local coordinates
  return refElement.checkInside(localPoint);
} // end function intersectionContains

template <class IntersectionType, class FieldType>
bool intersectionContains(const IntersectionType& intersection, const Dune::FieldVector<FieldType, 2>& globalPoint)
{
  return contains(intersection, globalPoint);
}

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_INTERSECTION_HH
