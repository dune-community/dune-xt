// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Andreas Buhr    (2014)
//   Felix Schindler (2012 - 2019)
//   René Fritze     (2012 - 2013, 2015 - 2016, 2018 - 2019)
//   Sven Kaulmann   (2013)
//   Tobias Leibner  (2014, 2020)

#ifndef DUNE_XT_GRID_INTERSECTION_HH
#define DUNE_XT_GRID_INTERSECTION_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/gridview.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/print.hh>
#include <dune/xt/common/ranges.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/xt/grid/type_traits.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

namespace Dune {
namespace XT::Grid {


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
template <class G, class I>
void print_intersection(const Intersection<G, I>& intersection,
                        const std::string name = Common::Typename<Intersection<G, I>>::value(),
                        std::ostream& out = std::cout,
                        const std::string prefix = "")
{
  if (!name.empty())
    out << prefix << name << ":\n";
  const auto& geometry = intersection.geometry();
  for (auto ii : Common::value_range(geometry.corners()))
    out << prefix << "  corner " + Common::to_string(ii) << " = " << geometry.corner(ii)
        << " (local: " << geometry.local(geometry.corner(ii)) << ")\n";
} // ... print_intersection(...)


template <class G, class I>
double diameter(const Intersection<G, I>& intersection)
{
  auto max_dist = std::numeric_limits<typename G::ctype>::min();
  const auto& geometry = intersection.geometry();
  for (auto i : Common::value_range(geometry.corners())) {
    const auto xi = geometry.corner(i);
    for (auto j : Common::value_range(i + 1, geometry.corners())) {
      auto xj = geometry.corner(j);
      xj -= xi;
      max_dist = std::max(max_dist, xj.two_norm());
    }
  }
  return max_dist;
} // diameter


/**
 * \brief Checks if intersection contains the given global_point (1d variant).
 *
 *        Returns true, if global_point and intersection coincide.
 */
template <class G, class I, class D>
bool contains(const Dune::Intersection<G, I>& intersection,
              const Dune::FieldVector<D, 1>& global_point,
              const D& tolerance = Common::FloatCmp::DefaultEpsilon<D>::value())
{
  return Common::FloatCmp::eq(intersection.geometry().center(), global_point, tolerance);
}


/**
 * \brief Checks if intersection contains the given global_point (2d variant).
 *
 *        Returns true if global_point lies on the line between the corners of intersection.
 */
template <class G, class I, class D>
bool contains(const Dune::Intersection<G, I>& intersection,
              const Dune::FieldVector<D, 2>& global_point,
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


/**
 * \brief Checks if intersection contains the given global_point (3d variant).
 *
 *        Returns true if global_point lies within the plane spanned by the first three corners of intersection.
 *
 * \sa
 * http://math.stackexchange.com/questions/684141/check-if-a-point-is-on-a-plane-minimize-the-use-of-multiplications-and-divisio
 */
template <class G, class I, class D>
bool contains(const Dune::Intersection<G, I>& intersection,
              const Dune::FieldVector<D, 3>& global_point,
              const D& tolerance = Common::FloatCmp::DefaultEpsilon<D>::value())
{
  const auto& geometry = intersection.geometry();
  // get the global coordinates of the intersections corners, there should be at least 3 (ignore the fourth if there is
  // one, 3 points is enough in 3d)
  assert(geometry.corners() >= 3);
  std::vector<Dune::FieldVector<D, 3>> points(4, Dune::FieldVector<D, 3>(0.));
  for (size_t ii = 0; ii < 3; ++ii)
    points[ii] = geometry.corner(boost::numeric_cast<int>(ii));
  points[3] = global_point;
  // form a matrix of these points, where the top three entries of each column are given by the three entries of each
  // point and the bottom row is given by one, i.e.
  // a_0 b_0 c_0 d_0
  // a_1 b_1 c_1 d_1
  // a_2 b_2 c_2 d_2
  //   1   1   1   1
  FieldMatrix<D, 4, 4> matrix(1.); // ensures the 1 on the last row
  for (size_t ii = 0; ii < 3; ++ii) // only set the first three rows
    for (size_t jj = 0; jj < 4; ++jj)
      matrix[ii][jj] = points[jj][ii];
  // the point lies on the plane given by the corners if the determinant of this matrix is zero
  const D det = matrix.determinant();
  return std::abs(det) < tolerance;
} // ... contains(...)


} // namespace XT::Grid


template <class G, class I>
std::ostream& DUNE_DEPRECATED_MSG("Use out << print(intersection) from <dune/xt/grid/print.hh> instead (05.07.2020)!")
operator<<(std::ostream& out, const Dune::Intersection<G, I>& intersection)
{
  XT::Common::print(intersection);
  return out;
}


} // namespace Dune

#endif // DUNE_XT_GRID_INTERSECTION_HH
