// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2018 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   Tobias Leibner  (2019)

#ifndef DUNE_XT_GRID_INTEGRALS_HH
#define DUNE_XT_GRID_INTEGRALS_HH

#include <functional>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/entity.hh>

namespace Dune {
namespace XT {
namespace Grid {


/**
 * \brief Computes the integral of a given function over a given element [most general variant].
 */
// We would have liked to use
//   std::function<F(const FieldVector<typename Entity<0, d, G, E>::Geometry::ctype, d>& point_in_reference_element)>
// (note the d), which did not compile for generic lambdas. Thus we use
//   G::dimension
// instead.
template <class RangeType, int d, class G, template <int, int, class> class E>
RangeType
element_integral(const Entity<0, d, G, E>& element,
                 std::function<RangeType(const FieldVector<typename Entity<0, d, G, E>::Geometry::ctype, G::dimension>&
                                             point_in_reference_element)> function,
                 const int polynomial_order_of_the_function)
{
  RangeType result(0.), local_result;
  for (auto&& quadrature_point : QuadratureRules<typename Entity<0, d, G, E>::Geometry::ctype, G::dimension>::rule(
           element.geometry().type(), polynomial_order_of_the_function)) {
    const auto point_in_reference_element = quadrature_point.position();
    const auto quadrature_weight = quadrature_point.weight();
    const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
    local_result = function(point_in_reference_element);
    local_result *= quadrature_weight * integration_factor;
    result += local_result;
  }
  return result;
} // ... element_integral(...)


/**
 * \brief Computes the integral of a given function over a given element [uses double as FieldType].
 */
template <int d, class G, template <int, int, class> class E>
double element_integral(const Entity<0, d, G, E>& element,
                        std::function<double(const FieldVector<typename Entity<0, d, G, E>::Geometry::ctype,
                                                               G::dimension>& point_in_reference_element)> function,
                        const int polynomial_order_of_the_function)
{
  return element_integral<double>(element, function, polynomial_order_of_the_function);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_INTEGRALS_HH
