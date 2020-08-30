// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2019)
//   Tobias Leibner  (2019 - 2020)

#ifndef DUNE_XT_GRID_INTEGRALS_HH
#define DUNE_XT_GRID_INTEGRALS_HH

#include <functional>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/entity.hh>

#include <dune/xt/grid/type_traits.hh>

namespace Dune {
namespace XT {
namespace Grid {


/**
 * \brief Computes the integral of a given function over a given element [most general variant].
 */
template <class RangeType, class Element>
RangeType element_integral(
    const Element& element,
    std::function<RangeType(
        const FieldVector<typename Element::Geometry::ctype, Element::dimension>& point_in_reference_element)> function,
    const int polynomial_order_of_the_function)
{
  static_assert(XT::Grid::is_entity<Element>::value, "element has to be a codim-0 grid entity");
  static_assert(Element::codimension == 0, "element has to be a codim-0 grid entity");
  RangeType result(0.), local_result;
  for (auto&& quadrature_point : QuadratureRules<typename Element::Geometry::ctype, Element::dimension>::rule(
           element.type(), polynomial_order_of_the_function)) {
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
template <class Element>
double element_integral(const Element& element,
                        std::function<double(const FieldVector<typename Element::Geometry::ctype, Element::dimension>&
                                                 point_in_reference_element)> function,
                        const int polynomial_order_of_the_function)
{
  return element_integral<double, Element>(element, function, polynomial_order_of_the_function);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_INTEGRALS_HH
