// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)
//   René Fritze     (2020)

#ifndef DUNE_XT_GRID_PRINT_HH
#define DUNE_XT_GRID_PRINT_HH

#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/entity.hh>

#include <dune/xt/common/print.hh>

namespace Dune::XT::Common {


/// \sa Common::Printer
template <class G, class I, bool use_repr>
class Printer<Dune::Intersection<G, I>, use_repr> : public internal::DefaultPrinter<Dune::Intersection<G, I>, use_repr>
{
public:
  using T = Dune::Intersection<G, I>;

  const std::string grid_name;

  Printer(const T& val, const Configuration& cfg = {})
    : internal::DefaultPrinter<T, use_repr>(val, cfg)
  {}

  void repr(std::ostream& out) const final
  {
    const auto& geometry = this->value.geometry();
    const auto num_corners = geometry.corners();
    out << G::dimension - 1 << "d-intersection{corners=" << print(geometry.corner(0));
    for (int ii = 1; ii < num_corners; ++ii)
      out << ", " << print(geometry.corner(ii));
    out << ", normal=" << this->value.centerUnitOuterNormal() << "}";
  }

  void str(std::ostream& out) const final
  {
    const auto& geometry = this->value.geometry();
    const auto num_corners = geometry.corners();
    out << G::dimension - 1 << "d-intersection{" << print(geometry.corner(0));
    for (int ii = 1; ii < num_corners; ++ii)
      out << ", " << print(geometry.corner(ii));
    out << "}";
  } // ... str(...)
}; // class Printer<Dune::Intersection<G, I>, ...>


/// \sa Common::Printer
template <int cd, int dim, class GridImp, template <int, int, class> class EntityImp, bool use_repr>
class Printer<Dune::Entity<cd, dim, GridImp, EntityImp>, use_repr>
  : public internal::DefaultPrinter<Dune::Entity<cd, dim, GridImp, EntityImp>, use_repr>
{
public:
  using T = Dune::Entity<cd, dim, GridImp, EntityImp>;

  const std::string grid_name;

  Printer(const T& val, const Configuration& cfg = {})
    : internal::DefaultPrinter<T, use_repr>(val, cfg)
  {}

  void repr(std::ostream& out) const final
  {
    const auto& geometry = this->value.geometry();
    const auto num_corners = geometry.corners();
    out << dim - cd << "d";
    if (cd != 0)
      out << "(" << dim << "d-" << cd << "d)";
    out << "-entity{corners=" << print(geometry.corner(0));
    for (int ii = 1; ii < num_corners; ++ii)
      out << ", " << print(geometry.corner(ii));
    out << "}";
  }

  void str(std::ostream& out) const final
  {
    const auto& geometry = this->value.geometry();
    const auto num_corners = geometry.corners();
    out << dim - cd << "d";
    if (cd != 0)
      out << "(" << dim << "d-" << cd << "d)";
    out << "-entity{corners=" << print(geometry.corner(0));
    for (int ii = 1; ii < num_corners; ++ii)
      out << ", " << print(geometry.corner(ii));
    out << "}";
  } // ... str(...)
}; // class Printer<Entity<cd, dim, GridImp, EntityImp>, ...>


} // namespace Dune::XT::Common

#endif // DUNE_XT_GRID_PRINT_HH
