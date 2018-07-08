// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   Rene Milk       (2016, 2018)
//   Tobias Leibner  (2016 - 2018)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH
#define DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH

#include <dune/xt/grid/exceptions.hh>
#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {

// We do not want to add a virtual destructor (to be able to use this as constexpr), so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#endif


class NoBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "no_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new NoBoundary();
  }
};


class UnknownBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "unknown_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new UnknownBoundary();
  }
};


class DirichletBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "dirichlet_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new DirichletBoundary();
  }
};


class NeumannBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "neumann_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new NeumannBoundary();
  }
};


class RobinBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "robin_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new RobinBoundary();
  }
};


class ReflectingBoundary : public BoundaryType
{
public:
  virtual std::string id() const override final
  {
    return "reflecting_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new ReflectingBoundary();
  }
};


class AbsorbingBoundary : public BoundaryType
{
public:
  virtual std::string id() const override final
  {
    return "absorbing_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new AbsorbingBoundary();
  }
};


class InflowBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "inflow_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new InflowBoundary();
  }
};


class OutflowBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "outflow_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new OutflowBoundary();
  }
};


class InflowOutflowBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "inflow_outflow_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new InflowOutflowBoundary();
  }
};


class ImpermeableBoundary : public BoundaryType
{
public:
  std::string id() const override final
  {
    return "impermeable_boundary";
  }

  BoundaryType* copy() const override final
  {
    return new ImpermeableBoundary();
  }
};


#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic pop
#endif


static inline BoundaryType* make_boundary_type(const std::string& id)
{
  if (id == NoBoundary().id())
    return new NoBoundary();
  else if (id == UnknownBoundary().id())
    return new UnknownBoundary();
  else if (id == DirichletBoundary().id())
    return new DirichletBoundary();
  else if (id == NeumannBoundary().id())
    return new NeumannBoundary();
  else if (id == RobinBoundary().id())
    return new RobinBoundary();
  else if (id == ReflectingBoundary().id())
    return new ReflectingBoundary();
  else if (id == InflowBoundary().id())
    return new InflowBoundary();
  else if (id == OutflowBoundary().id())
    return new OutflowBoundary();
  else if (id == InflowBoundary().id())
    return new InflowBoundary();
  else if (id == InflowOutflowBoundary().id())
    return new InflowOutflowBoundary();
  else if (id == ImpermeableBoundary().id())
    return new ImpermeableBoundary();
  else {
    DUNE_THROW(Exceptions::boundary_type_error, "id: " << id);
    return new UnknownBoundary();
  }
} // ... make_boundary_type(...)


static constexpr AbsorbingBoundary absorbing_boundary{};
static constexpr DirichletBoundary dirichlet_boundary{};
static constexpr NoBoundary no_boundary{};
static constexpr ReflectingBoundary reflecting_boundary{};


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH
