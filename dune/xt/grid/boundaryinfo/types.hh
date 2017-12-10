// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2017 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH
#define DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


// We do not want to add a virtual destructor (to be able to use this as constexpr),
// so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#endif

class NoBoundary : public BoundaryType
{
protected:
  virtual std::string id() const override final
  {
    return "no boundary";
  }
};

class UnknownBoundary : public BoundaryType
{
protected:
  virtual std::string id() const override final
  {
    return "unknown boundary";
  }
};

class DirichletBoundary : public BoundaryType
{
protected:
  virtual std::string id() const override final
  {
    return "dirichlet boundary";
  }
};

class NeumannBoundary : public BoundaryType
{
protected:
  virtual std::string id() const override final
  {
    return "neumann boundary";
  }
};

class RobinBoundary : public BoundaryType
{
protected:
  virtual std::string id() const override final
  {
    return "robin boundary";
  }
};

class ReflectingBoundary : public BoundaryType
{
protected:
  virtual std::string id() const override final
  {
    return "reflecting boundary";
  }
};

#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#pragma GCC diagnostic pop
#endif


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH
