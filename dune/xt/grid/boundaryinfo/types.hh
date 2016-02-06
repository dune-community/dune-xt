// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016)
//   Rene Milk       (2016)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH
#define DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {

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


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_TYPES_HH
