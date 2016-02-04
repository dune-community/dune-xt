// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_ALLNEUMANN_HH
#define DUNE_XT_GRID_BOUNDARYINFO_ALLNEUMANN_HH

#include "interfaces.hh"
#include "types.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class IntersectionImp>
class AllNeumannBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".allneumann";
  }

  virtual BoundaryType type(const IntersectionType& intersection) const override final
  {
    if (intersection.boundary())
      return NeumannBoundary();
    else
      return UnknownBoundary();
  }
}; // class AllNeumannBoundaryInfo


template <class I>
std::unique_ptr<AllNeumannBoundaryInfo<I>> make_allneumann_boundaryInfo()
{
  return AllNeumannBoundaryInfo<I>();
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_ALLNEUMANN_HH
