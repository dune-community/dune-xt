// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH
#define DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH

#include "interfaces.hh"
#include "types.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class IntersectionImp>
class NormalBasedBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;

public:
  using typename BaseType::IntersectionType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".normalbased";
  }

  virtual BoundaryType type(const IntersectionType& /*intersection*/) const override final
  {
    DUNE_THROW(NotImplemented, "");
  }

}; // class NormalBasedBoundaryInfo


template <class I>
std::unique_ptr<NormalBasedBoundaryInfo<I>> make_normalbased_boundaryInfo(const Common::Configuration& /*config*/)
{
  return XT::Common::make_unique<NormalBasedBoundaryInfo<I>>(/*...*/);
}

Common::Configuration normalbased_boundaryinfo_default_config()
{
  DUNE_THROW(NotImplemented, "");
  return Common::Configuration("type", "xt.grid.boundaryinfo.normalbased");
}

} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_NORMALBASED_HH
