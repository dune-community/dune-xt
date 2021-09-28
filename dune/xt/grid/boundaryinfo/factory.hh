// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018, 2020)
//   René Fritze     (2016, 2018 - 2020)
//   Tobias Leibner  (2020)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_FACTORY_HH
#define DUNE_XT_GRID_BOUNDARYINFO_FACTORY_HH

#include "interfaces.hh"
#include "alldirichlet.hh"
#include "allneumann.hh"
#include "boundarysegment.hh"
#include "normalbased.hh"

namespace Dune::XT::Grid {


template <class I>
class BoundaryInfoFactory
{
public:
  static std::vector<std::string> available()
  {
    return {AllDirichletBoundaryInfo<I>::static_id(),
            AllNeumannBoundaryInfo<I>::static_id(),
            NormalBasedBoundaryInfo<I>::static_id(),
            BoundarySegmentIndexBasedBoundaryInfo<I>::static_id()};
  }

  static Common::Configuration default_config(const std::string type)
  {
    if (type == AllDirichletBoundaryInfo<I>::static_id())
      return alldirichlet_boundaryinfo_default_config();
    else if (type == AllNeumannBoundaryInfo<I>::static_id())
      return allneumann_boundaryinfo_default_config();
    else if (type == BoundarySegmentIndexBasedBoundaryInfo<I>::static_id())
      return boundarysegment_boundaryinfo_default_config();
    else if (type == NormalBasedBoundaryInfo<I>::static_id())
      return normalbased_boundaryinfo_default_config<I::Entity::dimension>();
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << BoundaryInfo<I>::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<BoundaryInfo<I>> create(const Common::Configuration& config)
  {
    const auto type = config.get<std::string>("type");
    if (type == AllDirichletBoundaryInfo<I>::static_id())
      return make_alldirichlet_boundaryinfo<I>(config);
    else if (type == AllNeumannBoundaryInfo<I>::static_id())
      return make_allneumann_boundaryinfo<I>(config);
    else if (type == BoundarySegmentIndexBasedBoundaryInfo<I>::static_id())
      return make_boundarysegment_boundaryinfo<I>(config);
    else if (type == NormalBasedBoundaryInfo<I>::static_id())
      return make_normalbased_boundaryinfo<I>(config);
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << BoundaryInfo<I>::static_id() << "!");
  } // ... create(...)

  static std::unique_ptr<BoundaryInfo<I>> create(const std::string& type = available().at(0))
  {
    return create(default_config(type));
  }
}; // class BoundaryInfoFactory


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_BOUNDARYINFO_FACTORY_HH
