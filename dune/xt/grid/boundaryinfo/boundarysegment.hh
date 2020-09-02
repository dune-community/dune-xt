// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2017, 2019 - 2020)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_BOUNDARYSEGMENT_HH
#define DUNE_XT_GRID_BOUNDARYINFO_BOUNDARYSEGMENT_HH

#include <dune/xt/common/configuration.hh>

#include "interfaces.hh"
#include <dune/xt/grid/boundaryinfo/types.hh>

namespace Dune {
namespace XT {
namespace Grid {


static inline Common::Configuration boundarysegment_boundaryinfo_default_config()
{
  Common::Configuration config;
  config["type"] = "xt.grid.boundaryinfo.boundarysegmentindexbased";
  config["default"] = "dirichlet";
  config["neumann"] = "[5 7]"; // all ids from 5 to 7-1
  config["robin"] = "[1 3]";
  return config;
}


// We do not want to add a virtual destructor (to be able to use this as constexpr),
// so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#  pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
#endif
template <class IntersectionImp>
class BoundarySegmentIndexBasedBoundaryInfo : public BoundaryInfo<IntersectionImp>
{
  typedef BoundaryInfo<IntersectionImp> BaseType;
  typedef BoundarySegmentIndexBasedBoundaryInfo<IntersectionImp> ThisType;

public:
  using BaseType::dimDomain;
  using BaseType::dimWorld;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::WorldType;

  static std::string static_id()
  {
    return boundarysegment_boundaryinfo_default_config().template get<std::string>("type");
  }

  static std::unique_ptr<ThisType>
  create(const Common::Configuration cfg = boundarysegment_boundaryinfo_default_config())
  {
    const Common::Configuration default_cfg = boundarysegment_boundaryinfo_default_config();
    // get default
    const std::string default_type = cfg.get("default", default_cfg.get<std::string>("default"));
    if (default_type != "dirichlet" && default_type != "neumann" && default_type != "robin")
      DUNE_THROW(Common::Exceptions::configuration_error, "Wrong default '" << default_type << "' given!");
    // get ids
    std::map<std::string, std::pair<size_t, size_t>> id_map;
    for (auto&& key : {"dirichlet", "neumann", "robin"}) {
      if (cfg.has_key(key)) {
        auto id_ranges = cfg.get<std::vector<size_t>>(key);
        if (id_ranges.size() != 2)
          DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                     "Given id range has to be of lenggth 2, is " << id_ranges);
        id_map[key] = std::make_pair(id_ranges.at(0), id_ranges.at(1));
      }
    }
    // return
    return std::make_unique<ThisType>(default_type, id_map);
  } // ... create(...)

  BoundarySegmentIndexBasedBoundaryInfo(const std::string& def,
                                        const std::map<std::string, std::pair<size_t, size_t>>& id_map)
    : default_(def)
    , id_map_(id_map)
  {}

  const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary() || intersection.neighbor())
      return no_boundary;
    const size_t boundary_segment_index = intersection.boundarySegmentIndex();
    std::string id = default_;
    for (const auto& element : id_map_)
      if (element.second.first <= boundary_segment_index && boundary_segment_index < element.second.second)
        id = element.first;
    if (id == "dirichlet")
      return dirichlet_boundary;
    else if (id == "neumann")
      return neumann_boundary;
    else if (id == "robin")
      return robin_boundary;
    else
      return unknown_boundary;
  } // ... type(...)

private:
  const std::string default_;
  const std::map<std::string, std::pair<size_t, size_t>> id_map_;
}; // class BoundarySegmentIndexBasedBoundaryInfo
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic pop
#endif

template <class I>
std::unique_ptr<BoundarySegmentIndexBasedBoundaryInfo<I>>
make_boundarysegment_boundaryinfo(const Common::Configuration& cfg = boundarysegment_boundaryinfo_default_config())
{
  return BoundarySegmentIndexBasedBoundaryInfo<I>::create(cfg);
}


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_BOUNDARYSEGMENT_HH
