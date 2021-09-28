// This file is part of the dune-xt project:
//   https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-xt
// Copyright 2009-2021 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018, 2020)
//   René Fritze     (2016, 2018 - 2020)
//   Tobias Leibner  (2016 - 2020)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_INTERFACE_HH
#define DUNE_XT_GRID_BOUNDARYINFO_INTERFACE_HH

#include <string>

#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/grid/type_traits.hh>

namespace Dune::XT::Grid {

// We do not want to add a virtual destructor (to be able to use this as constexpr), so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#  pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
#endif


class BoundaryType
{
public:
  virtual std::string id() const = 0;

  virtual BoundaryType* copy() const = 0;

  virtual bool operator==(const BoundaryType& other) const
  {
    return id() == other.id();
  }

  virtual bool operator!=(const BoundaryType& other) const
  {
    return !operator==(other);
  }
}; // class BoundaryType


#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic pop
#endif


std::ostream& operator<<(std::ostream& out, const BoundaryType& type);


template <class I>
class BoundaryInfo;


template <class IntersectionImp>
class BoundaryInfo : public Common::WithLogger<BoundaryInfo<IntersectionImp>>
{
  static_assert(is_intersection<IntersectionImp>::value);
  using ThisType = BoundaryInfo<IntersectionImp>;
  using Logger = Common::WithLogger<BoundaryInfo<IntersectionImp>>;

public:
  using IntersectionType = IntersectionImp;
  using DomainFieldType = typename IntersectionType::ctype;
  static constexpr size_t dimDomain = IntersectionType::Entity::dimension;
  static constexpr size_t dimWorld = IntersectionType::dimensionworld;
  using DomainType = Common::FieldVector<DomainFieldType, dimDomain>;
  using WorldType = Common::FieldVector<DomainFieldType, dimWorld>;

  BoundaryInfo(const std::string& logging_prefix = "",
               const std::array<bool, 3>& logging_state = Common::default_logger_state())
    : Logger(logging_prefix.empty() ? "BoundaryInfo" : logging_prefix, logging_state)
  {}

  BoundaryInfo(const ThisType&) = default;

  BoundaryInfo(ThisType&&) = default;

  virtual ~BoundaryInfo() = default;

  virtual const BoundaryType& type(const IntersectionType& intersection) const = 0;

  virtual void repr(std::ostream& out) const
  {
    out << "BoundaryInfo(\?\?\?)";
  }

  virtual std::string str() const
  {
    std::stringstream ss;
    this->repr(ss);
    return ss.str();
  }

  static std::string static_id()
  {
    return "xt.grid.boundaryinfo";
  }
}; // class BoundaryInfo


template <class I>
std::ostream& operator<<(std::ostream& out, const BoundaryInfo<I>& bi)
{
  bi.repr(out);
  return out;
}


} // namespace Dune::XT::Grid

#endif // DUNE_XT_GRID_BOUNDARYINFO_INTERFACE_HH
