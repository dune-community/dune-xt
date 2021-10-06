// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2021)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_FUNCTIONBASED_HH
#define DUNE_XT_GRID_BOUNDARYINFO_FUNCTIONBASED_HH

#include <unordered_map>

#include <dune/geometry/referenceelements.hh>

#include <dune/xt/grid/boundaryinfo/types.hh>
#include <dune/xt/grid/print.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace XT {
namespace Grid {


// We do not want to add a virtual destructor (to be able to use this as constexpr),
// so just silence the warning.
#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wnon-virtual-dtor"
#  pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"
#endif

template <class I>
class FunctionBasedBoundaryInfo : public BoundaryInfo<I>
{
  using BaseType = BoundaryInfo<I>;
  using ThisType = FunctionBasedBoundaryInfo;
  using E = typename I::Entity;

public:
  using FunctionType = Functions::GridFunctionInterface<E>;
  using typename BaseType::DomainType;
  using typename BaseType::IntersectionType;
  using typename BaseType::WorldType;

  using BaseType::logger;

  FunctionBasedBoundaryInfo(const BoundaryType& default_boundary_type = NoBoundary(),
                            const double tol = 1e-10,
                            const std::string& logging_prefix = "")
    : BaseType(logging_prefix.empty() ? "FunctionBasedBoundaryInfo" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , tol_(tol)
    , default_boundary_type_(default_boundary_type.copy())
  {
    LOG_(debug) << "FunctionBasedBoundaryInfo(tol=" << tol_ << ", default_boundary_type=" << *default_boundary_type_
                << ")" << std::endl;
  }

  FunctionBasedBoundaryInfo(const ThisType&) = delete;

  FunctionBasedBoundaryInfo(ThisType&& source) = default;

  void repr(std::ostream& out) const override final
  {
    out << "FunctionBasedBoundaryInfo(tol=" << tol_ << ", default_boundary_type=" << *default_boundary_type_ << ")";
  }

  void register_new_function(const FunctionType& function, const BoundaryType& boundary_type)
  {
    auto fct_ptr = function.copy_as_grid_function();
    auto bt_ptr = std::unique_ptr<BoundaryType>(boundary_type.copy());
    auto local_fct_ptr = fct_ptr->local_function();
    function_boundarytype_localfunction_list_.push_back(
        {std::move(fct_ptr), std::move(bt_ptr), std::move(local_fct_ptr)});
  } // ... register_new_function(...)

  const BoundaryType& type(const IntersectionType& intersection) const override final
  {
    LOG_(debug) << "type(intersection=" << print(intersection) << "):" << std::endl;
    if (!intersection.boundary()) {
      LOG_(debug) << "  intersection.boundary() = " << intersection.boundary() << ", returning " << no_boundary
                  << std::endl;
      return no_boundary;
    }
    const auto& inside_element = intersection.inside();
    const auto x_intersection =
        ReferenceElements<typename I::ctype, I::mydimension>::general(intersection.type()).position(0, 0);
    const auto x_inside = intersection.geometryInInside().global(x_intersection);
    LOG_(debug) << "  x_intersection = " << x_inside << ", x_inside = " << x_inside << std::endl;
    size_t counter = 0;
    for (auto& function_boundarytype_localfunction : function_boundarytype_localfunction_list_) {
      auto& local_function = std::get<2>(function_boundarytype_localfunction);
      local_function->bind(inside_element);
      const auto value = local_function->evaluate(x_inside);
      if (value > tol_) {
        const auto& boundary_type = std::get<1>(function_boundarytype_localfunction);
        LOG_(debug) << "  function " << counter << " matches (value = " << value << "), returning " << *boundary_type
                    << std::endl;
        return *boundary_type;
      }
      counter += 1;
    }
    LOG_(debug) << "  no registered function matched, returning " << *default_boundary_type_ << std::endl;
    return *default_boundary_type_;
  } // ... type(...)

private:
  const std::unique_ptr<BoundaryType> default_boundary_type_;
  const double tol_;
  std::list<std::tuple<std::unique_ptr<FunctionType>,
                       std::unique_ptr<BoundaryType>,
                       std::unique_ptr<typename FunctionType::LocalFunctionType>>>
      function_boundarytype_localfunction_list_;
}; // class FunctionBasedBoundaryInfo

#if (defined(BOOST_CLANG) && BOOST_CLANG) || (defined(BOOST_GCC) && BOOST_GCC)
#  pragma GCC diagnostic pop
#endif


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_FUNCTIONBASED_HH
