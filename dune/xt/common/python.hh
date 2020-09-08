// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)
//   René Fritze     (2018 - 2019)
//   Tobias Leibner  (2020)
//
// Created by r_milk01 on 4/25/18.

#ifndef DUNE_XT_COMMON_PYTHON_HH
#define DUNE_XT_COMMON_PYTHON_HH

#include <functional>

#include <dune/pybindxi/pybind11.h>

namespace Dune {
namespace XT {
namespace Common {
namespace bindings {


void guarded_bind(const std::function<void()>& registrar);


[[deprecated("This is not required any more (08.08.2019)!")]] void
add_initialization(pybind11::module& /*m*/, std::string /*logger_name*/, std::string /*so_name*/ = "");


[[deprecated("use guarded_bind() instead (08.08.2019)!")]] void
try_register(pybind11::module& m, const std::function<void(pybind11::module&)>& registrar);


using Entry = std::unique_ptr<pybind11::object>;
using Registry = std::unordered_map<std::string, Entry>;

namespace internal {
inline static Registry& registry()
{
  static Registry registry_;
  return registry_;
}
} // namespace internal

template <class PyBindClass>
inline static PyBindClass& from_registry(pybind11::module& m, std::string name, std::string python_name)
{
  // name = typeid(PyBindClass)
  auto iter = internal::registry().emplace(name, new PyBindClass(m, name.c_str(), python_name.c_str()));
  // if( entry.second )
  //  echo a warning if already present
  return static_cast<PyBindClass&>(*iter.first->second);
}

} // namespace bindings
} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_COMMON_PYTHON_HH
