// This file is part of the dune-xt-grid project:
//   https://github.com/dune-community/dune-xt-grid
// Copyright 2009-2018 dune-xt-grid developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2012 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2013 - 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_XT_GRID_PROVIDER_HH
#define DUNE_XT_GRID_PROVIDER_HH

#include <memory>
#include <type_traits>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/type_traits.hh>

#include "cube.hh"
#include "dgf.hh"
#include "gmsh.hh"

namespace Dune {
namespace XT {
namespace Grid {


template <class GridType>
class GridProviderFactory
{
  template <class G, bool available = false>
  struct Helper
  {
    static std::vector<std::string> append(std::vector<std::string> in)
    {
      return in;
    }

    static bool compare(const std::string& /*type*/)
    {
      return false;
    }

    static Common::Configuration default_config()
    {
      DUNE_THROW(Common::Exceptions::internal_error, "This should not happen!");
      return Common::Configuration();
    }

    static GridProvider<GridType, none_t> create(const Common::Configuration& /*cfg*/)
    {
      DUNE_THROW(Common::Exceptions::internal_error, "This should not happen!");
      return GridProvider<GridType, none_t>(nullptr);
    }
  }; // struct Helper

  template <class G>
  struct Helper<G, true>
  {
    static std::vector<std::string> append(std::vector<std::string> in)
    {
      in.push_back(G::static_id());
      return in;
    }

    static bool compare(const std::string& type)
    {
      return type == G::static_id();
    }

    static Common::Configuration default_config()
    {
      return G::default_config();
    }

    static GridProvider<GridType, none_t> create(const Common::Configuration& cfg)
    {
      if (cfg.empty())
        return G::create();
      else
        return G::create(cfg);
    }
  }; // struct Helper< ..., true >

  template <class G>
  static std::vector<std::string> call_append(std::vector<std::string> in)
  {
    return Helper<G, G::available>::append(in);
  }

  template <class G>
  static bool call_compare(const std::string& type)
  {
    return Helper<G, G::available>::compare(type);
  }

  template <class G>
  static Common::Configuration call_default_config(const std::string sub_name)
  {
    return Helper<G, G::available>::default_config(sub_name);
  }

  template <class G>
  static GridProvider<GridType, none_t> call_create(const Common::Configuration& cfg)
  {
    return Helper<G, G::available>::create(cfg);
  }

  static std::string available_as_str()
  {
    std::string ret = "";
    const auto vals = available();
    if (vals.size() > 0) {
      ret += vals[0];
      for (size_t ii = 1; ii < vals.size(); ++ii)
        ret += "\n   " + vals[ii];
    }
    return ret;
  } // ... available_as_str(...)

  typedef CubeGridProviderFactory<GridType> CubeType;
  typedef DgfGridProviderFactory<GridType> DgfType;
  typedef GmshGridProviderFactory<GridType> GmshType;

public:
  static std::vector<std::string> available()
  {
    std::vector<std::string> ret;
    call_append<CubeType>(ret);
    call_append<DgfType>(ret);
    call_append<GmshType>(ret);
    return ret;
  }

  static Common::Configuration default_config(const std::string type)
  {
    if (call_compare<CubeType>(type))
      return call_default_config<CubeType>();
    else if (call_compare<DgfType>(type))
      return call_default_config<DgfType>();
    else if (call_compare<GmshType>(type))
      return call_default_config<GmshType>();
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value()
                                    << ":\n"
                                    << available_as_str());
  } // ... default_config(...)

  static GridProvider<GridType, none_t> create(const Common::Configuration& config)
  {
    return create(config.get<std::string>("type"), config);
  }

  static GridProvider<GridType, none_t> create(const std::string& type = available()[0],
                                               const Common::Configuration config = Common::Configuration())
  {
    if (call_compare<CubeType>(type))
      return call_create<CubeType>(config);
    else if (call_compare<DgfType>(type))
      return call_create<DgfType>(config);
    else if (call_compare<GmshType>(type))
      return call_create<GmshType>(config);
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value()
                                    << ":\n"
                                    << available_as_str());
  } // ... create(...)
}; // class GridProviderFactory


template <class GridType>
class DdSubdomainGridProviderFactory
{
  template <class G, bool available = false>
  struct Helper
  {
    static std::vector<std::string> append(std::vector<std::string> in)
    {
      return in;
    }

    static bool compare(const std::string& /*type*/)
    {
      return false;
    }

    static Common::Configuration default_config()
    {
      DUNE_THROW(Common::Exceptions::internal_error, "This should not happen!");
      return Common::Configuration();
    }

    static GridProvider<GridType, DD::SubdomainGrid<GridType>> create(const Common::Configuration& /*cfg*/)
    {
      DUNE_THROW(Common::Exceptions::internal_error, "This should not happen!");
      return GridProvider<GridType, DD::SubdomainGrid<GridType>>(nullptr);
    }
  }; // struct Helper

  template <class G>
  struct Helper<G, true>
  {
    static std::vector<std::string> append(std::vector<std::string> in)
    {
      in.push_back(G::static_id());
      return in;
    }

    static bool compare(const std::string& type)
    {
      return type == G::static_id();
    }

    static Common::Configuration default_config()
    {
      return G::default_config();
    }

    static GridProvider<GridType, DD::SubdomainGrid<GridType>> create(const Common::Configuration& cfg)
    {
      if (cfg.empty())
        return G::create();
      else
        return G::create(cfg);
    }
  }; // struct Helper< ..., true >

  template <class G>
  static std::vector<std::string> call_append(std::vector<std::string> in)
  {
    return Helper<G, G::available>::append(in);
  }

  template <class G>
  static bool call_compare(const std::string& type)
  {
    return Helper<G, G::available>::compare(type);
  }

  template <class G>
  static Common::Configuration call_default_config(const std::string sub_name)
  {
    return Helper<G, G::available>::default_config(sub_name);
  }

  template <class G>
  static GridProvider<GridType, DD::SubdomainGrid<GridType>> call_create(const Common::Configuration& cfg)
  {
    return Helper<G, G::available>::create(cfg);
  }

  static std::string available_as_str()
  {
    std::string ret = "";
    const auto vals = available();
    if (vals.size() > 0) {
      ret += vals[0];
      for (size_t ii = 1; ii < vals.size(); ++ii)
        ret += "\n   " + vals[ii];
    }
    return ret;
  } // ... available_as_str(...)

  typedef CubeDdSubdomainsGridProviderFactory<GridType> CubeType;

public:
  static std::vector<std::string> available()
  {
    std::vector<std::string> ret;
    return call_append<CubeType>(ret);
  }

  static Common::Configuration default_config(const std::string type)
  {
    if (call_compare<CubeType>(type))
      return call_default_config<CubeType>();
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value()
                                    << ":\n"
                                    << available_as_str());
  } // ... default_config(...)

  static GridProvider<GridType, DD::SubdomainGrid<GridType>> create(const Common::Configuration& config)
  {
    return create(config.get<std::string>("type"), config);
  }

  static GridProvider<GridType, DD::SubdomainGrid<GridType>>
  create(const std::string& type = available()[0], const Common::Configuration config = Common::Configuration())
  {
    if (call_compare<CubeType>(type))
      return call_create<CubeType>(config);
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no grid provider available for " << Common::Typename<GridType>::value() << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for "
                                    << Common::Typename<GridType>::value()
                                    << ":\n"
                                    << available_as_str());
  } // ... create(...)
}; // class DdSubdomainGridProviderFactory


} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_GRID_PROVIDER_HH
