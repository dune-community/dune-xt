// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2013 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_FACTORIES_GRID_FUNCTIONS_HH
#define DUNE_XT_FACTORIES_GRID_FUNCTIONS_HH

#include <dune/xt/common/configuration.hh>

#include <dune/xt/functions/interfaces/grid-function.hh>

#include "dune/xt/functions/checkerboard.hh"
#include "dune/xt/functions/ESV2007.hh"
#include "dune/xt/functions/spe10/model1.hh"
#include "dune/xt/functions/spe10/model2.hh"
#include "dune/xt/functions/indicator.hh"


namespace Dune {
namespace XT {
namespace Functions {


/**
 * \note If you want to add a new local function FooBar, do the following: provide a definition that is available for
all
 *       template arguments, like:
\code
template< class E, size_t r, size_t rC = 1, class R = double >
class FooBar
  : public GridFunctionInterface< E, r, rC, R >
{
  FooBar() { static_assert(AlwaysFalse< E >::value, "Not available for these dimensions!"); }
};
\endcode
 *       Every specialization that can be provided by the provider then has to define:
\code
static const bool available = true;
\endcode
 *       In addition you have to add the appropriate include
 *       in this file (of course) and the appropriate type below (just like the rest, should be obvious).
 */
template <class E, size_t r = 1, size_t rC = 1, class R = double>
class FunctionsFactory
{
public:
  using InterfaceType = GridFunctionInterface<E, r, rC, R>;

private:
  template <class FunctionType, bool available = false>
  struct Call
  {
    static std::vector<std::string> append(std::vector<std::string> in)
    {
      return in;
    }

    static bool compare(const std::string& /*type*/)
    {
      return false;
    }

    static Common::Configuration default_config(const std::string /*sub_name*/)
    {
      DUNE_THROW(Common::Exceptions::internal_error, "This should not happen!");
      return Common::Configuration();
    }

    static std::unique_ptr<FunctionType> create(const Common::Configuration& /*cfg*/)
    {
      DUNE_THROW(Common::Exceptions::internal_error, "This should not happen!");
      return std::unique_ptr<FunctionType>(nullptr);
    }
  }; // struct Call

  template <class FunctionType>
  struct Call<FunctionType, true>
  {
    static std::vector<std::string> append(std::vector<std::string> in)
    {
      in.push_back(FunctionType::static_id());
      return in;
    }

    static bool compare(const std::string& type)
    {
      return type == FunctionType::static_id();
    }

    static Common::Configuration default_config(const std::string sub_name)
    {
      return FunctionType::default_config(sub_name);
    }

    static std::unique_ptr<FunctionType> create(const Common::Configuration& cfg)
    {
      if (cfg.empty())
        return FunctionType::create();
      else
        return FunctionType::create(cfg);
    }
  }; // struct Call< ..., true >

  template <class F>
  static std::vector<std::string> call_append(std::vector<std::string> in)
  {
    return Call<F, F::available>::append(in);
  }

  template <class F>
  static bool call_compare(const std::string& type)
  {
    return Call<F, F::available>::compare(type);
  }

  template <class F>
  static Common::Configuration call_default_config(const std::string sub_name)
  {
    return Call<F, F::available>::default_config(sub_name);
  }

  template <class F>
  static std::unique_ptr<F> call_create(const Common::Configuration& cfg)
  {
    return Call<F, F::available>::create(cfg);
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

  using IndicatorType = Functions::IndicatorFunction<E, r, rC, R>;
  using CheckerboardType = Functions::CheckerboardFunction<E, r, rC, R>;
  using ESV2007CutoffType = Functions::ESV2007::CutoffFunction<E, R>;
  using Spe10Model1Type = Functions::Spe10::Model1Function<E, r, rC, R>;
  using Spe10Model2Type = Functions::Spe10::Model2Function<E, r, rC, R>;

public:
  static std::vector<std::string> available()
  {
    std::vector<std::string> ret;
    ret = call_append<CheckerboardType>(ret);
    ret = call_append<ESV2007CutoffType>(ret);
    ret = call_append<Spe10Model1Type>(ret);
    ret = call_append<Spe10Model2Type>(ret);
    ret = call_append<IndicatorType>(ret);
    return ret;
  } // ... available(...)

  static Common::Configuration default_config(const std::string type, const std::string sub_name = "")
  {
    if (call_compare<IndicatorType>(type))
      return call_default_config<IndicatorType>(sub_name);
    else if (call_compare<CheckerboardType>(type))
      return call_default_config<CheckerboardType>(sub_name);
    else if (call_compare<ESV2007CutoffType>(type))
      return call_default_config<ESV2007CutoffType>(sub_name);
    else if (call_compare<Spe10Model1Type>(type))
      return call_default_config<Spe10Model1Type>(sub_name);
    else if (call_compare<Spe10Model2Type>(type))
      return call_default_config<Spe10Model2Type>(sub_name);

    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no " << InterfaceType::static_id() << " available for dimensions " << size_t(E::d) << " -> "
                                << size_t(r)
                                << " x "
                                << size_t(rC)
                                << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for dimensions " << size_t(E::d) << " -> "
                                    << size_t(r)
                                    << " x "
                                    << size_t(rC)
                                    << ":\n"
                                    << available_as_str());
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const std::string type = available()[0],
                                               const Common::Configuration cfg = Common::Configuration())
  {
    if (call_compare<IndicatorType>(type))
      return call_create<IndicatorType>(cfg);
    else if (call_compare<CheckerboardType>(type))
      return call_create<CheckerboardType>(cfg);
    else if (call_compare<ESV2007CutoffType>(type))
      return call_create<ESV2007CutoffType>(cfg);
    else if (call_compare<Spe10Model1Type>(type))
      return call_create<Spe10Model1Type>(cfg);
    else if (call_compare<Spe10Model2Type>(type))
      return call_create<Spe10Model2Type>(cfg);
    else if (available().empty())
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "There is no " << InterfaceType::static_id() << " available for dimensions " << size_t(E::d) << " -> "
                                << size_t(r)
                                << " x "
                                << size_t(rC)
                                << "!");
    else
      DUNE_THROW(Common::Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not one of those avaible for dimensions " << size_t(E::d) << " -> "
                                    << size_t(r)
                                    << " x "
                                    << size_t(rC)
                                    << ":\n"
                                    << available_as_str());
  } // ... create(...)

  static std::unique_ptr<InterfaceType> create(const Common::Configuration& cfg)
  {
    if (!cfg.has_key("type"))
      DUNE_THROW(Common::Exceptions::wrong_input_given, "Missing 'type' in given cfg (see below)!\n\n" << cfg);
    return create(cfg.get<std::string>("type"), cfg);
  }

}; // class FunctionsProvider


} // namespace Functions {
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_HH
