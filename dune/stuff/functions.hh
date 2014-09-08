// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_HH
#define DUNE_STUFF_FUNCTIONS_HH

#include <string>
#include <vector>
#include <memory>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/configuration.hh>

#include "functions/interfaces.hh"
#include "functions/checkerboard.hh"
#include "functions/constant.hh"
#include "functions/ESV2007.hh"
#include "functions/expression.hh"
#include "functions/spe10.hh"

namespace Dune {
namespace Stuff {
namespace internal {


/**
 *  \brief  Base class containing functions that are available for all dimension.
 */
template <class E, class D, int d, class R, int r, int rC = 1>
class FunctionsProviderBase
{
public:
  typedef LocalizableFunctionInterface<E, D, d, R, r, rC> InterfaceType;

protected:
  template <class FunctionType>
  static std::unique_ptr<InterfaceType> call_create(const Common::Configuration& config)
  {
    if (config.empty())
      return FunctionType::create();
    else
      return FunctionType::create(config);
  } // ... call_create(...)

public:
  static std::vector<std::string> available()
  {
    return {Functions::Constant<E, D, d, R, r, rC>::static_id()};
  } // ... available(...)

  static Common::Configuration default_config(const std::string type = available()[0], const std::string sub_name = "")
  {
    if (type == Functions::Constant<E, D, d, R, r, rC>::static_id())
      return Functions::Constant<E, D, d, R, r, rC>::default_config(sub_name);
    else
      DUNE_THROW(Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const std::string type = available()[0],
                                               const Common::Configuration config = Common::Configuration())
  {
    if (type == Functions::Constant<E, D, d, R, r, rC>::static_id())
      return call_create<Functions::Constant<E, D, d, R, r, rC>>(config);
    else
      DUNE_THROW(Exceptions::wrong_input_given,
                 "Requested type '" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
}; // class FunctionsProviderBase


} // namespace internal


template <class E, class D, int d, class R, int r, int rC = 1>
class FunctionsProvider : public internal::FunctionsProviderBase<E, D, d, R, r, rC>
{
};


template <class E, class D, int d, class R, int r>
class FunctionsProvider<E, D, d, R, r, 1> : public internal::FunctionsProviderBase<E, D, d, R, r, 1>
{
  static const unsigned int rC = 1;
  typedef internal::FunctionsProviderBase<E, D, d, R, r, rC> BaseType;

public:
  using typename BaseType::InterfaceType;

  static std::vector<std::string> available()
  {
    auto base = BaseType::available();
    base.push_back(Functions::Expression<E, D, d, R, r, rC>::static_id());
    base.push_back(Functions::Checkerboard<E, D, d, R, r, rC>::static_id());
    return base;
  } // ... available(...)

  static Common::Configuration default_config(const std::string type = available()[0], const std::string sub_name = "")
  {
    if (type == Functions::Expression<E, D, d, R, r, rC>::static_id())
      return Functions::Expression<E, D, d, R, r, rC>::default_config(sub_name);
    else if (type == Functions::Checkerboard<E, D, d, R, r, rC>::static_id())
      return Functions::Checkerboard<E, D, d, R, r, rC>::default_config(sub_name);
    else
      return BaseType::default_config(type, sub_name);
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const std::string type = available()[0],
                                               const Common::Configuration config = Common::Configuration())
  {
    if (type == Functions::Expression<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::Expression<E, D, d, R, r, rC>>(config);
    else if (type == Functions::Checkerboard<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::Checkerboard<E, D, d, R, r, rC>>(config);
    else
      return BaseType::create(type, config);
  } // ... create(...)
}; // class FunctionsProvider< ..., 1 >


template <class E, class D, class R>
class FunctionsProvider<E, D, 2, R, 1, 1> : public internal::FunctionsProviderBase<E, D, 2, R, 1, 1>
{
  static const unsigned int d  = 2;
  static const unsigned int r  = 1;
  static const unsigned int rC = 1;
  typedef internal::FunctionsProviderBase<E, D, d, R, r, rC> BaseType;

public:
  using typename BaseType::InterfaceType;

  static std::vector<std::string> available()
  {
    auto base = BaseType::available();
    base.push_back(Functions::Expression<E, D, d, R, r, rC>::static_id());
    base.push_back(Functions::Checkerboard<E, D, d, R, r, rC>::static_id());
    base.push_back(Functions::ESV2007::Testcase1Force<E, D, d, R, r, rC>::static_id());
    base.push_back(Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r, rC>::static_id());
    base.push_back(Functions::Spe10Model1<E, D, d, R, r, rC>::static_id());
    return base;
  } // ... available(...)

  static Common::Configuration default_config(const std::string type = available()[0], const std::string sub_name = "")
  {
    if (type == Functions::Expression<E, D, d, R, r, rC>::static_id())
      return Functions::Expression<E, D, d, R, r, rC>::default_config(sub_name);
    else if (type == Functions::Checkerboard<E, D, d, R, r, rC>::static_id())
      return Functions::Checkerboard<E, D, d, R, r, rC>::default_config(sub_name);
    else if (type == Functions::ESV2007::Testcase1Force<E, D, d, R, r, rC>::static_id())
      return Functions::ESV2007::Testcase1Force<E, D, d, R, r, rC>::default_config(sub_name);
    else if (type == Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r, rC>::static_id())
      return Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r, rC>::default_config(sub_name);
    else if (type == Functions::Spe10Model1<E, D, d, R, r, rC>::static_id())
      return Functions::Spe10Model1<E, D, d, R, r, rC>::default_config(sub_name);
    else
      return BaseType::default_config(type, sub_name);
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const std::string type = available()[0],
                                               const Common::Configuration config = Common::Configuration())
  {
    if (type == Functions::Expression<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::Expression<E, D, d, R, r, rC>>(config);
    else if (type == Functions::Checkerboard<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::Checkerboard<E, D, d, R, r, rC>>(config);
    else if (type == Functions::ESV2007::Testcase1Force<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::ESV2007::Testcase1Force<E, D, d, R, r, rC>>(config);
    else if (type == Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::ESV2007::Testcase1ExactSolution<E, D, d, R, r, rC>>(config);
    else if (type == Functions::Spe10Model1<E, D, d, R, r, rC>::static_id())
      return BaseType::template call_create<Functions::Spe10Model1<E, D, d, R, r, rC>>(config);
    else
      return BaseType::create(type, config);
  } // ... create(...)
}; // class FunctionsProvider< ..., 2, 1, 1 >


} // namespace Stuff
} // namespace Dune


#endif // DUNE_STUFF_FUNCTIONS_HH
