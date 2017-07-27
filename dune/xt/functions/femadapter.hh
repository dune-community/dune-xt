// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2017 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2013 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_FUNCTIONS_FEMADAPTER
#define DUNE_XT_FUNCTIONS_FEMADAPTER

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

#if HAVE_DUNE_FEM

namespace Dune {
namespace XT {
namespace Functions {


template <class DiscreteFunctionType>
class FemAdapterFunction
    : public LocalizableFunctionInterface<typename DiscreteFunctionType::EntityType,
                                          typename DiscreteFunctionType::DomainFieldType,
                                          DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain,
                                          typename DiscreteFunctionType::RangeFieldType,
                                          DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange,
                                          1>
{
  typedef LocalizableFunctionInterface<typename DiscreteFunctionType::EntityType,
                                       typename DiscreteFunctionType::DomainFieldType,
                                       DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain,
                                       typename DiscreteFunctionType::RangeFieldType,
                                       DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange,
                                       1>
      BaseType;
  typedef FemAdapterFunction<DiscreteFunctionType> ThisType;

  class Localfunction : public LocalfunctionInterface<typename DiscreteFunctionType::EntityType,
                                                      typename DiscreteFunctionType::DomainFieldType,
                                                      DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain,
                                                      typename DiscreteFunctionType::RangeFieldType,
                                                      DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange,
                                                      1>
  {
    typedef LocalfunctionInterface<typename DiscreteFunctionType::EntityType,
                                   typename DiscreteFunctionType::DomainFieldType,
                                   DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain,
                                   typename DiscreteFunctionType::RangeFieldType,
                                   DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange,
                                   1>
        BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const DiscreteFunctionType& df, const EntityType& ent)
      : BaseType(ent)
      , wrapped_localfunction_(df.localFunction(ent))
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override
    {
      return wrapped_localfunction_.order();
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override
    {
      wrapped_localfunction_.evaluate(xx, ret);
    }

    virtual void
    jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& /*mu*/ = {}) const override
    {
      wrapped_localfunction_.jacobian(xx, ret);
    }

  private:
    typedef typename DiscreteFunctionType::LocalFunctionType WrappedLocalfunctionType;
    WrappedLocalfunctionType wrapped_localfunction_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  FemAdapterFunction(const DiscreteFunctionType& df)
    : df_(df)
  {
  }

  static std::string static_id()
  {
    return BaseType::static_id() + ".femadapter";
  }

  virtual ThisType* copy() const override
  {
    return new ThisType(*this);
  }

  virtual std::string type() const override
  {
    return BaseType::static_id() + ".femadapter";
  }

  virtual std::string name() const override
  {
    return df_.name();
  }

  //! this intentionally hides
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return Common::make_unique<Localfunction>(df_, entity);
  } // ... local_function(...)

private:
  const DiscreteFunctionType& df_;
}; // class CheckerboardFunction


} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // if HAVE_DUNE_FEM

#endif // DUNE_XT_FUNCTIONS_FEMADAPTER
