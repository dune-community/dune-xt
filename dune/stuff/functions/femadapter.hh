// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTION_FEMADAPTER
#define DUNE_STUFF_FUNCTION_FEMADAPTER

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/aliases.hh>

#if HAVE_DUNE_FEM

namespace Dune {
namespace Stuff {
namespace Functions {


template <class DiscreteFunctionType>
class FemAdapter
    : public LocalizableFunctionInterface<
          typename DiscreteFunctionType::EntityType, typename DiscreteFunctionType::DomainFieldType,
          DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain, typename DiscreteFunctionType::RangeFieldType,
          DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange, 1>
{
  typedef LocalizableFunctionInterface<
      typename DiscreteFunctionType::EntityType, typename DiscreteFunctionType::DomainFieldType,
      DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain, typename DiscreteFunctionType::RangeFieldType,
      DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange, 1> BaseType;
  typedef FemAdapter<DiscreteFunctionType> ThisType;

  class Localfunction
      : public LocalfunctionInterface<
            typename DiscreteFunctionType::EntityType, typename DiscreteFunctionType::DomainFieldType,
            DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain, typename DiscreteFunctionType::RangeFieldType,
            DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange, 1>
  {
    typedef LocalfunctionInterface<
        typename DiscreteFunctionType::EntityType, typename DiscreteFunctionType::DomainFieldType,
        DiscreteFunctionType::DiscreteFunctionSpaceType::dimDomain, typename DiscreteFunctionType::RangeFieldType,
        DiscreteFunctionType::DiscreteFunctionSpaceType::dimRange, 1> BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const unsigned int dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const unsigned int dimRange     = BaseType::dimRange;
    static const unsigned int dimRangeCols = BaseType::dimRangeCols;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

    Localfunction(const DiscreteFunctionType& df, const EntityType& ent)
      : BaseType(ent)
      , wrapped_localfunction_(df.localFunction(ent))
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE
    {
      return wrapped_localfunction_.order();
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
    {
      wrapped_localfunction_.evaluate(xx, ret);
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
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

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const int dimRange     = BaseType::dimRange;
  static const int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;

  FemAdapter(const DiscreteFunctionType& df)
    : df_(df)
  {
  }

  static std::string static_id()
  {
    return BaseType::static_id() + ".femadapter";
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(*this);
  }

  virtual std::string type() const DS_OVERRIDE
  {
    return BaseType::static_id() + ".femadapter";
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return df_.name();
  }

  //! this intentionally hides
  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const
  {
    return DSC::make_unique<Localfunction>(df_, entity);
  } // ... local_function(...)

private:
  const DiscreteFunctionType& df_;
}; // class Checkerboard


} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // if HAVE_DUNE_FEM

#endif // DUNE_STUFF_FUNCTION_FEMADAPTER
