// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff
// Copysubtrahend holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_COMBINED_HH
#define DUNE_STUFF_FUNCTIONS_COMBINED_HH

#include <type_traits>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


template <class MinuendType, class SubtrahendType>
class Difference
    : public LocalizableFunctionInterface<typename MinuendType::EntityType, typename MinuendType::DomainFieldType,
                                          MinuendType::dimDomain, typename MinuendType::RangeFieldType,
                                          MinuendType::dimRange, MinuendType::dimRangeCols>
{
  typedef LocalizableFunctionInterface<typename MinuendType::EntityType, typename MinuendType::DomainFieldType,
                                       MinuendType::dimDomain, typename MinuendType::RangeFieldType,
                                       MinuendType::dimRange, MinuendType::dimRangeCols> BaseType;
  typedef Difference<MinuendType, SubtrahendType> ThisType;

  class LocalFunction
      : public LocalfunctionInterface<typename MinuendType::EntityType, typename MinuendType::DomainFieldType,
                                      MinuendType::dimDomain, typename MinuendType::RangeFieldType,
                                      MinuendType::dimRange, MinuendType::dimRangeCols>
  {
    typedef LocalfunctionInterface<typename MinuendType::EntityType, typename MinuendType::DomainFieldType,
                                   MinuendType::dimDomain, typename MinuendType::RangeFieldType, MinuendType::dimRange,
                                   MinuendType::dimRangeCols> BaseType;

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

  private:
    static_assert(std::is_base_of<LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                               dimRange, dimRangeCols>,
                                  MinuendType>::value,
                  "MinuendType has to be derived from LocalizableFunctionInterface!");
    static_assert(std::is_same<EntityType, typename SubtrahendType::EntityType>::value, "Types do not match!");
    static_assert(std::is_same<DomainFieldType, typename SubtrahendType::DomainFieldType>::value,
                  "Types do not match!");
    static_assert(dimDomain == SubtrahendType::dimDomain, "Dimensions do not match!");
    static_assert(std::is_same<RangeFieldType, typename SubtrahendType::RangeFieldType>::value, "Types do not match!");
    static const unsigned int rRS = SubtrahendType::dimRange;
    static_assert(dimRange == SubtrahendType::dimRange, "Dimensions do not match!");
    static const unsigned int rCS = SubtrahendType::dimRangeCols;
    static_assert(dimRangeCols == SubtrahendType::dimRangeCols, "Dimensions do not match!");
    static_assert(
        std::is_base_of<LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, rRS, rCS>,
                        SubtrahendType>::value,
        "SubtrahendType has to be derived from LocalizableFunctionInterface!");
    typedef typename MinuendType::LocalfunctionType LocalMinuendType;
    typedef typename SubtrahendType::LocalfunctionType LocalSubtrahendType;

  public:
    LocalFunction(const MinuendType& minuend, const SubtrahendType& subtrahend, const EntityType& ent)
      : BaseType(ent)
      , local_minuend_(minuend.local_function(this->entity()))
      , local_subtrahend_(subtrahend.local_function(this->entity()))
      , tmp_range_(RangeFieldType(0))
      , tmp_jacobian_range_(RangeFieldType(0))
    {
    }

    LocalFunction(const LocalFunction& /*other*/) = delete;

    LocalFunction& operator=(const LocalFunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE
    {
      return std::max(local_minuend_->order(), local_subtrahend_->order());
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
    {
      local_minuend_->evaluate(xx, ret);
      local_subtrahend_->evaluate(xx, tmp_range_);
      ret -= tmp_range_;
    } // ... evaluate(...)

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
    {
      local_minuend_->jacobian(xx, ret);
      local_subtrahend_->jacobian(xx, tmp_jacobian_range_);
      ret -= tmp_jacobian_range_;
    } // ... jacobian(...)

  private:
    const std::shared_ptr<const LocalMinuendType> local_minuend_;
    const std::shared_ptr<const LocalSubtrahendType> local_subtrahend_;
    mutable RangeType tmp_range_;
    mutable JacobianRangeType tmp_jacobian_range_;
  }; // class LocalFunction

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

  typedef typename BaseType::LocalfunctionType LocalfunctionType;

public:
  Difference(const MinuendType& minuend, const SubtrahendType& subtrahend)
    : minuend_(minuend)
    , subtrahend_(subtrahend)
  {
  }

  Difference(const ThisType& other)
    : minuend_(other.minuend_)
    , subtrahend_(other.subtrahend_)
  {
  }

  ThisType& operator=(const ThisType& /*other*/) = delete;

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE
  {
    return std::unique_ptr<LocalFunction>(new LocalFunction(minuend_, subtrahend_, entity));
  }

  virtual ThisType* copy() const DS_OVERRIDE
  {
    return new ThisType(minuend_, subtrahend_);
  }

  virtual std::string name() const DS_OVERRIDE
  {
    return "difference between '" + minuend_.name() + "' and '" + subtrahend_.name() + "'";
  }

private:
  const MinuendType& minuend_;
  const SubtrahendType& subtrahend_;
}; // class Difference


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_COMBINED_HH
