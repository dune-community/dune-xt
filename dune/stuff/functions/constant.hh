#ifndef DUNE_STUFF_FUNCTIONS_CONSTANT_HH
#define DUNE_STUFF_FUNCTIONS_CONSTANT_HH

#include <memory>

#include <dune/stuff/common/parameter/tree.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Function {


/**
 *  \note Only implemented scalar and vector valued at the moment!
 */
template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class Constant;

template <class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class Constant<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
  typedef LocalizableFunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;
  typedef Constant<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

  class Localfunction : public LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
  {
    typedef LocalfunctionInterface<EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

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

    Localfunction(const EntityType& ent, const std::shared_ptr<const RangeType> value)
      : BaseType(ent)
      , value_(value)
    {
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret = *value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    const std::shared_ptr<const RangeType> value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;

  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  static std::string static_id();

  static Dune::ParameterTree defaultSettings(const std::string subName = "");

  static ThisType* create(const DSC::ExtendedParameterTree settings = defaultSettings());

  Constant(const RangeFieldType& val, const std::string nm = static_id());

  Constant(const RangeType& val, const std::string nm = static_id());

  Constant(const ThisType& other);

  ThisType& operator=(const ThisType& other);

  virtual ThisType* copy() const DS_OVERRIDE;

  virtual std::string name() const DS_OVERRIDE;

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE;

private:
  std::shared_ptr<const RangeType> value_;
  std::string name_;
}; // class Constant


} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_CONSTANT_HH
