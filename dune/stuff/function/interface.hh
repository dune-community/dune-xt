#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#include <vector>
#include <memory>
#include <string>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/function/common/function.hh>
#include <dune/fem/space/common/functionspace.hh>
#endif

#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter.hh>

#include <dune/stuff/localfunction/interface.hh>

#include "local.hh"

namespace Dune {
namespace Stuff {


// forward, needed in the interface, included below
template <class RangeFieldImp>
class AffineParametricCoefficientFunction;


// forward, needed in the interface
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class FunctionInterface;


/**
 *  \brief  Interface for stationary matrix valued functions, nonparametric, parametric and affineparametric.
 *
 *          See spezialization for rangeDimCols = 1 for scalar and vector valued functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class GenericStationaryFunctionInterface : public LocalizableFunction
{
public:
  typedef GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
      ThisType;

  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  typedef Common::Parameter::FieldType ParamFieldType;
  static const unsigned int maxParamDim = Common::Parameter::maxDim;
  typedef Common::Parameter::Type ParamType;

  typedef FunctionInterface<DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols> ComponentType;
  typedef AffineParametricCoefficientFunction<RangeFieldType> CoefficientType;

  virtual ~GenericStationaryFunctionInterface()
  {
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return "function";
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup type ´´This method has to be implemented for parametric functions and determines which evaluate() is
   *            callable.''
   */
  /* @{ */
  virtual bool parametric() const
  {
    return false;
  }
  /* @} */

  /** \defgroup nonparametric-must ´´This method has to be implemented, if parametric() == false.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == false!");
  }
  /* @} */

  /** \defgroup parametric-must ´´This method has to be implemented, if parametric() == true.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const ParamType& /*_mu*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual size_t paramSize() const
  {
    if (!parametric())
      return 0;
    else
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual const std::vector<ParamType>& paramRange() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual const std::vector<std::string>& paramExplanation() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual bool affineparametric() const
  {
    return false;
  }
  /* @} */

  /** \defgroup affineparametric ´´These methods have to be implemented, if affineparametric() == true.'' */
  /* @{ */
  virtual const std::vector<std::shared_ptr<const ComponentType>>& components() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if affineparametric() == true!");
  }

  virtual const std::vector<std::shared_ptr<const CoefficientType>>& coefficients() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if affineparametric() == true!");
  }

  virtual bool hasAffinePart() const
  {
    return false;
  }
  /* @} */

  /** \defgroup affineparametric-affinepart ´´This method has to be implemented, if hasAffinePart() == true.'' */
  /* @{ */
  virtual const std::shared_ptr<const ComponentType>& affinePart() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if hasAffinePart() == true!");
  }
  /* @} */

  /** \defgroup provided ´´These methods are provided by the interface itself, but may not be implemented optimal.'' */
  /* @{ */
  virtual RangeType evaluate(const DomainType& _x) const
  {
    assert(!parametric());
    RangeType ret;
    evaluate(_x, ret);
    return ret;
  }

  virtual RangeType evaluate(const DomainType& _x, const ParamType& _mu) const
  {
    assert(parametric());
    RangeType ret;
    evaluate(_x, _mu, ret);
    return ret;
  }

  template <class EntityType>
  LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
  localFunction(const EntityType& entity) const
  {
    return LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>(
        *this, entity);
  }
  /* @} */
}; // class GenericStationaryFunctionInterface


/**
 *  \brief  Interface for stationary scalar and vector valued functions, nonparametric, parametric and affineparametric.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public LocalizableFunction
{
  typedef GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  typedef Common::Parameter::FieldType ParamFieldType;
  static const unsigned int maxParamDim = Common::Parameter::maxDim;
  typedef Common::Parameter::Type ParamType;

  typedef FunctionInterface<DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols> ComponentType;
  typedef AffineParametricCoefficientFunction<RangeFieldType> CoefficientType;

  virtual ~GenericStationaryFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return id();
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup type ´´This method has to be implemented for parametric functions and determines which evaluate() is
   *            callable.''
   */
  /* @{ */
  virtual bool parametric() const
  {
    return false;
  }
  /* @} */

  /** \defgroup nonparametric-must ´´This method has to be implemented, if parametric() == false.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == false!");
  }
  /* @} */

  /** \defgroup parametric-must ´´This method has to be implemented, if parametric() == true.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const ParamType& /*_mu*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual size_t paramSize() const
  {
    if (!parametric())
      return 0;
    else
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual const std::vector<ParamType>& paramRange() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual const std::vector<std::string>& paramExplanation() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }

  virtual bool affineparametric() const
  {
    return false;
  }
  /* @} */

  /** \defgroup affineparametric ´´These methods have to be implemented, if affineparametric() == true.'' */
  /* @{ */
  virtual const std::vector<std::shared_ptr<const ComponentType>>& components() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if affineparametric() == true!");
  }

  virtual const std::vector<std::shared_ptr<const CoefficientType>>& coefficients() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if affineparametric() == true!");
  }

  virtual bool hasAffinePart() const
  {
    return false;
  }
  /* @} */

  /** \defgroup affineparametric-affinepart ´´This method has to be implemented, if hasAffinePart() == true.'' */
  /* @{ */
  virtual const std::shared_ptr<const ComponentType>& affinePart() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if hasAffinePart() == true!");
  }
  /* @} */

  /** \defgroup provided ´´These methods are provided by the interface itself, but may not be implemented optimal.'' */
  /* @{ */
  virtual RangeType evaluate(const DomainType& _x) const
  {
    assert(!parametric());
    RangeType ret;
    evaluate(_x, ret);
    return ret;
  }

  virtual RangeType evaluate(const DomainType& _x, const ParamType& _mu) const
  {
    assert(parametric());
    RangeType ret;
    evaluate(_x, _mu, ret);
    return ret;
  }

  template <class EntityType>
  struct LocalFunction
  {
    typedef LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols> Type;
  };

  template <class EntityType>
  LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
  localFunction(const EntityType& entity) const
  {
    return LocalizedFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>(
        *this, entity);
  }
  /* @} */
}; // class GenericStationaryFunctionInterface< ..., 1 >


/**
 *  \brief  Interface for matrix valued nonparametric stationary function.
 *
 *          See specialization for rangeDimRows = 1 for scalar and vector valued functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class FunctionInterface
    : public GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
{
  typedef GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
      BaseType;

public:
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRangeRows = BaseType::dimRangeRows;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  virtual ~FunctionInterface()
  {
  }

  static std::string id()
  {
    return "function";
  }

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, RangeType& /*_ret*/) const = 0;
  /* @} */

  using BaseType::localFunction;
}; // class FunctionInterface


/**
 *  \brief  Interface for scalar and vector valued nonparametric stationary function.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
    : public GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
#if HAVE_DUNE_FEM
#if DUNE_VERSION_NEWER(DUNE_FEM, 1, 4)
      ,
      public Dune::Fem::Function<Dune::Fem::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
#else
      ,
      public Dune::Fem::Function<Dune::FunctionSpace<DomainFieldImp, RangeFieldImp, domainDim, rangeDim>,
#endif
                                 FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>>
#endif // HAVE_DUNE_FEM
{
  typedef GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1> BaseType;

public:
  typedef typename BaseType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange     = BaseType::dimRange;
  static const unsigned int dimRangeRows = BaseType::dimRangeRows;
  static const unsigned int dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  virtual ~FunctionInterface()
  {
  }

  static std::string id()
  {
    return "function";
  }

  /** \defgroup must This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, RangeType& /*_ret*/) const = 0;
  /* @} */

  using BaseType::evaluate;
  using BaseType::localFunction;
}; // class FunctionInterface


/**
 *  \brief  Interface for parametric functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class ParametricFunctionInterface
    : public GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
{
  typedef GenericStationaryFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
      BaseType;

public:
  typedef typename BaseType::DomainType DomainType;
  typedef typename BaseType::RangeType RangeType;
  typedef typename BaseType::ParamType ParamType;

  virtual ~ParametricFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function.parametric";
  }

  virtual bool parametric() const
  {
    return true;
  }

  /** \defgroup must ´´This methods have to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const ParamType& /*_mu*/, RangeType& /*_ret*/) const = 0;

  virtual size_t paramSize() const = 0;

  virtual const std::vector<ParamType>& paramRange() const = 0;

  virtual const std::vector<std::string>& paramExplanation() const = 0;
  /* @} */
}; // class ParametricFunctionInterface


/**
 *  \brief  Interface for affine parametric functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class AffineParametricFunctionInterface
    : public ParametricFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols>
{
  typedef ParametricFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDimRows, rangeDimCols> BaseType;

public:
  typedef typename BaseType::ComponentType ComponentType;
  typedef typename BaseType::CoefficientType CoefficientType;

  virtual ~AffineParametricFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function.affineparametric";
  }

  virtual bool affineparametric() const
  {
    return true;
  }

  /** \defgroup must ´´These methods have to be implemented.'' */
  /* @{ */
  virtual const std::vector<std::shared_ptr<const ComponentType>>& components() const = 0;

  virtual const std::vector<std::shared_ptr<const CoefficientType>>& coefficients() const = 0;
  /* @} */

  virtual bool hasAffinePart() const
  {
    return false;
  }

  /** \defgroup affineparametric-affinepart ´´This method has to be implemented, if hasAffinePart() == true.'' */
  /* @{ */
  virtual const std::shared_ptr<const ComponentType>& affinePart() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if hasAffinePart() == true!");
  }
  /* @} */
}; // class AffineParametricFunctionInterface


/**
 *  \brief  Interface for matrix valued nonparametric timedependent function.
 *
 *          See spezialization for rangeDimCols = 1 for scalar and vector valued functions.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDimRows, int rangeDimCols>
class TimedependentFunctionInterface
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRangeRows = rangeDimRows;
  static const unsigned int dimRangeCols = rangeDimCols;
  typedef Dune::FieldMatrix<RangeFieldType, dimRangeRows, dimRangeCols> RangeType;

  virtual ~TimedependentFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function.timedependent";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return "function.timedependent";
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const double& /*_t*/, RangeType& /*_ret*/) const = 0;
  /* @} */
}; // class TimedependentFunctionInterface


/**
 *  \brief  Interface for scalar and vector valued nonparametric timedependent function.
 */
template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class TimedependentFunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1>
{
public:
  typedef DomainFieldImp DomainFieldType;
  static const unsigned int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const unsigned int dimRange     = rangeDim;
  static const unsigned int dimRangeRows = dimRange;
  static const unsigned int dimRangeCols = 1;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  virtual ~TimedependentFunctionInterface()
  {
  }

  static std::string id()
  {
    return "function.timedependent";
  }

  /** \defgroup info ´´These methods should be implemented in order to identify the function.'' */
  /* @{ */
  virtual std::string name() const
  {
    return "function.timedependent";
  }

  virtual int order() const
  {
    return -1;
  }
  /* @} */

  /** \defgroup must ´´This method has to be implemented.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, const double& /*_t*/, RangeType& /*_ret*/) const = 0;
  /* @} */
}; // class TimedependentFunctionInterface


} // namespace Stuff
} // namespace Dune

#include "affineparametric/coefficient.hh"

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
