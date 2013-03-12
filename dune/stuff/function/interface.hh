#ifndef DUNE_STUFF_FUNCTION_INTERFACE_HH
#define DUNE_STUFF_FUNCTION_INTERFACE_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <vector>

#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/function.hh>

#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/parameter.hh>

namespace Dune {
namespace Stuff {


//! forward
template <class RangeFieldImp>
class FunctionAffineSeparablCoefficient;


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim>
class FunctionInterface
{
public:
  typedef FunctionInterface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim> ThisType;

  typedef DomainFieldImp DomainFieldType;
  static const int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const int dimRange = rangeDim;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  typedef Common::Parameter::FieldType ParamFieldType;
  static const int maxParamDim = Common::Parameter::maxDim;
  typedef Common::Parameter::Type ParamType;

  typedef ThisType ComponentType;
  typedef FunctionAffineSeparablCoefficient<RangeFieldType> CoefficientType;

  virtual ~FunctionInterface()
  {
  }

  static const std::string id()
  {
    return "function";
  }

  /** \defgroup type ´´Theis method has to be implemented for parametric functions and determines,
   *                   which evaluate() is callable.''
   */
  /* @{ */
  virtual bool parametric() const
  {
    return false;
  }
  /* @} */

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

  /** \defgroup nonparametric-must ´´These methods have to be implemented, if parametric() == false.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == false!");
  }
  /* @} */

  /** \defgroup parametric-must ´´These methods have to be implemented, if parametric() == true.'' */
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

  virtual bool separable() const
  {
    return false;
  }
  /* @} */

  /** \defgroup separable ´´These methods have to be implemented, if separable() == true.'' */
  /* @{ */
  virtual size_t numComponents() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if separable() == true!");
  }

  virtual const std::vector<Dune::shared_ptr<const ComponentType>>& components() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if separable() == true!");
  }

  virtual size_t numCoefficients() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if separable() == true!");
  }

  virtual const std::vector<Dune::shared_ptr<const CoefficientType>>& coefficients() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if separable() == true!");
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
}; // class FunctionInterface


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
