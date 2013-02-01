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

#include <dune/stuff/common/color.hh>

//#if HAVE_EIGEN
//  #include <Eigen/Core>
//#endif // HAVE_EIGEN

namespace Dune {
namespace Stuff {
namespace Function {


// forward
template <class ParamFieldImp, int maxParamDim, class RangeFieldImp>
class Coefficient;


template <class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, class ParamFieldImp = double,
          int maxParamDim                                                                             = 0>
class Interface
{
public:
  typedef Interface<DomainFieldImp, domainDim, RangeFieldImp, rangeDim, ParamFieldImp, maxParamDim> ThisType;

  typedef DomainFieldImp DomainFieldType;
  static const int dimDomain = domainDim;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef RangeFieldImp RangeFieldType;
  static const int dimRange = rangeDim;
  typedef Dune::FieldVector<RangeFieldType, dimRange> RangeType;

  typedef ParamFieldImp ParamFieldType;
  static const int maxDimParam = maxParamDim;
  typedef Dune::DynamicVector<ParamFieldType> ParamType;
  typedef typename ParamType::size_type size_type;

  typedef Interface<DomainFieldType, dimDomain, RangeFieldType, dimRange> ComponentType;
  typedef Coefficient<ParamFieldType, maxParamDim, RangeFieldType> CoefficientType;

  /** \defgroup purevirtual-type ´´These methods have to be implemented and determine the type of the function (and
   * also, which of the below methods have to be implemented).'' */
  /* @{ */
  virtual bool parametric() const = 0;

  virtual bool separable() const = 0;
  /* @} */

  /** \defgroup purevirtual ´´These methods have to be implemented.'' */
  /* @{ */
  virtual std::string name() const = 0;

  virtual int order() const = 0;
  /* @} */

  /** \defgroup nonparametric-must ´´These methods have to be implemented, if parametric() == false.'' */
  /* @{ */
  virtual void evaluate(const DomainType& /*_x*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == false!");
  }
  /* @} */

  /** \defgroup nonparametric-can ´´These methods are provided by the interface, but not optimal. You can provide an
   * implementation if parametric() == false'' */
  /* @{ */
  virtual void evaluate(const ParamType& /*_mu*/, RangeType& /*_ret*/) const
  {
    assert(false && "I should be implemented!");
  }
  /* @} */

  /** \defgroup parametric-must ´´These methods have to be implemented, if parametric() == true.'' */
  /* @{ */
  virtual size_type paramSize() const
  {
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

  virtual void evaluate(const DomainType& /*_x*/, const ParamType& /*_mu*/, RangeType& /*_ret*/) const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if parametric() == true!");
  }
  /* @} */

  /** \defgroup separable ´´These methods have to be implemented, if separable() == true.'' */
  /* @{ */
  virtual size_type numComponents() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if separable() == true!");
  }

  virtual const std::vector<Dune::shared_ptr<const ComponentType>>& components() const
  {
    DUNE_THROW(Dune::NotImplemented,
               "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " implement me if separable() == true!");
  }

  virtual unsigned int numCoefficients() const
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
    RangeType ret;
    evaluate(_x, ret);
    return ret;
  }

  //  virtual RangeType evaluate(const DomainType& x, const ParamType& mu) const
  //  {
  //    RangeType ret;
  //    evaluate(x, mu, ret);
  //    return ret;
  //  }

  //  virtual void evaluate(const ParamType& x, const ParamType& mu, RangeType& ret) const
  //  {
  //    // process input
  //    assert(x.size() == dimDomain);
  //    DomainType x_fvector;
  //    for (int i = 0; i < dimDomain; ++i)
  //      x_fvector[i] = x(i);
  //    // evaluate
  //    evaluate(x_fvector, mu, ret);
  //  }

  //  virtual void evaluate(const DomainType& x, const ParamType& mu, ParamType& ret) const
  //  {
  //    // evaluate
  //    RangeType ret_fvector;
  //    evaluate(x, mu, ret_fvector);
  //    // process output
  //    assert(ret.size() == dimRange);
  //    for (int i = 0; i < dimRange; ++i)
  //      ret(i) = ret_fvector[i];
  //  }

  //  virtual void evaluate(const ParamType& x, const ParamType& mu, ParamType& ret) const
  //  {
  //    // process input
  //    assert(x.size() == dimDomain);
  //    DomainType x_fvector;
  //    for (int i = 0; i < dimDomain; ++i)
  //      x_fvector[i] = x(i);
  //    // evaluate
  //    RangeType ret_fvector;
  //    evaluate(x_fvector, mu, ret_fvector);
  //    // process output
  //    assert(ret.size() == dimRange);
  //    for (int i = 0; i < dimRange; ++i)
  //      ret(i) = ret_fvector[i];
  //  }

  //  virtual ParamType evaluate(const ParamType& x, const ParamType& mu) const
  //  {
  //    ParamType ret;
  //    evaluate(x, mu, ret);
  //    return ret;
  //  }

  //  void report(std::ostream& out = std::cout, std::string prefix = "") const
  //  {
  //    out << prefix << "parameter explanation:" << std::endl;
  //    assert(paramExplanation().size() == paramSize());
  //    assert(paramRange().size() == 2);
  //    assert(paramRange()[0].size() == paramSize());
  //    assert(paramRange()[1].size() == paramSize());
  //    for (unsigned int pp = 0; pp < paramSize(); ++pp)
  //      out << prefix << "  " << paramExplanation()[pp] << ", between " << paramRange()[0](pp) << " and " <<
  //      paramRange()[1](pp) << std::endl;
  //  }
  /* @} */
}; // class Interface

} // namespace Function
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTION_INTERFACE_HH
