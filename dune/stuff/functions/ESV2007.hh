// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_FUNCTIONS_ESV2007_HH
#define DUNE_STUFF_FUNCTIONS_ESV2007_HH

#include <cmath>

#include <dune/geometry/referenceelements.hh>

#include "interfaces.hh"

namespace Dune {
namespace Stuff {
namespace Functions {


template <class DiffusionImp>
class DUNE_DEPRECATED_MSG("Use ESV2007::Cutoff instead!") ESV2007Cutoff
    : public LocalizableFunctionInterface<typename DiffusionImp::EntityType, typename DiffusionImp::DomainFieldType,
                                          DiffusionImp::dimDomain, typename DiffusionImp::RangeFieldType, 1, 1>
{
  static_assert(std::is_base_of<Tags::LocalizableFunction, DiffusionImp>::value,
                "DiffusionImp has to be tagged as a LocalizableFunction!");
  typedef typename DiffusionImp::EntityType E_;
  typedef typename DiffusionImp::DomainFieldType D_;
  static const unsigned int d_ = DiffusionImp::dimDomain;
  typedef typename DiffusionImp::RangeFieldType R_;
  typedef LocalizableFunctionInterface<E_, D_, d_, R_, 1> BaseType;
  typedef ESV2007Cutoff<DiffusionImp> ThisType;

  class Localfunction : public LocalfunctionInterface<E_, D_, d_, R_, 1, 1>
  {
    typedef LocalfunctionInterface<E_, D_, d_, R_, 1, 1> BaseType;

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
    template <class D, int r, int rR>
    struct Compute
    {
      static_assert(AlwaysFalse<D>::value, "Not implemented for these dimensions!");
    };

    template <class D>
    struct Compute<D, 1, 1>
    {
      static RangeFieldType min_eigenvalue_of(const D& diffusion, const EntityType& ent)
      {
        const auto local_diffusion = diffusion.local_function(ent);
        assert(local_diffusion->order() == 0);
        const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(ent.type());
        return local_diffusion->evaluate(reference_element.position(0, 0))[0];
      } // ... min_eigenvalue_of_(...)
    }; // class Compute< ..., 1, 1 >

  public:
    Localfunction(const EntityType& ent, const DiffusionImp& diffusion, const RangeFieldType poincare_constant)
      : BaseType(ent)
      , value_(0)
    {
      const RangeFieldType min_eigen_value =
          Compute<DiffusionImp, DiffusionImp::dimRange, DiffusionImp::dimRangeCols>::min_eigenvalue_of(diffusion, ent);
      assert(min_eigen_value > 0.0);
      const DomainFieldType hh = compute_diameter_of_(ent);
      value_                   = (poincare_constant * hh * hh) / min_eigen_value;
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE DS_FINAL
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret[0] = value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    static DomainFieldType compute_diameter_of_(const EntityType& ent)
    {
      DomainFieldType ret(0);
      for (int cc = 0; cc < ent.template count<dimDomain>(); ++cc) {
        const auto vertex = ent.template subEntity<dimDomain>(cc)->geometry().center();
        for (int dd = cc + 1; dd < ent.template count<dimDomain>(); ++dd) {
          const auto other_vertex = ent.template subEntity<dimDomain>(dd)->geometry().center();
          const auto diff         = vertex - other_vertex;
          ret                     = std::max(ret, diff.two_norm());
        }
      }
      return ret;
    } // ... compute_diameter_of_(...)

    RangeFieldType value_;
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

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007Cutoff";
  }

  ESV2007Cutoff(const DiffusionImp& diffusion, const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl),
                const std::string nm = static_id())
    : diffusion_(diffusion)
    , poincare_constant_(poincare_constant)
    , name_(nm)
  {
  }

  ESV2007Cutoff(const ThisType& other)
    : diffusion_(other.diffusion_)
    , poincare_constant_(other.poincare_constant_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  virtual ThisType* copy() const DS_OVERRIDE DS_FINAL
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, diffusion_, poincare_constant_));
  }

private:
  const DiffusionImp& diffusion_;
  const RangeFieldType poincare_constant_;
  std::string name_;
}; // class ESV2007Cutoff


namespace ESV2007 {


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Testcase1Force : public GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef Testcase1Force<EntityImp, DomainFieldImp, RangeFieldImp> ThisType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> BaseType;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.testcase1.force";
  }

  Testcase1Force(const size_t ord = 3, const std::string nm = static_id())
    : order_(ord)
    , name_(nm)
  {
  }

  Testcase1Force(const ThisType& /*other*/) = default;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

  virtual ThisType* copy() const DS_OVERRIDE DS_FINAL
  {
    return new ThisType(*this);
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return order_;
  }

  /**
   * \brief "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   */
  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    ret[0] = M_PI_2l * M_PIl * cos(M_PI_2l * xx[0]) * cos(M_PI_2l * xx[1]);
  }

  /**
   * \brief ["-0.25 * pi * pi * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   *         "-0.25 * pi * pi * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"]
   */
  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    const DomainFieldType pre   = -0.25 * M_PIl * M_PIl * M_PIl;
    const DomainFieldType x_arg = M_PI_2l * xx[0];
    const DomainFieldType y_arg = M_PI_2l * xx[1];
    ret[0][0]                   = pre * sin(x_arg) * cos(y_arg);
    ret[0][1]                   = pre * cos(x_arg) * sin(y_arg);
  } // ... jacobian(...)

private:
  const size_t order_;
  const std::string name_;
}; // class Testcase1Force


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Testcase1ExactSolution : public GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef Testcase1ExactSolution<EntityImp, DomainFieldImp, RangeFieldImp> ThisType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1> BaseType;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.testcase1.exactsolution";
  }

  Testcase1ExactSolution(const size_t ord = 3, const std::string nm = static_id())
    : order_(ord)
    , name_(nm)
  {
  }

  Testcase1ExactSolution(const ThisType& /*other*/) = default;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

  virtual ThisType* copy() const DS_OVERRIDE DS_FINAL
  {
    return new ThisType(*this);
  }

  virtual size_t order() const DS_OVERRIDE DS_FINAL
  {
    return order_;
  }

  /**
   * \brief "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   */
  virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    ret[0] = cos(M_PI_2l * xx[0]) * cos(M_PI_2l * xx[1]);
  }

  /**
   * \brief ["-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   *         "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"]
   */
  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
  {
    const DomainFieldType pre   = -0.5 * M_PIl;
    const DomainFieldType x_arg = M_PI_2l * xx[0];
    const DomainFieldType y_arg = M_PI_2l * xx[1];
    ret[0][0]                   = pre * sin(x_arg) * cos(y_arg);
    ret[0][1]                   = pre * cos(x_arg) * sin(y_arg);
  } // ... jacobian(...)

private:
  const size_t order_;
  const std::string name_;
}; // class Testcase1ExactSolution


template <class DiffusionImp>
class Cutoff
    : public LocalizableFunctionInterface<typename DiffusionImp::EntityType, typename DiffusionImp::DomainFieldType,
                                          DiffusionImp::dimDomain, typename DiffusionImp::RangeFieldType, 1, 1>
{
  static_assert(std::is_base_of<Tags::LocalizableFunction, DiffusionImp>::value,
                "DiffusionImp has to be tagged as a LocalizableFunction!");
  typedef typename DiffusionImp::EntityType E_;
  typedef typename DiffusionImp::DomainFieldType D_;
  static const unsigned int d_ = DiffusionImp::dimDomain;
  typedef typename DiffusionImp::RangeFieldType R_;
  typedef LocalizableFunctionInterface<E_, D_, d_, R_, 1> BaseType;
  typedef Cutoff<DiffusionImp> ThisType;

  class Localfunction : public LocalfunctionInterface<E_, D_, d_, R_, 1, 1>
  {
    typedef LocalfunctionInterface<E_, D_, d_, R_, 1, 1> BaseType;

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
    template <class D, int r, int rR>
    struct Compute
    {
      static_assert(AlwaysFalse<D>::value, "Not implemented for these dimensions!");
    };

    template <class D>
    struct Compute<D, 1, 1>
    {
      static RangeFieldType min_eigenvalue_of(const D& diffusion, const EntityType& ent)
      {
        const auto local_diffusion = diffusion.local_function(ent);
        assert(local_diffusion->order() == 0);
        const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(ent.type());
        return local_diffusion->evaluate(reference_element.position(0, 0))[0];
      } // ... min_eigenvalue_of_(...)
    }; // class Compute< ..., 1, 1 >

  public:
    Localfunction(const EntityType& ent, const DiffusionImp& diffusion, const RangeFieldType poincare_constant)
      : BaseType(ent)
      , value_(0)
    {
      const RangeFieldType min_eigen_value =
          Compute<DiffusionImp, DiffusionImp::dimRange, DiffusionImp::dimRangeCols>::min_eigenvalue_of(diffusion, ent);
      assert(min_eigen_value > 0.0);
      const DomainFieldType hh = compute_diameter_of_(ent);
      value_                   = (poincare_constant * hh * hh) / min_eigen_value;
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order() const DS_OVERRIDE DS_FINAL
    {
      return 0;
    }

    virtual void evaluate(const DomainType& xx, RangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret[0] = value_;
    }

    virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const DS_OVERRIDE DS_FINAL
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    static DomainFieldType compute_diameter_of_(const EntityType& ent)
    {
      DomainFieldType ret(0);
      for (int cc = 0; cc < ent.template count<dimDomain>(); ++cc) {
        const auto vertex = ent.template subEntity<dimDomain>(cc)->geometry().center();
        for (int dd = cc + 1; dd < ent.template count<dimDomain>(); ++dd) {
          const auto other_vertex = ent.template subEntity<dimDomain>(dd)->geometry().center();
          const auto diff         = vertex - other_vertex;
          ret                     = std::max(ret, diff.two_norm());
        }
      }
      return ret;
    } // ... compute_diameter_of_(...)

    RangeFieldType value_;
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

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.cutoff";
  }

  Cutoff(const DiffusionImp& diffusion, const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl),
         const std::string nm = static_id())
    : diffusion_(diffusion)
    , poincare_constant_(poincare_constant)
    , name_(nm)
  {
  }

  Cutoff(const ThisType& other)
    : diffusion_(other.diffusion_)
    , poincare_constant_(other.poincare_constant_)
    , name_(other.name_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  virtual ThisType* copy() const DS_OVERRIDE DS_FINAL
  {
    return new ThisType(*this);
  }

  virtual std::string name() const DS_OVERRIDE DS_FINAL
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const DS_OVERRIDE DS_FINAL
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, diffusion_, poincare_constant_));
  }

private:
  const DiffusionImp& diffusion_;
  const RangeFieldType poincare_constant_;
  std::string name_;
}; // class Cutoff


} // namespace ESV2007
} // namespace Functions
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_FUNCTIONS_ESV2007_HH
