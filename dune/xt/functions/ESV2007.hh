// This file is part of the dune-xt-functions project:
//   https://github.com/dune-community/dune-xt-functions
// Copyright 2009-2018 dune-xt-functions developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014 - 2016, 2018)
//   Tobias Leibner  (2014 - 2017)

#ifndef DUNE_XT_FUNCTIONS_ESV2007_HH
#define DUNE_XT_FUNCTIONS_ESV2007_HH

#include <cmath>
#include <type_traits>

#if HAVE_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <dune/geometry/referenceelements.hh>

#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/debug.hh>
#include <dune/xt/common/ranges.hh>
#if HAVE_DUNE_XT_LA
#include <dune/xt/la/container/eigen.hh>
#endif

#include <dune/xt/functions/interfaces.hh>

namespace Dune {
namespace XT {
namespace Functions {
namespace ESV2007 {


template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class Testcase1Force : public LocalizableFunctionInterface<E, D, d, R, r, rC>
{
  Testcase1Force()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Testcase1Force<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1>
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1>
{
  typedef Testcase1Force<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1> ThisType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1> BaseType;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.testcase1.force";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["integration_order"] = "3";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // create
    return Common::make_unique<ThisType>(cfg.get("integration_order", default_cfg.get<size_t>("integration_order")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  Testcase1Force(const size_t ord = default_config().template get<size_t>("integration_order"),
                 const std::string nm = static_id())
    : order_(ord)
    , name_(nm)
  {
  }

  Testcase1Force(const ThisType& /*other*/) = default;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".ESV2007.testcase1.force";
  }

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return order_;
  }

  /**
   * \brief "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   */
  virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    ret[0] = M_PI_2l * M_PIl * cos(M_PI_2l * xx[0]) * cos(M_PI_2l * xx[1]);
  }

  /**
   * \brief ["-0.25 * pi * pi * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   *         "-0.25 * pi * pi * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"]
   */
  virtual void
  jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    const DomainFieldType pre = -0.25 * M_PIl * M_PIl * M_PIl;
    const DomainFieldType x_arg = M_PI_2l * xx[0];
    const DomainFieldType y_arg = M_PI_2l * xx[1];
    ret[0][0] = pre * sin(x_arg) * cos(y_arg);
    ret[0][1] = pre * cos(x_arg) * sin(y_arg);
  } // ... jacobian(...)

private:
  const size_t order_;
  const std::string name_;
}; // class Testcase1Force


template <class E, class D, size_t d, class R, size_t r, size_t rC = 1>
class Testcase1ExactSolution : public LocalizableFunctionInterface<E, D, d, R, r, rC>
{
  Testcase1ExactSolution()
  {
    static_assert(AlwaysFalse<E>::value, "Not available for these dimensions!");
  }
};


template <class EntityImp, class DomainFieldImp, class RangeFieldImp>
class Testcase1ExactSolution<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1>
    : public GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1>
{
  typedef Testcase1ExactSolution<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1> ThisType;
  typedef GlobalFunctionInterface<EntityImp, DomainFieldImp, 2, RangeFieldImp, 1, 1> BaseType;

public:
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::JacobianRangeType;

  static const bool available = true;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.testcase1.exactsolution";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config;
    config["integration_order"] = "3";
    config["name"] = static_id();
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    // get correct config
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // create
    return Common::make_unique<ThisType>(cfg.get("integration_order", default_cfg.get<size_t>("integration_order")),
                                         cfg.get("name", default_cfg.get<std::string>("name")));
  } // ... create(...)

  Testcase1ExactSolution(const size_t ord = default_config().template get<size_t>("integration_order"),
                         const std::string nm = static_id())
    : order_(ord)
    , name_(nm)
  {
  }

  Testcase1ExactSolution(const ThisType& /*other*/) = default;

  ThisType& operator=(const ThisType& /*other*/) = delete;

  virtual std::string type() const override final
  {
    return BaseType::static_id() + ".ESV2007.testcase1.exactsolution";
  }

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    return order_;
  }

  /**
   * \brief "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   */
  virtual void evaluate(const DomainType& xx, RangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    ret[0] = cos(M_PI_2l * xx[0]) * cos(M_PI_2l * xx[1]);
  }

  /**
   * \brief ["-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])"
   *         "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"]
   */
  virtual void
  jacobian(const DomainType& xx, JacobianRangeType& ret, const Common::Parameter& /*mu*/ = {}) const override final
  {
    const DomainFieldType pre = -0.5 * M_PIl;
    const DomainFieldType x_arg = M_PI_2l * xx[0];
    const DomainFieldType y_arg = M_PI_2l * xx[1];
    ret[0][0] = pre * sin(x_arg) * cos(y_arg);
    ret[0][1] = pre * cos(x_arg) * sin(y_arg);
  } // ... jacobian(...)

private:
  const size_t order_;
  const std::string name_;
}; // class Testcase1ExactSolution


template <class DiffusionFactorType, class DiffusionTensorType = void>
class CutoffFunction;

template <class DiffusionType>
class CutoffFunction<DiffusionType, void> : public LocalizableFunctionInterface<typename DiffusionType::EntityType,
                                                                                typename DiffusionType::DomainFieldType,
                                                                                DiffusionType::dimDomain,
                                                                                typename DiffusionType::RangeFieldType,
                                                                                1,
                                                                                1>
{
  static_assert(is_localizable_function<DiffusionType>::value,
                "DiffusionType has to be tagged as a LocalizableFunction!");
  typedef typename DiffusionType::EntityType E_;
  typedef typename DiffusionType::DomainFieldType D_;
  static const size_t d_ = DiffusionType::dimDomain;
  typedef typename DiffusionType::RangeFieldType R_;
  typedef LocalizableFunctionInterface<E_, D_, d_, R_, 1> BaseType;
  typedef CutoffFunction<DiffusionType> ThisType;

  class Localfunction : public LocalfunctionInterface<E_, D_, d_, R_, 1, 1>
  {
    typedef LocalfunctionInterface<E_, D_, d_, R_, 1, 1> BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const size_t dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const size_t dimRange = BaseType::dimRange;
    static const size_t dimRangeCols = BaseType::dimRangeCols;
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
    Localfunction(const EntityType& ent, const DiffusionType& diffusion, const RangeFieldType poincare_constant)
      : BaseType(ent)
      , value_(0)
    {
      const RangeFieldType min_eigen_value =
          Compute<DiffusionType, DiffusionType::dimRange, DiffusionType::dimRangeCols>::min_eigenvalue_of(diffusion,
                                                                                                          ent);
      assert(min_eigen_value > 0.0);
      const DomainFieldType hh = compute_diameter_of_(ent);
      value_ = (poincare_constant * hh * hh) / min_eigen_value;
    }

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
    {
      return 0;
    }

    virtual void evaluate(const DomainType& DXTC_DEBUG_ONLY(xx),
                          RangeType& ret,
                          const Common::Parameter& /*mu*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret[0] = value_;
    }

    virtual void jacobian(const DomainType& DXTC_DEBUG_ONLY(xx),
                          JacobianRangeType& ret,
                          const Common::Parameter& /*mu*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    static DomainFieldType compute_diameter_of_(const EntityType& ent)
    {
      DomainFieldType ret(0);
      for (auto cc : Common::value_range(ent.subEntities(dimDomain))) {
        const auto vertex = ent.template subEntity<dimDomain>(cc).geometry().center();
        for (auto dd : Common::value_range(cc + 1, ent.subEntities(dimDomain))) {
          const auto other_vertex = ent.template subEntity<dimDomain>(dd).geometry().center();
          const auto diff = vertex - other_vertex;
          ret = std::max(ret, diff.two_norm());
        }
      }
      return ret;
    } // ... compute_diameter_of_(...)

    RangeFieldType value_;
  }; // class Localfunction

public:
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::LocalfunctionType LocalfunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.cutoff";
  }

  CutoffFunction(const DiffusionType& diffusion,
                 const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl),
                 const std::string nm = static_id())
    : diffusion_(diffusion)
    , poincare_constant_(poincare_constant)
    , name_(nm)
  {
  }

  CutoffFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return std::unique_ptr<Localfunction>(new Localfunction(entity, diffusion_, poincare_constant_));
  }

private:
  const DiffusionType& diffusion_;
  const RangeFieldType poincare_constant_;
  std::string name_;
}; // class Cutoff

template <class DiffusionFactorType, class DiffusionTensorType>
class CutoffFunction : public LocalizableFunctionInterface<typename DiffusionFactorType::EntityType,
                                                           typename DiffusionFactorType::DomainFieldType,
                                                           DiffusionFactorType::dimDomain,
                                                           typename DiffusionFactorType::RangeFieldType,
                                                           1,
                                                           1>
{
  static_assert(is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be tagged as a LocalizableFunction!");
  static_assert(is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be tagged as a LocalizableFunction!");
  typedef typename DiffusionFactorType::EntityType E_;
  typedef typename DiffusionFactorType::DomainFieldType D_;
  static const size_t d_ = DiffusionFactorType::dimDomain;
  typedef typename DiffusionFactorType::RangeFieldType R_;
  typedef LocalizableFunctionInterface<E_, D_, d_, R_, 1> BaseType;
  typedef CutoffFunction<DiffusionFactorType, DiffusionTensorType> ThisType;

  static_assert(DiffusionFactorType::dimRange == 1, "The diffusion factor has to be scalar!");
  static_assert(DiffusionFactorType::dimRangeCols == 1, "The diffusion factor has to be scalar!");

  static_assert(std::is_same<typename DiffusionTensorType::EntityType, E_>::value, "Types do not match!");
  static_assert(std::is_same<typename DiffusionTensorType::DomainFieldType, D_>::value, "Types do not match!");
  static_assert(DiffusionTensorType::dimDomain == d_, "Dimensions do not match!");
  static_assert(std::is_same<typename DiffusionTensorType::RangeFieldType, R_>::value, "Types do not match!");

  static_assert(DiffusionTensorType::dimRange == d_, "The diffusion tensor has to be a matrix!");
  static_assert(DiffusionTensorType::dimRangeCols == d_, "The diffusion tensor has to be a matrix!");

  class Localfunction : public LocalfunctionInterface<E_, D_, d_, R_, 1, 1>
  {
    typedef LocalfunctionInterface<E_, D_, d_, R_, 1, 1> BaseType;

  public:
    typedef typename BaseType::EntityType EntityType;

    typedef typename BaseType::DomainFieldType DomainFieldType;
    static const size_t dimDomain = BaseType::dimDomain;
    typedef typename BaseType::DomainType DomainType;

    typedef typename BaseType::RangeFieldType RangeFieldType;
    static const size_t dimRange = BaseType::dimRange;
    static const size_t dimRangeCols = BaseType::dimRangeCols;
    typedef typename BaseType::RangeType RangeType;

    typedef typename BaseType::JacobianRangeType JacobianRangeType;

  private:
    template <class DF, size_t r, size_t rR>
    struct ComputeDiffusionFactor
    {
      static_assert(AlwaysFalse<DF>::value, "Not implemented for these dimensions!");
    };

    template <class DF>
    struct ComputeDiffusionFactor<DF, 1, 1>
    {
      /**
       * We try to find the minimum of a polynomial of given order by evaluating it at the points of a quadrature that
       * would integrate this polynomial exactly.
       * \todo These are just some heuristics and should be replaced by something proper.
       */
      static RangeFieldType min_of(const DF& diffusion_factor, const EntityType& ent)
      {
        typename DF::RangeType tmp_value(0);
        RangeFieldType minimum = std::numeric_limits<RangeFieldType>::max();
        const auto local_diffusion_factor = diffusion_factor.local_function(ent);
        const size_t ord = local_diffusion_factor->order();
        const auto& quadrature =
            QuadratureRules<DomainFieldType, dimDomain>::rule(ent.type(), boost::numeric_cast<int>(ord));
        const auto quad_point_it_end = quadrature.end();
        for (auto quad_point_it = quadrature.begin(); quad_point_it != quad_point_it_end; ++quad_point_it) {
          local_diffusion_factor->evaluate(quad_point_it->position(), tmp_value);
          minimum = std::min(minimum, tmp_value[0]);
        }
        return minimum;
      } // ... min_of(...)
    }; // class ComputeDiffusionFactor< ..., 1, 1 >

    template <class DT, size_t r, size_t rR>
    struct ComputeDiffusionTensor
    {
      static_assert(AlwaysFalse<DT>::value, "Not implemented for these dimensions!");
    };

    template <class DT, size_t d>
    struct ComputeDiffusionTensor<DT, d, d>
    {
      static RangeFieldType min_eigenvalue_of(const DT& diffusion_tensor, const EntityType& ent)
      {
#if !HAVE_DUNE_XT_LA
        static_assert(AlwaysFalse<DT>::value, "You are missing dune-xt-la!");
#else
#if !HAVE_EIGEN
        static_assert(AlwaysFalse<DT>::value, "You are missing eigen!");
#else
        const auto local_diffusion_tensor = diffusion_tensor.local_function(ent);
        assert(local_diffusion_tensor->order() == 0);
        const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(ent.type());
        const XT::LA::EigenDenseMatrix<RangeFieldType> tensor =
            local_diffusion_tensor->evaluate(reference_element.position(0, 0));
        ::Eigen::EigenSolver<typename LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver(tensor.backend());

        assert(eigen_solver.info() == ::Eigen::Success);
        const auto eigenvalues = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of

        RangeFieldType min_ev = std::numeric_limits<RangeFieldType>::max();
        for (size_t ii = 0; ii < boost::numeric_cast<size_t>(eigenvalues.size()); ++ii) {
          // assert this is real
          assert(std::abs(eigenvalues[ii].imag()) < 1e-15);
          // assert that this eigenvalue is positive
          const RangeFieldType eigenvalue = eigenvalues[ii].real();
          assert(eigenvalue > 1e-15);
          min_ev = std::min(min_ev, eigenvalue);
        }
        return min_ev;
#endif // HAVE_EIGEN
#endif // HAVE_DUNE_XT_LA
      } // ... min_eigenvalue_of_(...)
    }; // class Compute< ..., d, d >

  public:
    Localfunction(const EntityType& ent,
                  const DiffusionFactorType& diffusion_factor,
                  const DiffusionTensorType& diffusion_tensor,
                  const RangeFieldType poincare_constant)
      : BaseType(ent)
      , value_(0)
    {
      const RangeFieldType min_diffusion_factor =
          ComputeDiffusionFactor<DiffusionFactorType,
                                 DiffusionFactorType::dimRange,
                                 DiffusionFactorType::dimRangeCols>::min_of(diffusion_factor, ent);
      const RangeFieldType min_eigen_value_diffusion_tensor =
          ComputeDiffusionTensor<DiffusionTensorType,
                                 DiffusionTensorType::dimRange,
                                 DiffusionTensorType::dimRangeCols>::min_eigenvalue_of(diffusion_tensor, ent);
      assert(min_diffusion_factor > RangeFieldType(0));
      assert(min_eigen_value_diffusion_tensor > RangeFieldType(0));
      const DomainFieldType hh = compute_diameter_of_(ent);
      value_ = (poincare_constant * hh * hh) / (min_diffusion_factor * min_eigen_value_diffusion_tensor);
    } // Localfunction(...)

    Localfunction(const Localfunction& /*other*/) = delete;

    Localfunction& operator=(const Localfunction& /*other*/) = delete;

    virtual size_t order(const XT::Common::Parameter& /*mu*/ = {}) const override final
    {
      return 0;
    }

    virtual void evaluate(const DomainType& DXTC_DEBUG_ONLY(xx),
                          RangeType& ret,
                          const Common::Parameter& /*mu*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret[0] = value_;
    }

    virtual void jacobian(const DomainType& DXTC_DEBUG_ONLY(xx),
                          JacobianRangeType& ret,
                          const Common::Parameter& /*mu*/ = {}) const override final
    {
      assert(this->is_a_valid_point(xx));
      ret *= RangeFieldType(0);
    }

  private:
    static DomainFieldType compute_diameter_of_(const EntityType& ent)
    {
      DomainFieldType ret(0);
      for (auto cc : Common::value_range(ent.subEntities(dimDomain))) {
        const auto vertex = ent.template subEntity<dimDomain>(cc).geometry().center();
        for (auto dd : Common::value_range(cc + 1, ent.subEntities(dimDomain))) {
          const auto other_vertex = ent.template subEntity<dimDomain>(dd).geometry().center();
          const auto diff = vertex - other_vertex;
          ret = std::max(ret, diff.two_norm());
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
  static const size_t dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType RangeType;

  static std::string static_id()
  {
    return BaseType::static_id() + ".ESV2007.cutoff";
  }

  CutoffFunction(const DiffusionFactorType& diffusion_factor,
                 const DiffusionTensorType& diffusion_tensor,
                 const RangeFieldType poincare_constant = 1.0 / (M_PIl * M_PIl),
                 const std::string nm = static_id())
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , poincare_constant_(poincare_constant)
    , name_(nm)
  {
  }

  CutoffFunction(const ThisType& other) = default;

  ThisType& operator=(const ThisType& other) = delete;

  virtual std::string name() const override final
  {
    return name_;
  }

  virtual std::unique_ptr<LocalfunctionType> local_function(const EntityType& entity) const override final
  {
    return std::unique_ptr<Localfunction>(
        new Localfunction(entity, diffusion_factor_, diffusion_tensor_, poincare_constant_));
  }

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const RangeFieldType poincare_constant_;
  std::string name_;
}; // class Cutoff

} // namespace ESV2007
} // namespace Functions
} // namespace XT
} // namespace Dune

#endif // DUNE_XT_FUNCTIONS_ESV2007_HH
