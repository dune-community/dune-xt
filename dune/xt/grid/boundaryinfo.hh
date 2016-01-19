// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff
// The copyright lies with the authors of this file (see below).
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2012 - 2015)
//   Rene Milk       (2012 - 2015)
//   Sven Kaulmann   (2013 - 2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_XT_GRID_BOUNDARYINFO_HH
#define DUNE_XT_GRID_BOUNDARYINFO_HH

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#if HAVE_DUNE_PDELAB
#include <dune/typetree/leafnode.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/fvector.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/timedlogging.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace internal {

inline static std::string boundary_info_static_id()
{
  return "stuff.grid.boundaryinfo";
}

} // namespace internal

template <class IntersectionImp>
class BoundaryInfoInterface
#if HAVE_DUNE_PDELAB
    : public TypeTree::LeafNode
#endif
{
#if HAVE_DUNE_PDELAB
  typedef PDELab::IntersectionGeometry<IntersectionImp> IntersectionGeometryType;
  typedef FieldVector<typename IntersectionGeometryType::ctype, IntersectionGeometryType::dimension - 1>
      Codim1DomainType;
#endif // HAVE_DUNE_PDELAB
public:
  typedef IntersectionImp IntersectionType;
  typedef typename IntersectionType::ctype DomainFieldType;
  static const size_t dimDomain = IntersectionType::dimension;
  static const size_t dimWorld  = IntersectionType::dimensionworld;
  typedef Common::FieldVector<DomainFieldType, dimDomain> DomainType;
  typedef Common::FieldVector<DomainFieldType, dimWorld> WorldType;

  static const std::string static_id()
  {
    return internal::boundary_info_static_id();
  }

  /** \note throw spec here and in derived classes are only there to work around icc-mic bug
   **/
  virtual ~BoundaryInfoInterface() throw()
  {
  }

  virtual bool has_dirichlet() const
  {
    return true;
  }

  virtual bool has_neumann() const
  {
    return true;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const = 0;

  virtual bool neumann(const IntersectionType& intersection) const = 0;

#if HAVE_DUNE_PDELAB
  //! default implementation for signature used in PDELab
  virtual bool isDirichlet(const IntersectionGeometryType& intersection_geometry,
                           const Codim1DomainType& /*coord*/) const
  {
    return dirichlet(intersection_geometry.intersection());
  }

  //! default implementation for signature used in PDELab
  virtual bool isNeumann(const IntersectionGeometryType& intersection_geometry, const Codim1DomainType& /*coord*/) const
  {
    return neumann(intersection_geometry.intersection());
  }
#endif // HAVE_DUNE_PDELAB
}; // class BoundaryInfoInterface

namespace BoundaryInfoConfigs {

class AllDirichlet
{
public:
  static std::string static_id()
  {
    return internal::boundary_info_static_id() + ".alldirichlet";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    if (sub_name.empty())
      return Common::Configuration("type", static_id());
    else
      return Common::Configuration(sub_name + ".type", static_id());
  }
}; // class AllDirichlet

class AllNeumann
{
public:
  static std::string static_id()
  {
    return internal::boundary_info_static_id() + ".allneumann";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    if (sub_name.empty())
      return Common::Configuration("type", static_id());
    else
      return Common::Configuration(sub_name + ".type", static_id());
  }
}; // class AllNeumann

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

class IdBased
{
public:
  static const std::string static_id()
  {
    return internal::boundary_info_static_id() + ".idbased";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config("type", static_id());
    config["default_to_dirichlet"] = "true";
    config["neumann"] = "[4]";
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)
}; // class IdBased

#endif // #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

class NormalBased
{
public:
  static const std::string static_id()
  {
    return internal::boundary_info_static_id() + ".normalbased";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config("type", static_id());
    config["default"]           = "dirichlet";
    config["compare_tolerance"] = "1e-10";
    config["neumann.0"]         = "[1.0 0.0]";
    config["dirichlet.0"] = "[0.0 1.0]";
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)
}; // class NormalBased

class FunctionBased
{
public:
  static const std::string static_id()
  {
    return internal::boundary_info_static_id() + ".functionbased";
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    Common::Configuration config("type", static_id());
    config["default"]                       = "dirichlet";
    config["dirichlet.function.type"]       = "stuff.function.constant";
    config["dirichlet.function.expression"] = "1";
    config["dirichlet.value_range.1"]       = "[0 0.1]";
    config["neumann.function.type"]         = "stuff.function.constant";
    config["neumann.function.expression"]   = "1";
    config["neumann.value_range.1"] = "[0 0.1]";
    if (sub_name.empty())
      return config;
    else {
      Common::Configuration tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)
}; // class FunctionBased

} // namespace BoundaryInfoConfigs

namespace BoundaryInfos {

template <class IntersectionImp>
class AllDirichlet : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef AllDirichlet<IntersectionImp> ThisType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string static_id()
  {
    return BoundaryInfoConfigs::AllDirichlet::static_id();
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return BoundaryInfoConfigs::AllDirichlet::default_config(sub_name);
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration /*config*/ = default_config(),
                                          const std::string /*sub_name*/ = static_id())
  {
    return Common::make_unique<ThisType>();
  }

  AllDirichlet()
  {
  }

  virtual ~AllDirichlet() throw()
  {
  }

  virtual bool has_dirichlet() const override final
  {
    return true;
  }

  virtual bool has_neumann() const override final
  {
    return false;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }

  virtual bool neumann(const IntersectionType& /*intersection*/) const override final
  {
    return false;
  }
}; // class AllDirichlet

template <class IntersectionImp>
class AllNeumann : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef AllNeumann<IntersectionImp> ThisType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string static_id()
  {
    return BoundaryInfoConfigs::AllNeumann::static_id();
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return BoundaryInfoConfigs::AllNeumann::default_config(sub_name);
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration /*config*/ = default_config(),
                                          const std::string /*sub_name*/ = static_id())
  {
    return Common::make_unique<ThisType>();
  }

  AllNeumann()
  {
  }

  virtual ~AllNeumann() throw(){};

  virtual bool has_dirichlet() const override final
  {
    return false;
  }

  virtual bool has_neumann() const override final
  {
    return true;
  }

  virtual bool dirichlet(const IntersectionType& /*intersection*/) const override final
  {
    return false;
  }

  virtual bool neumann(const IntersectionType& intersection) const override final
  {
    return intersection.boundary();
  }
}; // class AllNeumann

#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

template <class IntersectionImp>
class IdBased : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef IdBased<IntersectionImp> ThisType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string static_id()
  {
    return BoundaryInfoConfigs::IdBased::static_id();
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return BoundaryInfoConfigs::IdBased::default_config(sub_name);
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    std::map<std::string, std::set<int>> id_to_type_map;
    for (const std::string& type : {"dirichlet", "neumann"})
      if (cfg.has_key(type)) {
        const auto ids = cfg.get<std::vector<int>>(type);
        id_to_type_map.insert(std::make_pair(type, std::set<int>(ids.begin(), ids.end())));
      }
    return Common::make_unique<ThisType>(
        id_to_type_map, cfg.get("default_to_dirichlet", default_config().template get<bool>("default_to_dirichlet")));
  }

  IdBased(const std::map<std::string, std::set<int>> id_to_type_map_in, const bool dirichlet_default = true)
    : id_to_type_map_(id_to_type_map_in)
    , dirichletDefault_(dirichlet_default)
    , hasDirichlet_(id_to_type_map_.find("dirichlet") != id_to_type_map_.end() || dirichletDefault_)
    , hasNeumann_(id_to_type_map_.find("neumann") != id_to_type_map_.end() || !dirichletDefault_)
  {
  }

  virtual ~IdBased() = default;

  const std::map<std::string, std::set<int>>& id_to_type_map() const
  {
    return id_to_type_map_;
  }

  virtual bool has_dirichlet() const override final
  {
    return hasDirichlet_;
  }

  virtual bool has_neumann() const override final
  {
    return hasNeumann_;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const override final
  {
    if (hasDirichlet_ && intersection.boundary()) {
#include <dune/stuff/common/disable_warnings.hh>
      const int boundaryId = intersection.boundaryId();
#include <dune/stuff/common/reenable_warnings.hh>
      // get set of dirichlet ids
      const auto result = id_to_type_map_.find("dirichlet");
      if (result == id_to_type_map_.end())
        return dirichletDefault_;
      else {
        const auto& idSet = result->second;
        return (idSet.count(boundaryId) > 0) || dirichletDefault_;
      }
    } else
      return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType& intersection) const override final
  {
    if (hasNeumann_ && intersection.boundary()) {
// get boundary id
#include <dune/stuff/common/disable_warnings.hh>
      const int boundaryId = intersection.boundaryId();
#include <dune/stuff/common/reenable_warnings.hh>
      // get set of neumann ids
      const auto result = id_to_type_map_.find("neumann");
      if (result == id_to_type_map_.end())
        return !dirichletDefault_;
      else {
        const auto& idSet = result->second;
        return (idSet.count(boundaryId) > 0) || !dirichletDefault_;
      }
    } else
      return false;
  } // bool neumann(const IntersectionType& intersection) const

private:
  std::map<std::string, std::set<int>> id_to_type_map_;
  bool dirichletDefault_;
  bool hasDirichlet_;
  bool hasNeumann_;
}; // class IdBased

#endif // #if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS

template <class IntersectionImp>
class NormalBased : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef NormalBased<IntersectionImp> ThisType;

public:
  using typename BaseType::IntersectionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::WorldType;
  using BaseType::dimDomain;
  using BaseType::dimWorld;

  static const std::string static_id()
  {
    return BoundaryInfoConfigs::NormalBased::static_id();
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return BoundaryInfoConfigs::NormalBased::default_config(sub_name);
  }

  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // get default
    const std::string default_type = cfg.get("default", default_cfg.get<std::string>("default"));
    if (default_type != "dirichlet" && default_type != "neumann")
      DUNE_THROW(Exceptions::configuration_error, "Wrong default '" << default_type << "' given!");
    const bool default_to_dirichlet = default_type == "dirichlet";
    // get tolerance
    const DomainFieldType tol = cfg.get("compare_tolerance", default_cfg.get<DomainFieldType>("compare_tolerance"));
    // get dirichlet and neumann
    std::vector<WorldType> dirichlets = getVectors(cfg, "dirichlet");
    std::vector<WorldType> neumanns   = getVectors(cfg, "neumann");
    // return
    return Common::make_unique<ThisType>(default_to_dirichlet, dirichlets, neumanns, tol);
  } // ... create(...)

  NormalBased(const bool default_to_dirichlet = true,
              const std::vector<WorldType> dirichlet_normals = std::vector<WorldType>(),
              const std::vector<WorldType> neumann_normals   = std::vector<WorldType>(),
              const DomainFieldType tol = 1e-10)
    : default_to_dirichlet_(default_to_dirichlet)
    , dirichlet_normals_(dirichlet_normals)
    , neumann_normals_(neumann_normals)
    , tol_(tol)
  {
    // normalize
    for (auto& normal : dirichlet_normals_)
      normal /= normal.two_norm();
    for (auto& normal : neumann_normals_)
      normal /= normal.two_norm();
    // sanity check
    for (auto& dirichletNormal : dirichlet_normals_) {
      if (contains(dirichletNormal, neumann_normals_))
        DUNE_THROW(Exceptions::wrong_input_given, "Given normals are too close for given tolerance '" << tol << "'!");
    }
  } // NormalBased(...)

  virtual ~NormalBased() throw()
  {
  }

  virtual bool has_dirichlet() const override final
  {
    return default_to_dirichlet_ || (dirichlet_normals_.size() > 0);
  }

  virtual bool has_neumann() const override final
  {
    return !default_to_dirichlet_ || (neumann_normals_.size() > 0);
  }

  virtual bool dirichlet(const IntersectionType& intersection) const override final
  {
    if (intersection.boundary()) {
      const WorldType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, dirichlet_normals_))
        return true;
      else if (contains(outerNormal, neumann_normals_))
        return false;
      else
        return default_to_dirichlet_;
    }
    return false;
  } // ... dirichlet(...)

  virtual bool neumann(const IntersectionType& intersection) const override final
  {
    if (intersection.boundary()) {
      const WorldType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, neumann_normals_))
        return true;
      else if (contains(outerNormal, dirichlet_normals_))
        return false;
      else
        return !default_to_dirichlet_;
    }
    return false;
  } // ... neumann(...)

private:
  static std::vector<WorldType> getVectors(const Common::Configuration& config, const std::string key)
  {
    std::vector<WorldType> ret;
    if (config.has_sub(key)) {
      bool found     = true;
      size_t counter = 0;
      while (found) {
        const std::string localKey = key + "." + Dune::Stuff::Common::toString(counter);
        if (config.has_key(localKey))
          ret.push_back(config.get<WorldType>(localKey, dimWorld));
        else
          found = false;
        ++counter;
      }
    } else if (config.has_key(key))
      ret.push_back(config.get<WorldType>(key, dimWorld));
    return ret;
  } // ... getVectors(...)

  bool contains(const WorldType& normal, const std::vector<WorldType>& vectors) const
  {
    for (auto& vector : vectors)
      if (Dune::Stuff::Common::FloatCmp::eq(normal, vector, tol_))
        return true;
    return false;
  }

  const bool default_to_dirichlet_;
  std::vector<WorldType> dirichlet_normals_;
  std::vector<WorldType> neumann_normals_;
  const DomainFieldType tol_;
}; // class NormalBased

template <class IntersectionImp>
class FunctionBased : public BoundaryInfoInterface<IntersectionImp>
{
  typedef BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef FunctionBased<IntersectionImp> ThisType;

public:
  using typename BaseType::IntersectionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using BaseType::dimDomain;

private:
  typedef typename IntersectionType::Entity EntityType;
  typedef LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, double, 1> FunctionType;

public:
  static const std::string static_id()
  {
    return BoundaryInfoConfigs::FunctionBased::static_id();
  }

  static Common::Configuration default_config(const std::string sub_name = "")
  {
    return BoundaryInfoConfigs::FunctionBased::default_config(sub_name);
  }

private:
  static std::vector<Common::FieldVector<double, 2>>
  get_value_range(const Common::Configuration& cfg, const Common::Configuration& default_cfg, const std::string& id)
  {
    auto logger = Common::TimedLogger().get("stuff.grid.boundaryinfo.get_value_range");
    std::vector<Common::FieldVector<double, 2>> value_range;
    if (cfg.has_key(id) && !cfg.has_sub(id))
      value_range.push_back(cfg.get<Common::FieldVector<double, 2>>(id));
    else {
      Common::Configuration sub_cfg = cfg.has_sub(id) ? cfg.sub(id) : default_cfg.sub(id);
      size_t counter = 0;
      while (sub_cfg.has_key(Common::toString(counter))) {
        value_range.push_back(sub_cfg.get<Common::FieldVector<double, 2>>(Common::toString(counter)));
        ++counter;
      }
    }
    logger.debug() << "found " << value_range.size() << " value ranges for id '" << id << "'" << std::endl;
    return value_range;
  } // ... get_value_range(...)

public:
  static std::unique_ptr<ThisType> create(const Common::Configuration config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::Configuration cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::Configuration default_cfg = default_config();
    // get default
    const std::string default_type = cfg.get("default", default_cfg.get<std::string>("default"));
    if (default_type != "dirichlet" && default_type != "neumann")
      DUNE_THROW(Exceptions::configuration_error, "Wrong default '" << default_type << "' given!");
    const bool default_to_dirichlet = default_type == "dirichlet";
    // create
    typedef FunctionsProvider<EntityType, DomainFieldType, dimDomain, double, 1> FunctionProviderType;
    return Common::make_unique<ThisType>(
        default_to_dirichlet,
        FunctionProviderType::create(cfg.has_sub("dirichlet.function") ? cfg.sub("dirichlet.function")
                                                                       : default_cfg.sub("dirichlet.function")),
        FunctionProviderType::create(cfg.has_sub("dirichlet.function") ? cfg.sub("dirichlet.function")
                                                                       : default_cfg.sub("dirichlet.function")),
        get_value_range(cfg, default_cfg, "dirichlet.value_range"),
        get_value_range(cfg, default_cfg, "neumann.value_range"));
  } // ... create(...)

  FunctionBased(const bool default_is_dirichlet, std::shared_ptr<const FunctionType> dirichlet_function,
                std::shared_ptr<const FunctionType> neumann_function,
                const std::vector<Common::FieldVector<double, 2>>& dirichlet_value_range,
                const std::vector<Common::FieldVector<double, 2>>& neumann_value_range)
    : default_is_dirichlet_(default_is_dirichlet)
    , dirichlet_function_(dirichlet_function)
    , neumann_function_(neumann_function)
    , dirichlet_value_range_(dirichlet_value_range)
    , neumann_value_range_(neumann_value_range)
  {
  }

  virtual ~FunctionBased() throw()
  {
  }

  virtual bool dirichlet(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary())
      return false;
    const auto entity_ptr = intersection.inside();
    const auto& entity    = *entity_ptr;
    const auto xx_global  = intersection.geometry().center();
    const auto xx_entity  = entity.geometry().local(xx_global);
    auto logger = Common::TimedLogger().get("stuff.grid.boundaryinfo.dirichlet");
    logger.debug() << "xx_global = " << xx_global << std::endl;
    logger.debug() << "xx_entity = " << xx_entity << std::endl;
    const auto local_dirichlet_function = dirichlet_function_->local_function(entity);
    const auto local_neumann_function = neumann_function_->local_function(entity);
    for (const auto& value_range : dirichlet_value_range_)
      if (matches(*local_dirichlet_function, value_range, xx_entity))
        return true;
    for (const auto& value_range : neumann_value_range_)
      if (matches(*local_neumann_function, value_range, xx_entity))
        return false;
    return default_is_dirichlet_;
  } // ... dirichlet(...)

  virtual bool neumann(const IntersectionType& intersection) const override final
  {
    if (!intersection.boundary())
      return false;
    const auto entity_ptr               = intersection.inside();
    const auto& entity                  = *entity_ptr;
    const auto xx_global                = intersection.geometry().center();
    const auto xx_entity                = entity.geometry().local(xx_global);
    const auto local_dirichlet_function = dirichlet_function_->local_function(entity);
    const auto local_neumann_function = neumann_function_->local_function(entity);
    for (const auto& value_range : neumann_value_range_)
      if (matches(*local_neumann_function, value_range, xx_entity))
        return true;
    for (const auto& value_range : dirichlet_value_range_)
      if (matches(*local_dirichlet_function, value_range, xx_entity))
        return false;
    return !default_is_dirichlet_;
  } // ... neumann(...)

private:
  bool matches(const typename FunctionType::LocalfunctionType& local_function,
               const Common::FieldVector<double, 2>& value_range, const DomainType& xx) const
  {
    auto logger      = Common::TimedLogger().get("stuff.grid.boundaryinfo.matches");
    const auto value = local_function.evaluate(xx)[0];
    bool does = Common::FloatCmp::ge(value, value_range[0]) && Common::FloatCmp::le(value, value_range[1]);
    logger.debug() << "value = " << value << std::endl;
    logger.debug() << "match: " << does << std::endl;
    return does;
  }

  const bool default_is_dirichlet_;
  const std::shared_ptr<const FunctionType> dirichlet_function_;
  const std::shared_ptr<const FunctionType> neumann_function_;
  const std::vector<Common::FieldVector<double, 2>> dirichlet_value_range_;
  const std::vector<Common::FieldVector<double, 2>> neumann_value_range_;
}; // class FunctionBased

} // namespace BoundaryInfos

template <class I>
class BoundaryInfoProvider
{
public:
  typedef BoundaryInfoInterface<I> InterfaceType;

  static std::vector<std::string> available()
  {
    using namespace Stuff::Grid::BoundaryInfos;
    return
    {
      BoundaryInfos::AllDirichlet<I>::static_id(), BoundaryInfos::AllNeumann<I>::static_id()
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
                                                       ,
          BoundaryInfos::IdBased<I>::static_id()
#endif
              ,
          BoundaryInfos::NormalBased<I>::static_id(), BoundaryInfos::FunctionBased<I>::static_id()
    };
  } // ... available(...)

  static Common::Configuration default_config(const std::string type, const std::string subname = "")
  {
    using namespace Stuff::Grid::BoundaryInfos;
    if (type == BoundaryInfos::AllDirichlet<I>::static_id())
      return BoundaryInfos::AllDirichlet<I>::default_config(subname);
    else if (type == BoundaryInfos::AllNeumann<I>::static_id())
      return BoundaryInfos::AllNeumann<I>::default_config(subname);
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    else if (type == BoundaryInfos::IdBased<I>::static_id())
      return BoundaryInfos::IdBased<I>::default_config(subname);
#endif
    else if (type == BoundaryInfos::NormalBased<I>::static_id())
      return BoundaryInfos::NormalBased<I>::default_config(subname);
    else if (type == BoundaryInfos::FunctionBased<I>::static_id())
      return BoundaryInfos::FunctionBased<I>::default_config(subname);
    else
      DUNE_THROW(Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const Common::Configuration& config)
  {
    return create(config.get<std::string>("type"), config);
  }

  static std::unique_ptr<InterfaceType> create(const std::string& type = available()[0],
                                               const Common::Configuration config = default_config(available()[0]))
  {
    using namespace Stuff::Grid::BoundaryInfos;
    if (type == BoundaryInfos::AllDirichlet<I>::static_id())
      return BoundaryInfos::AllDirichlet<I>::create(config);
    else if (type == BoundaryInfos::AllNeumann<I>::static_id())
      return BoundaryInfos::AllNeumann<I>::create(config);
#if DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS
    else if (type == BoundaryInfos::IdBased<I>::static_id())
      return BoundaryInfos::IdBased<I>::create(config);
#endif
    else if (type == BoundaryInfos::NormalBased<I>::static_id())
      return BoundaryInfos::NormalBased<I>::create(config);
    else if (type == BoundaryInfos::FunctionBased<I>::static_id())
      return BoundaryInfos::FunctionBased<I>::create(config);
    else
      DUNE_THROW(Exceptions::wrong_input_given,
                 "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
}; // class BoundaryInfoProvider

} // namespace Grid

} // namespace Stuff
} // namespace Dune

#endif // DUNE_XT_GRID_BOUNDARYINFO_HH
