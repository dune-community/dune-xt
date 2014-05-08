// This file is part of the dune-stuff project:
//   https://users.dune-project.org/projects/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Sven Kaulmann

#ifndef DUNE_STUFF_GRID_BOUNDARYINFO_HH
#define DUNE_STUFF_GRID_BOUNDARYINFO_HH

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/configtree.hh>
#include <dune/stuff/common/memory.hh>

#if HAVE_DUNE_PDELAB
#include <dune/typetree/leafnode.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#endif

namespace Dune {
namespace Stuff {
namespace Grid {


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
#endif
public:
  typedef IntersectionImp IntersectionType;

  static const std::string static_id()
  {
    return "stuff.grid.boundaryinfo";
  }

  virtual ~BoundaryInfoInterface()
  {
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


namespace BoundaryInfos {


template <class IntersectionImp>
class AllDirichlet : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef AllDirichlet<IntersectionImp> ThisType;

public:
  using typename BaseType::IntersectionType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".alldirichlet";
  }

  static Common::ConfigTree default_config(const std::string sub_name = "")
  {
    if (sub_name.empty())
      return Common::ConfigTree("type", static_id());
    else
      return Common::ConfigTree(sub_name + ".type", static_id());
  }

  static std::unique_ptr<ThisType> create(const Common::ConfigTree /*config*/ = default_config(),
                                          const std::string /*sub_name*/ = static_id())
  {
    return Common::make_unique<ThisType>();
  }

  AllDirichlet()
  {
  }

  virtual ~AllDirichlet()
  {
  }

  virtual bool dirichlet(const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return intersection.boundary();
  }

  virtual bool neumann(const IntersectionType& /*intersection*/) const DS_OVERRIDE DS_FINAL
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
  using typename BaseType::IntersectionType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".allneumann";
  }

  static Common::ConfigTree default_config(const std::string sub_name = "")
  {
    if (sub_name.empty())
      return Common::ConfigTree("type", static_id());
    else
      return Common::ConfigTree(sub_name + ".type", static_id());
  }

  static std::unique_ptr<ThisType> create(const Common::ConfigTree /*config*/ = default_config(),
                                          const std::string /*sub_name*/ = static_id())
  {
    return Common::make_unique<ThisType>();
  }

  AllNeumann()
  {
  }

  virtual ~AllNeumann()
  {
  }

  virtual bool dirichlet(const IntersectionType& /*intersection*/) const DS_OVERRIDE DS_FINAL
  {
    return false;
  }

  virtual bool neumann(const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    return intersection.boundary();
  }
}; // class AllNeumann


template <class IntersectionImp>
class IdBased : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef IdBased<IntersectionImp> ThisType;

public:
  using typename BaseType::IntersectionType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".idbased";
  }

  static Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Common::ConfigTree config("type", static_id());
    config["dirichlet"] = "[1 2 3]";
    config["neumann"] = "[4]";
    if (sub_name.empty())
      return config;
    else {
      Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::ConfigTree config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::ConfigTree cfg = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const auto dirichlet_ids     = cfg.get<std::vector<int>>("dirichlet", 0);
    const auto neumann_ids       = cfg.get<std::vector<int>>("neumann", 0);
    std::map<std::string, std::set<int>> id_to_type_map;
    id_to_type_map.insert(std::make_pair("dirichlet", std::set<int>(dirichlet_ids.begin(), dirichlet_ids.end())));
    id_to_type_map.insert(std::make_pair("neumann", std::set<int>(neumann_ids.begin(), neumann_ids.end())));
    return Common::make_unique<ThisType>(id_to_type_map);
  }

  IdBased(const std::map<std::string, std::set<int>> id_to_type_map)
    : id_to_type_map_(id_to_type_map)
    , hasDirichlet_(id_to_type_map_.find("dirichlet") != id_to_type_map_.end())
    , hasNeumann_(id_to_type_map_.find("neumann") != id_to_type_map_.end())
  {
  }

  virtual ~IdBased()
  {
  }

  const std::map<std::string, std::set<int>>& id_to_type_map() const
  {
    return id_to_type_map_;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    if (hasDirichlet_ && intersection.boundary()) {
#include <dune/stuff/common/disable_warnings.hh>
      const int boundaryId = intersection.boundaryId();
#include <dune/stuff/common/reenable_warnings.hh>
      // get set of dirichlet ids (has to be found, otherwise hasDirichlet_ would be false)
      const auto result = id_to_type_map_.find("dirichlet");
      assert(result != id_to_type_map_.end());
      const auto& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    if (hasNeumann_ && intersection.boundary()) {
// get boundary id
#include <dune/stuff/common/disable_warnings.hh>
      const int boundaryId = intersection.boundaryId();
#include <dune/stuff/common/reenable_warnings.hh>
      // get set of neumann ids (has to be found, otherwise hasNeumann_ would be false)
      const auto result = id_to_type_map_.find("neumann");
      assert(result != id_to_type_map_.end());
      const auto& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool neumann(const IntersectionType& intersection) const

private:
  std::map<std::string, std::set<int>> id_to_type_map_;
  bool hasDirichlet_;
  bool hasNeumann_;
}; // class IdBased


template <class IntersectionImp>
class NormalBased : public Stuff::Grid::BoundaryInfoInterface<IntersectionImp>
{
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionImp> BaseType;
  typedef NormalBased<IntersectionImp> ThisType;

public:
  using typename BaseType::IntersectionType;

  typedef typename IntersectionType::ctype DomainFieldType;
  static const unsigned int dimDomain = IntersectionType::dimension;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  static const std::string static_id()
  {
    return BaseType::static_id() + ".normalbased";
  }

  static Common::ConfigTree default_config(const std::string sub_name = "")
  {
    Common::ConfigTree config("type", static_id());
    config["default"]           = "dirichlet";
    config["compare_tolerance"] = "1e-10";
    config["neumann.0"]         = "[1.0 0.0]";
    config["dirichlet.0"] = "[0.0 1.0]";
    if (sub_name.empty())
      return config;
    else {
      Common::ConfigTree tmp;
      tmp.add(config, sub_name);
      return tmp;
    }
  } // ... default_config(...)

  static std::unique_ptr<ThisType> create(const Common::ConfigTree config = default_config(),
                                          const std::string sub_name = static_id())
  {
    const Common::ConfigTree cfg         = config.has_sub(sub_name) ? config.sub(sub_name) : config;
    const Common::ConfigTree default_cfg = default_config();
    // get default
    const std::string default_type = cfg.get("default", default_cfg.get<std::string>("default"));
    if (default_type != "dirichlet" && default_type != "neumann")
      DUNE_THROW_COLORFULLY(Exceptions::configuration_error, "Wrong default '" << default_type << "' given!");
    const bool default_to_dirichlet = default_type == "dirichlet";
    // get tolerance
    const DomainFieldType tol = cfg.get("compare_tolerance", default_cfg.get<DomainFieldType>("compare_tolerance"));
    // get dirichlet and neumann
    std::vector<DomainType> dirichlets = getVectors(cfg, "dirichlet");
    std::vector<DomainType> neumanns   = getVectors(cfg, "neumann");
    // return
    return Common::make_unique<ThisType>(default_to_dirichlet, dirichlets, neumanns, tol);
  } // ... create(...)

  NormalBased(const bool default_to_dirichlet = true,
              const std::vector<DomainType> dirichlet_normals = std::vector<DomainType>(),
              const std::vector<DomainType> neumann_normals   = std::vector<DomainType>(),
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
        DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                              "Given normals are too close for given tolerance '" << tol << "'!");
    }
  } // NormalBased(...)

  virtual ~NormalBased()
  {
  }

  virtual bool dirichlet(const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    if (intersection.boundary()) {
      const DomainType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, dirichlet_normals_))
        return true;
      else if (contains(outerNormal, neumann_normals_))
        return false;
      else
        return default_to_dirichlet_;
    }
    return false;
  } // ... dirichlet(...)

  virtual bool neumann(const IntersectionType& intersection) const DS_OVERRIDE DS_FINAL
  {
    if (intersection.boundary()) {
      const DomainType outerNormal = intersection.centerUnitOuterNormal();
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
  static std::vector<DomainType> getVectors(const Common::ConfigTree& config, const std::string key)
  {
    std::vector<DomainType> ret;
    if (config.has_sub(key)) {
      bool found           = true;
      unsigned int counter = 0;
      while (found) {
        const std::string localKey = key + "." + Dune::Stuff::Common::toString(counter);
        if (config.has_key(localKey))
          ret.push_back(config.get<DomainType>(localKey, dimDomain));
        else
          found = false;
        ++counter;
      }
    } else if (config.has_key(key))
      ret.push_back(config.get<DomainType>(key, dimDomain));
    return ret;
  } // ... getVectors(...)

  bool contains(const DomainType& normal, const std::vector<DomainType>& vectors) const
  {
    for (auto& vector : vectors)
      if (Dune::Stuff::Common::FloatCmp::eq(normal, vector, tol_))
        return true;
    return false;
  }

  const bool default_to_dirichlet_;
  std::vector<DomainType> dirichlet_normals_;
  std::vector<DomainType> neumann_normals_;
  const DomainFieldType tol_;
}; // class NormalBased


} // namespace BoundaryInfos


template <class I>
class BoundaryInfoProvider
{
public:
  typedef BoundaryInfoInterface<I> InterfaceType;

  static std::vector<std::string> available()
  {
    using namespace Stuff::Grid::BoundaryInfos;
    return {BoundaryInfos::AllDirichlet<I>::static_id(),
            BoundaryInfos::AllNeumann<I>::static_id(),
            BoundaryInfos::IdBased<I>::static_id(),
            BoundaryInfos::NormalBased<I>::static_id()};
  } // ... available(...)

  static Common::ConfigTree default_config(const std::string type, const std::string subname = "")
  {
    using namespace Stuff::Grid::BoundaryInfos;
    if (type == BoundaryInfos::AllDirichlet<I>::static_id())
      return BoundaryInfos::AllDirichlet<I>::default_config(subname);
    else if (type == BoundaryInfos::AllNeumann<I>::static_id())
      return BoundaryInfos::AllNeumann<I>::default_config(subname);
    else if (type == BoundaryInfos::IdBased<I>::static_id())
      return BoundaryInfos::IdBased<I>::default_config(subname);
    else if (type == BoundaryInfos::NormalBased<I>::static_id())
      return BoundaryInfos::NormalBased<I>::default_config(subname);
    else
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... default_config(...)

  static std::unique_ptr<InterfaceType> create(const std::string& type = available()[0],
                                               const Common::ConfigTree config = default_config(available()[0]))
  {
    using namespace Stuff::Grid::BoundaryInfos;
    if (type == BoundaryInfos::AllDirichlet<I>::static_id())
      return BoundaryInfos::AllDirichlet<I>::create(config);
    else if (type == BoundaryInfos::AllNeumann<I>::static_id())
      return BoundaryInfos::AllNeumann<I>::create(config);
    else if (type == BoundaryInfos::IdBased<I>::static_id())
      return BoundaryInfos::IdBased<I>::create(config);
    else if (type == BoundaryInfos::NormalBased<I>::static_id())
      return BoundaryInfos::NormalBased<I>::create(config);
    else
      DUNE_THROW_COLORFULLY(Exceptions::wrong_input_given,
                            "'" << type << "' is not a valid " << InterfaceType::static_id() << "!");
  } // ... create(...)
}; // class BoundaryInfoProvider


} // namespace Grid


template <class IntersectionImp>
class DUNE_DEPRECATED_MSG("Use Stuff::Grid::BoundaryInfoInterface instead!") GridboundaryInterface
#if HAVE_DUNE_PDELAB
    : public TypeTree::LeafNode // makes this usable as BoundaryCondition type directly in PDElab
#endif
{
#if HAVE_DUNE_PDELAB
  typedef PDELab::IntersectionGeometry<IntersectionImp> IntersectionGeometryType;
  typedef FieldVector<typename IntersectionGeometryType::ctype, IntersectionGeometryType::dimension - 1>
      Codim1DomainType;
#endif

public:
  typedef IntersectionImp IntersectionType;

  static const std::string id()
  {
    return "boundaryinfo";
  }

  virtual ~GridboundaryInterface()
  {
  }

  virtual bool dirichlet(const IntersectionType&) const = 0;
  virtual bool neumann(const IntersectionType&) const = 0;

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
}; // class GridboundaryInterface


template <class IntersectionImp>
class DUNE_DEPRECATED_MSG("Use Stuff::Grid::BoundaryInfos::AllDirichlet instead!") GridboundaryAllDirichlet
    : public GridboundaryInterface<IntersectionImp>
{
  typedef GridboundaryInterface<IntersectionImp> BaseType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".alldirichlet";
  }

  static Dune::ParameterTree defaultSettings(const std::string /*subName*/ = "")
  {
    return Dune::ParameterTree();
  }

  GridboundaryAllDirichlet()
  {
  }

  virtual bool dirichlet(const IntersectionType& intersection) const DS_OVERRIDE
  {
    return intersection.boundary();
  } // virtual bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType&) const DS_OVERRIDE
  {
    return false;
  }
}; // class GridboundaryAllDirichlet


template <class IntersectionImp>
class DUNE_DEPRECATED_MSG("Use Stuff::Grid::BoundaryInfos::AllNeumann instead!") GridboundaryAllNeumann
    : public GridboundaryInterface<IntersectionImp>
{
  typedef GridboundaryInterface<IntersectionImp> BaseType;

public:
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".allneumann";
  }

  static Dune::ParameterTree defaultSettings(const std::string /*subName*/ = "")
  {
    return Dune::ParameterTree();
  }

  GridboundaryAllNeumann()
  {
  }

  virtual bool dirichlet(const IntersectionType& /*intersection*/) const DS_OVERRIDE
  {
    return false;
  }

  virtual bool neumann(const IntersectionType& intersection) const DS_OVERRIDE
  {
    return intersection.boundary();
  } // virtual bool neumann(const IntersectionType& intersection) const
}; // class GridboundaryAllNeumann


template <class IntersectionImp>
class DUNE_DEPRECATED_MSG("Use Stuff::Grid::BoundaryInfos::IdBased instead!") GridboundaryIdBased
    : public GridboundaryInterface<IntersectionImp>
{
  typedef GridboundaryInterface<IntersectionImp> BaseType;

public:
  typedef GridboundaryIdBased<IntersectionImp> ThisType;
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".idbased";
  }

  typedef int IdType;
  typedef std::set<IdType> IdSetType;
  typedef std::map<std::string, IdSetType> IdSetMapType;

  GridboundaryIdBased(const IdSetMapType boundaryInfoMap)
    : boundaryInfoMap_(boundaryInfoMap)
    , hasDirichlet_(boundaryInfoMap_.find("dirichlet") != boundaryInfoMap_.end())
    , hasNeumann_(boundaryInfoMap_.find("neumann") != boundaryInfoMap_.end())
  {
  }

  GridboundaryIdBased(const Dune::ParameterTree& settings, const std::string subName = id())
  {
    const Common::ExtendedParameterTree ext_settings = settings.hasSub(subName) ? settings.sub(subName) : settings;
    const std::vector<int> dirichletIds              = ext_settings.getVector<int>("dirichlet", 0, 0);
    const std::vector<int> neumannIds = ext_settings.getVector<int>("neumann", 0, 0);
    boundaryInfoMap_.insert(
        std::pair<std::string, IdSetType>("dirichlet", IdSetType(dirichletIds.begin(), dirichletIds.end())));
    boundaryInfoMap_.insert(
        std::pair<std::string, IdSetType>("neumann", IdSetType(neumannIds.begin(), neumannIds.end())));
    hasDirichlet_ = boundaryInfoMap_.find("dirichlet") != boundaryInfoMap_.end();
    hasNeumann_   = boundaryInfoMap_.find("neumann") != boundaryInfoMap_.end();
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["dirichlet"] = "[1 2 3]";
    description["neumann"] = "[4]";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  }

  static ThisType* create(const Dune::ParameterTree& settings, const std::string subName = id())
  {
    return new ThisType(settings, subName);
  }

  const IdSetMapType& boundaryInfoMap() const
  {
    return boundaryInfoMap_;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const DS_OVERRIDE
  {
    if (hasDirichlet_ && intersection.boundary()) {
      const IdType boundaryId = intersection.boundaryId();
      // get set of dirichlet ids (has to be found, otherwise hasDirichlet_ would be false)
      const auto result = boundaryInfoMap_.find("dirichlet");
      assert(result != boundaryInfoMap_.end());
      const IdSetType& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType& intersection) const DS_OVERRIDE
  {
    if (hasNeumann_ && intersection.boundary()) {
      // get boundary id
      const IdType boundaryId = intersection.boundaryId();
      // get set of neumann ids (has to be found, otherwise hasNeumann_ would be false)
      const typename IdSetMapType::const_iterator result = boundaryInfoMap_.find("neumann");
      assert(result != boundaryInfoMap_.end());
      const IdSetType& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool neumann(const IntersectionType& intersection) const

private:
  IdSetMapType boundaryInfoMap_;
  bool hasDirichlet_;
  bool hasNeumann_;
}; // class GridboundaryIdBased


template <class IntersectionImp>
class DUNE_DEPRECATED_MSG("Use Stuff::Grid::BoundaryInfos::NormalBased instead!") GridboundaryNormalBased
    : public GridboundaryInterface<IntersectionImp>
{
  typedef GridboundaryInterface<IntersectionImp> BaseType;

public:
  typedef GridboundaryNormalBased<IntersectionImp> ThisType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef typename IntersectionType::ctype DomainFieldType;
  static const unsigned int dimDomain = IntersectionType::dimension;
  typedef Dune::FieldVector<DomainFieldType, dimDomain> DomainType;

  static const std::string id()
  {
    return BaseType::id() + ".normalbased";
  }

  GridboundaryNormalBased(const bool defaultIsDirichlet = true,
                          const std::vector<DomainType> dirichletNormals = std::vector<DomainType>(),
                          const std::vector<DomainType> neumannNormals   = std::vector<DomainType>(),
                          const DomainFieldType tol = 1e-10)
    : defaultIsDirichlet_(defaultIsDirichlet)
    , dirichletNormals_(dirichletNormals)
    , neumannNormals_(neumannNormals)
    , tol_(tol)
  {
    // normalize
    for (auto& normal : dirichletNormals_)
      normal /= normal.two_norm();
    for (auto& normal : neumannNormals_)
      normal /= normal.two_norm();
    // sanity check
    for (auto& dirichletNormal : dirichletNormals_) {
      if (contains(dirichletNormal, neumannNormals_))
        DUNE_THROW(Dune::InvalidStateException,
                   "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " normals are too close!");
    }
  }

  GridboundaryNormalBased(const ThisType& other)
    : defaultIsDirichlet_(other.defaultIsDirichlet_)
    , dirichletNormals_(other.dirichletNormals_)
    , neumannNormals_(other.neumannNormals_)
    , tol_(other.tol_)
  {
  }

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      defaultIsDirichlet_ = other.defaultIsDirichlet_;
      dirichletNormals_   = other.dirichletNormals_;
      neumannNormals_     = other.neumannNormals_;
      tol_                = other.tol_;
    }
    return this;
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["default"]           = "dirichlet";
    description["compare_tolerance"] = "1e-10";
    description["neumann.0"]         = "[1.0 0.0]";
    description["dirichlet.0"] = "[0.0 1.0]";
    if (subName.empty())
      return description;
    else {
      Dune::Stuff::Common::ExtendedParameterTree extendedDescription;
      extendedDescription.add(description, subName);
      return extendedDescription;
    }
  }

  static ThisType* create(const Dune::ParameterTree& _settings, const std::string subName = id())
  {
    // get correct _settings
    Common::ExtendedParameterTree settings;
    if (_settings.hasSub(subName))
      settings = _settings.sub(subName);
    else
      settings = _settings;
    // get default
    bool dirichletDef     = false;
    const std::string def = settings.get("default", "dirichlet");
    if (def == "dirichlet")
      dirichletDef = true;
    else if (def == "neumann")
      dirichletDef = false;
    else
      DUNE_THROW(Dune::IOError, "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " wrong 'default' given!");
    // get tolerance
    const DomainFieldType tol = settings.get("compare_tolerance", 1e-10);
    // get dirichlet and neumann
    std::vector<DomainType> dirichlets = getVectors(settings, "dirichlet");
    std::vector<DomainType> neumanns   = getVectors(settings, "neumann");
    // return
    return new ThisType(dirichletDef, dirichlets, neumanns, tol);
  }

  virtual bool dirichlet(const IntersectionType& intersection) const DS_OVERRIDE
  {
    if (intersection.boundary()) {
      const DomainType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, dirichletNormals_))
        return true;
      else if (contains(outerNormal, neumannNormals_))
        return false;
      else
        return defaultIsDirichlet_;
    }
    return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType& intersection) const DS_OVERRIDE
  {
    if (intersection.boundary()) {
      const DomainType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, neumannNormals_))
        return true;
      else if (contains(outerNormal, dirichletNormals_))
        return false;
      else
        return !defaultIsDirichlet_;
    }
    return false;
  } // bool neumann(const IntersectionType& intersection) const

private:
  static std::vector<DomainType> getVectors(const Common::ExtendedParameterTree& paramTree, const std::string key)
  {
    std::vector<DomainType> ret;
    if (paramTree.hasSub(key)) {
      bool found           = true;
      unsigned int counter = 0;
      while (found) {
        const std::string localKey = key + "." + Dune::Stuff::Common::toString(counter);
        if (paramTree.hasKey(localKey)) {
          std::vector<DomainFieldType> vec = paramTree.getVector<DomainFieldType>(localKey, dimDomain);
          DomainType fvec(0);
          for (size_t dd = 0; dd < dimDomain; ++dd)
            fvec[dd] = vec[dd];
          ret.push_back(fvec);
        } else
          found = false;
        ++counter;
      }
    } else if (paramTree.hasKey(key)) {
      std::vector<DomainFieldType> vec = paramTree.getVector<DomainFieldType>(key, dimDomain);
      DomainType fvec(0);
      for (size_t dd = 0; dd < dimDomain; ++dd)
        fvec[dd] = vec[dd];
      ret.push_back(fvec);
    }
    return ret;
  }

  bool contains(const DomainType& normal, const std::vector<DomainType>& vectors) const
  {
    for (auto& vector : vectors)
      if (Dune::Stuff::Common::FloatCmp::eq(normal, vector, tol_))
        return true;
    return false;
  }

  const bool defaultIsDirichlet_;
  std::vector<DomainType> dirichletNormals_;
  std::vector<DomainType> neumannNormals_;
  const DomainFieldType tol_;
}; // class GridboundaryNormalBased


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
