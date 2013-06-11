#ifndef DUNE_STUFF_GRID_BOUNDARYINFO_HH
#define DUNE_STUFF_GRID_BOUNDARYINFO_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#include <map>
#include <set>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/vector.hh>

namespace Dune {
namespace Stuff {


template <class GridViewImp>
class GridboundaryInterface
{
public:
  typedef GridViewImp GridViewType;
  typedef typename GridViewType::Intersection IntersectionType;

  static const std::string id()
  {
    return "boundaryinfo";
  }

  virtual ~GridboundaryInterface()
  {
  }

  virtual bool dirichlet(const IntersectionType&) const = 0;

  virtual bool neumann(const IntersectionType&) const = 0;
}; // class GridboundaryInterface


template <class GridViewImp>
class GridboundaryAllDirichlet : public GridboundaryInterface<GridViewImp>
{
public:
  typedef GridboundaryInterface<GridViewImp> BaseType;
  typedef GridboundaryAllDirichlet<GridViewImp> ThisType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".alldirichlet";
  }

  static Dune::ParameterTree createDefaultSettings(const std::string /*subName*/ = "")
  {
    return Dune::ParameterTree();
  }

  virtual bool dirichlet(const IntersectionType& intersection) const
  {
    if (intersection.boundary())
      return true;
    else
      return false;
  } // virtual bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType&) const
  {
    return false;
  }
}; // class GridboundaryAllDirichlet


template <class GridViewImp>
class GridboundaryAllNeumann : public GridboundaryInterface<GridViewImp>
{
public:
  typedef GridboundaryInterface<GridViewImp> BaseType;
  typedef GridboundaryAllNeumann<GridViewImp> ThisType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".allneumann";
  }

  static Dune::ParameterTree createDefaultSettings(const std::string /*subName*/ = "")
  {
    return Dune::ParameterTree();
  }

  virtual bool dirichlet(const IntersectionType& /*intersection*/) const
  {
    return false;
  }

  virtual bool neumann(const IntersectionType& intersection) const
  {
    if (intersection.boundary())
      return true;
    else
      return false;
  } // virtual bool neumann(const IntersectionType& intersection) const
}; // class GridboundaryAllNeumann


template <class GridViewImp>
class GridboundaryIdBased : public GridboundaryInterface<GridViewImp>
{
public:
  typedef GridboundaryInterface<GridViewImp> BaseType;
  typedef GridboundaryIdBased<GridViewImp> ThisType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".idbased";
  }

  typedef int IdType;
  typedef std::set<IdType> IdSetType;
  typedef std::map<std::string, IdSetType> IdSetMapType;

  GridboundaryIdBased(const Dune::shared_ptr<const IdSetMapType> _boundaryInfoMap)
    : boundaryInfoMap_(_boundaryInfoMap)
  {
    setup();
  }

  GridboundaryIdBased(const ThisType& other)
    : boundaryInfoMap_(other.boundaryInfoMap_)
    , hasDirichlet_(other.hasDirichlet_)
    , hasNeumann_(other.hasNeumann_)
  {
  }

  ThisType& operator=(ThisType& other)
  {
    if (this != &other) {
      boundaryInfoMap_ = other.boundaryInfoMap();
      setup();
    }
    return *this;
  }

  static Dune::ParameterTree createDefaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["dirichlet"] = "[1; 2; 3]";
    description["neumann"] = "[4]";
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
    // get correct settings
    Common::ExtendedParameterTree settings;
    if (_settings.hasSub(subName))
      settings = _settings.sub(subName);
    else
      settings = _settings;
    // get dirichlet
    const std::vector<int> dirichletIds = settings.getVector<int>("dirichlet", 0, 0);
    IdSetType dirichletSet;
    for (unsigned int i = 0; i < dirichletIds.size(); ++i)
      dirichletSet.insert(dirichletIds[i]);
    // get neumann
    const std::vector<int> neumannIds = settings.getVector<int>("neumann", 0, 0);
    IdSetType neumannSet;
    for (unsigned int i = 0; i < neumannIds.size(); ++i)
      neumannSet.insert(neumannIds[i]);
    // create map and return
    shared_ptr<IdSetMapType> idSetMap(new IdSetMapType());
    idSetMap->insert(std::pair<std::string, IdSetType>("dirichlet", dirichletSet));
    idSetMap->insert(std::pair<std::string, IdSetType>("neumann", neumannSet));
    return new ThisType(idSetMap);
  }

  const Dune::shared_ptr<const IdSetMapType> boundaryInfoMap()
  {
    return boundaryInfoMap_;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const
  {
    if (hasDirichlet_ && intersection.boundary()) {
// get boundary id
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      const IdType boundaryId = intersection.boundaryId();
#pragma GCC diagnostic pop
      // get set of dirichlet ids (has to be found, otherwise hasDirichlet_ would be false)
      const typename IdSetMapType::const_iterator result = boundaryInfoMap_->find("dirichlet");
      assert(result != boundaryInfoMap_->end());
      const IdSetType& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType& intersection) const
  {
    if (hasNeumann_ && intersection.boundary()) {
// get boundary id
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
      const IdType boundaryId = intersection.boundaryId();
#pragma GCC diagnostic pop
      // get set of neumann ids (has to be found, otherwise hasNeumann_ would be false)
      const typename IdSetMapType::const_iterator result = boundaryInfoMap_->find("neumann");
      assert(result != boundaryInfoMap_->end());
      const IdSetType& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool neumann(const IntersectionType& intersection) const

private:
  void setup()
  {
    hasDirichlet_ = boundaryInfoMap_->find("dirichlet") != boundaryInfoMap_->end();
    hasNeumann_   = boundaryInfoMap_->find("neumann") != boundaryInfoMap_->end();
  }

  Dune::shared_ptr<const IdSetMapType> boundaryInfoMap_;
  bool hasDirichlet_;
  bool hasNeumann_;
}; // class GridboundaryIdBased


template <class GridViewImp>
class GridboundaryNormalBased : public GridboundaryInterface<GridViewImp>
{
public:
  typedef GridboundaryInterface<GridViewImp> BaseType;
  typedef GridboundaryNormalBased<GridViewImp> ThisType;

  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::IntersectionType IntersectionType;

  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
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

  static Dune::ParameterTree createDefaultSettings(const std::string subName = "")
  {
    Dune::ParameterTree description;
    description["default"]           = "dirichlet";
    description["compare_tolerance"] = "1e-10";
    description["neumann.0"]         = "[1.0; 0.0]";
    description["dirichlet.0"] = "[0.0; 1.0]";
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
      DUNE_THROW(Dune::IOError, "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " wrong 'default'' given!");
    // get tolerance
    const DomainFieldType tol = settings.get("compare_tolerance", 1e-10);
    // get dirichlet and neumann
    std::vector<DomainType> dirichlets = getVectors(settings, "dirichlet");
    std::vector<DomainType> neumanns   = getVectors(settings, "neumann");
    // return
    return new ThisType(dirichletDef, dirichlets, neumanns, tol);
  }

  virtual bool dirichlet(const IntersectionType& intersection) const
  {
    if (intersection.boundary()) {
      const DomainType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, dirichletNormals_))
        return true;
      else if (contains(outerNormal, neumannNormals_))
        return false;
      else if (defaultIsDirichlet_)
        return true;
    }
    return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  virtual bool neumann(const IntersectionType& intersection) const
  {
    if (intersection.boundary()) {
      const DomainType outerNormal = intersection.centerUnitOuterNormal();
      if (contains(outerNormal, neumannNormals_))
        return true;
      else if (contains(outerNormal, dirichletNormals_))
        return false;
      else if (!defaultIsDirichlet_)
        return true;
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
          std::vector<DomainFieldType> vec = paramTree.getVector<DomainFieldType>(key, dimDomain);
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
      if (Dune::Stuff::Common::float_cmp(normal, vector, tol_))
        return true;
    return false;
  }

  const bool defaultIsDirichlet_;
  std::vector<DomainType> dirichletNormals_;
  std::vector<DomainType> neumannNormals_;
  const DomainFieldType tol_;
}; // class GridboundaryNormalBased


template <class GridViewType>
class Gridboundaries
{
public:
  static std::vector<std::string> available()
  {
    return {"boundaryinfo.alldirichlet", "boundaryinfo.allneumann", "boundaryinfo.idbased", "boundaryinfo.normalbased"};
  } // ... available(...)

  static Dune::ParameterTree createDefaultSettings(const std::string type, const std::string subname = "")
  {
    if (type == "boundaryinfo.alldirichlet")
      return GridboundaryAllDirichlet<GridViewType>::createDefaultSettings(subname);
    else if (type == "boundaryinfo.allneumann")
      return GridboundaryAllNeumann<GridViewType>::createDefaultSettings(subname);
    else if (type == "boundaryinfo.idbased")
      return GridboundaryIdBased<GridViewType>::createDefaultSettings(subname);
    else if (type == "boundaryinfo.normalbased")
      return GridboundaryNormalBased<GridViewType>::createDefaultSettings(subname);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown boundaryinfo '" << type
                      << "' requested!");
  } // ... createDefaultSettings(...)

  static GridboundaryInterface<GridViewType>* create(const std::string& type = available()[0],
                                                     const Dune::ParameterTree settings = Dune::ParameterTree())
  {
    if (type == "boundaryinfo.alldirichlet") {
      return new GridboundaryAllDirichlet<GridViewType>();
    } else if (type == "boundaryinfo.allneumann") {
      return new GridboundaryAllNeumann<GridViewType>();
    } else if (type == "boundaryinfo.idbased") {
      return GridboundaryIdBased<GridViewType>::create(settings);
    } else if (type == "boundaryinfo.normalbased") {
      return GridboundaryNormalBased<GridViewType>::create(settings);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown boundaryinfo '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Gridboundaries


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
