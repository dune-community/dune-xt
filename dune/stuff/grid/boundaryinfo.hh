// This file is part of the dune-stuff project:
//   http://users.dune-project.org/projects/dune-stuff/
// Copyright Holders: Felix Albrecht, Rene Milk
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_GRID_BOUNDARYINFO_HH
#define DUNE_STUFF_GRID_BOUNDARYINFO_HH

#include <map>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/common/float_cmp.hh>

namespace Dune {
namespace Stuff {


template <class IntersectionImp>
class GridboundaryInterface
{
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
}; // class GridboundaryInterface


template <class IntersectionImp>
class GridboundaryAllDirichlet : public GridboundaryInterface<IntersectionImp>
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


template <class IntersectionImp>
class GridboundaryAllNeumann : public GridboundaryInterface<IntersectionImp>
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

  virtual bool dirichlet(const IntersectionType& /*intersection*/) const
  {
    return false;
  }

  virtual bool neumann(const IntersectionType& intersection) const
  {
    return intersection.boundary();
  } // virtual bool neumann(const IntersectionType& intersection) const
}; // class GridboundaryAllNeumann


template <class IntersectionImp>
class GridboundaryIdBased : public GridboundaryInterface<IntersectionImp>
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

  static ThisType* create(const Dune::ParameterTree& settings, const std::string subName = id())
  {
    return new ThisType(settings, subName);
  }

  const IdSetMapType& boundaryInfoMap() const
  {
    return boundaryInfoMap_;
  }

  virtual bool dirichlet(const IntersectionType& intersection) const
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

  virtual bool neumann(const IntersectionType& intersection) const
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
class GridboundaryNormalBased : public GridboundaryInterface<IntersectionImp>
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
    normalize();
  }

  GridboundaryNormalBased(const Dune::ParameterTree& _settings, const std::string subName = id())
  {
    const Common::ExtendedParameterTree settings = _settings.hasSub(subName) ? _settings.sub(subName) : _settings;
    const std::string def = settings.get("default", "dirichlet");
    if (def == "dirichlet")
      defaultIsDirichlet_ = true;
    else if (def == "neumann")
      defaultIsDirichlet_ = false;
    else
      DUNE_THROW(Dune::IOError, "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " wrong 'default' given!");
    tol_              = settings.get("compare_tolerance", 1e-10);
    dirichletNormals_ = getVectors(settings, "dirichlet");
    neumannNormals_ = getVectors(settings, "neumann");
    normalize();
  }

  static Dune::ParameterTree defaultSettings(const std::string subName = "")
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

  static ThisType* create(const Dune::ParameterTree& settings, const std::string subName = id())
  {
    return new ThisType(settings, subName);
  }

  virtual bool dirichlet(const IntersectionType& intersection) const
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

  virtual bool neumann(const IntersectionType& intersection) const
  {
    return !dirichlet(intersection);
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

  void normalize()
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

  bool defaultIsDirichlet_;
  std::vector<DomainType> dirichletNormals_;
  std::vector<DomainType> neumannNormals_;
  DomainFieldType tol_;
}; // class GridboundaryNormalBased


template <class IntersectionType>
class Gridboundaries
{
public:
  static std::vector<std::string> available()
  {
    return {"boundaryinfo.alldirichlet", "boundaryinfo.allneumann", "boundaryinfo.idbased", "boundaryinfo.normalbased"};
  } // ... available(...)

  static Dune::ParameterTree defaultSettings(const std::string type, const std::string subname = "")
  {
    if (type == "boundaryinfo.alldirichlet")
      return GridboundaryAllDirichlet<IntersectionType>::defaultSettings(subname);
    else if (type == "boundaryinfo.allneumann")
      return GridboundaryAllNeumann<IntersectionType>::defaultSettings(subname);
    else if (type == "boundaryinfo.idbased")
      return GridboundaryIdBased<IntersectionType>::defaultSettings(subname);
    else if (type == "boundaryinfo.normalbased")
      return GridboundaryNormalBased<IntersectionType>::defaultSettings(subname);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown boundaryinfo '" << type
                      << "' requested!");
  } // ... createDefaultSettings(...)

  static GridboundaryInterface<IntersectionType>* create(const std::string& type = available()[0],
                                                         const Dune::ParameterTree settings = Dune::ParameterTree())
  {
    if (type == "boundaryinfo.alldirichlet") {
      return new GridboundaryAllDirichlet<IntersectionType>();
    } else if (type == "boundaryinfo.allneumann") {
      return new GridboundaryAllNeumann<IntersectionType>();
    } else if (type == "boundaryinfo.idbased") {
      return GridboundaryIdBased<IntersectionType>::create(settings);
    } else if (type == "boundaryinfo.normalbased") {
      return GridboundaryNormalBased<IntersectionType>::create(settings);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown boundaryinfo '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Gridboundaries


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
