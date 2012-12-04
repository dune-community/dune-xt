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

#include <dune/common/shared_ptr.hh>

#include <dune/stuff/common/parameter/tree.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace BoundaryInfo {

class AllDirichlet
{
public:
  template <class IntersectionType>
  bool dirichlet(const IntersectionType& intersection) const
  {
    if (intersection.boundary())
      return true;
    else
      return false;
  }

  template <class IntersectionType>
  bool neumann(const IntersectionType&) const
  {
    return false;
  }
}; // class AllDirichlet

class AllNeumann
{
public:
  template <class IntersectionType>
  bool dirichlet(const IntersectionType& intersection) const
  {
    return false;
  }

  template <class IntersectionType>
  bool neumann(const IntersectionType& intersection) const
  {
    if (intersection.boundary())
      return true;
    else
      return false;
  }
}; // class AllNeumann

class IdBased
{
public:
  typedef IdBased ThisType;

  typedef int IdType;

  typedef std::set<IdType> IdSetType;

  typedef std::map<std::string, IdSetType> IdSetMapType;

  static const std::string id()
  {
    return "stuff.boundaryinfo.idbased";
  }

  IdBased(const Dune::shared_ptr<const IdSetMapType> _boundaryInfoMap)
    : boundaryInfoMap_(_boundaryInfoMap)
  {
  }

  IdBased(const ThisType& other)
    : boundaryInfoMap_(other.boundaryInfoMap_)
    , hasDirichlet_(other.hasDirichlet_)
    , hasNeumann_(other.hasNeumann_)
  {
  }

  static ThisType createFromParamTree(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Common::ParameterTreeX paramTreeX;
    if (paramTree.hasSub(subName))
      paramTreeX = paramTree.sub(subName);
    else
      paramTreeX = paramTree;
    // get dirichlet
    const std::vector<int> dirichletIds = paramTreeX.getVector<int>("dirichlet", 0, 0);
    IdSetType dirichletSet;
    for (unsigned int i = 0; i < dirichletIds.size(); ++i)
      dirichletSet.insert(dirichletIds[i]);
    // get neumann
    const std::vector<int> neumannIds = paramTreeX.getVector<int>("neumann", 0, 0);
    IdSetType neumannSet;
    for (unsigned int i = 0; i < neumannIds.size(); ++i)
      neumannSet.insert(neumannIds[i]);
    // create map and return
    shared_ptr<IdSetMapType> idSetMap(new IdSetMapType());
    idSetMap->insert(std::pair<std::string, IdSetType>("dirichlet", dirichletSet));
    idSetMap->insert(std::pair<std::string, IdSetType>("neumann", neumannSet));
    return ThisType(idSetMap);
  }

  ThisType& operator=(ThisType& other)
  {
    if (this != &other) {
      boundaryInfoMap_ = other.boundaryInfoMap();
      setup();
    }
    return *this;
  }

  const Dune::shared_ptr<const IdSetMapType> boundaryInfoMap()
  {
    return boundaryInfoMap_;
  }

  template <class IntersectionType>
  bool dirichlet(const IntersectionType& intersection) const
  {
    if (hasDirichlet_ && intersection.boundary()) {
      // get boundary id
      const IdType boundaryId = intersection.boundaryId();
      // get set of dirichlet ids (has to be found, otherwise hasDirichlet_ would be false)
      const typename IdSetMapType::const_iterator result = boundaryInfoMap_->find("dirichlet");
      assert(result != boundaryInfoMap_->end());
      const IdSetType& idSet = result->second;
      return (idSet.count(boundaryId) > 0);
    } else
      return false;
  } // bool dirichlet(const IntersectionType& intersection) const

  template <class IntersectionType>
  bool neumann(const IntersectionType& intersection) const
  {
    if (hasNeumann_ && intersection.boundary()) {
      // get boundary id
      const IdType boundaryId = intersection.boundaryId();
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
}; // class IdBased

} // namespace BoundaryInfo
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
