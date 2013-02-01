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
#include <dune/common/exceptions.hh>

#include <dune/stuff/common/parameter/tree.hh>

namespace Dune {
namespace Stuff {
namespace Grid {
namespace BoundaryInfo {

template <class GridViewImp>
class Interface
{
public:
  typedef GridViewImp GridViewType;

  typedef typename GridViewType::Intersection IntersectionType;

  static const std::string id()
  {
    return "stuff.grid.boundaryinfo";
  }

  virtual bool dirichlet(const IntersectionType&) const = 0;

  virtual bool neumann(const IntersectionType&) const = 0;
}; // class Interface


template <class GridViewImp>
class AllDirichlet : public Interface<GridViewImp>
{
public:
  typedef Interface<GridViewImp> BaseType;

  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".alldirichlet";
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
}; // class AllDirichlet

template <class GridViewImp>
class AllNeumann : public Interface<GridViewImp>
{
public:
  typedef Interface<GridViewImp> BaseType;

  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".allneumann";
  }

  virtual bool dirichlet(const IntersectionType& intersection) const
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
}; // class AllNeumann

template <class GridViewImp>
class IdBased : public Interface<GridViewImp>
{
public:
  typedef Interface<GridViewImp> BaseType;

  typedef typename BaseType::IntersectionType IntersectionType;

  static const std::string id()
  {
    return BaseType::id() + ".idbased";
  }

  typedef IdBased ThisType;

  typedef int IdType;

  typedef std::set<IdType> IdSetType;

  typedef std::map<std::string, IdSetType> IdSetMapType;

  IdBased(const Dune::shared_ptr<const IdSetMapType> _boundaryInfoMap)
    : boundaryInfoMap_(_boundaryInfoMap)
  {
    setup();
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
}; // class IdBased

template <class GridViewType>
Dune::shared_ptr<Interface<GridViewType>> create(const std::string& type = "stuff.grid.boundaryinfo.alldirichlet",
                                                 const Dune::ParameterTree paramTree = Dune::ParameterTree())
{
  if (type == "stuff.grid.boundaryinfo.alldirichlet") {
    typedef AllDirichlet<GridViewType> AllDirichletType;
    return Dune::make_shared<AllDirichletType>();
  } else if (type == "stuff.grid.boundaryinfo.allneumann") {
    typedef AllNeumann<GridViewType> AllNeumannType;
    return Dune::make_shared<AllNeumannType>();
  } else if (type == "stuff.grid.boundaryinfo.idbased") {
    typedef IdBased<GridViewType> IdBasedType;
    return Dune::make_shared<IdBasedType>(IdBasedType::createFromParamTree(paramTree));
  } else
    DUNE_THROW(Dune::RangeError, "\nERROR: unknown boundaryinfo '" << type << "' requested!");
} // Interface< GridViewType >* create(...)

} // namespace BoundaryInfo
} // namespace Grid
} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
