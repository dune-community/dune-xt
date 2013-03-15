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

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/color.hh>

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

  static Dune::ParameterTree createSampleDescription(const std::string /*subName*/ = "")
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

  static Dune::ParameterTree createSampleDescription(const std::string /*subName*/ = "")
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

  static Dune::ParameterTree createSampleDescription(const std::string subName = "")
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

  static ThisType* create(const Dune::ParameterTree& paramTree, const std::string subName = id())
  {
    // get correct paramTree
    Common::ExtendedParameterTree paramTreeX;
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


template <class GridViewType>
class Gridboundaries
{
public:
  static std::vector<std::string> available()
  {
    return {"boundaryinfo.alldirichlet", "boundaryinfo.allneumann", "boundaryinfo.idbased"};
  } // ... available(...)

  static Dune::ParameterTree createSampleDescription(const std::string type, const std::string subname = "")
  {
    if (type == "boundaryinfo.alldirichlet")
      return GridboundaryAllDirichlet<GridViewType>::createSampleDescription(subname);
    else if (type == "boundaryinfo.allneumann")
      return GridboundaryAllNeumann<GridViewType>::createSampleDescription(subname);
    else if (type == "boundaryinfo.idbased")
      return GridboundaryIdBased<GridViewType>::createSampleDescription(subname);
    else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown boundaryinfo '" << type
                      << "' requested!");
  } // ... createSampleDescription(...)

  static GridboundaryInterface<GridViewType>* create(const std::string& type = available()[0],
                                                     const Dune::ParameterTree description = Dune::ParameterTree())
  {
    if (type == "boundaryinfo.alldirichlet") {
      return new GridboundaryAllDirichlet<GridViewType>();
    } else if (type == "boundaryinfo.allneumann") {
      return new GridboundaryAllNeumann<GridViewType>();
    } else if (type == "boundaryinfo.idbased") {
      return GridboundaryIdBased<GridViewType>::create(description);
    } else
      DUNE_THROW(Dune::RangeError,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:") << " unknown boundaryinfo '" << type
                      << "' requested!");
  } // ... create(...)
}; // class Gridboundaries


} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
