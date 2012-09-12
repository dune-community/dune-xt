
#ifndef DUNE_STUFF_GRID_BOUNDARYINFO_HH
#define DUNE_STUFF_GRID_BOUNDARYINFO_HH

// system
#include <map>
#include <set>

// dune-common
#include <dune/common/shared_ptr.hh>

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
  bool neumann(const IntersectionType& intersection) const
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
  typedef int IdType;

  typedef std::set<IdType> IdSetType;

  typedef std::map<std::string, IdSetType> IdSetMapType;

  IdBased(const Dune::shared_ptr<const IdSetMapType> boundaryInfoMap)
    : boundaryInfoMap_(boundaryInfoMap)
    , hasDirichlet_(boundaryInfoMap_->find("dirichlet") != boundaryInfoMap_->end())
    , hasNeumann_(boundaryInfoMap_->find("neumann") != boundaryInfoMap_->end())
  {
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
  const Dune::shared_ptr<const IdSetMapType> boundaryInfoMap_;
  const bool hasDirichlet_;
  const bool hasNeumann_;
}; // class IdBased

} // namespace BoundaryInfo

} // namespace Grid

} // namespace Stuff

} // namespace Dune

#endif // DUNE_STUFF_GRID_BOUNDARYINFO_HH
