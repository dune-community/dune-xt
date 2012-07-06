#ifndef DUNE_STUFF_GRID_INFORMATION_HH
#define DUNE_STUFF_GRID_INFORMATION_HH

namespace Dune {

namespace Stuff {

namespace Grid {

namespace Information {

/**
  * \attention Not optimal, does a whole grid walk!
  **/
template <class GridPartType>
unsigned int maxNumberOfNeighbors(const GridPartType& gridPart)
{
  // some preparations
  unsigned int maxNeighbours = 0;
  unsigned int neighbours    = 0;
  typedef typename GridPartType::template Codim<0>::IteratorType EntityIteratorType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
  // walk over all entities
  const EntityIteratorType entityIteratorEnd = gridPart.template end<0>();
  for (EntityIteratorType entityIterator = gridPart.template begin<0>(); entityIterator != entityIteratorEnd;
       ++entityIterator) {
    const EntityType& entity = *entityIterator;
    neighbours               = 0;
    // walk over all neighbors
    const IntersectionIteratorType intersectionIteratorEnd = gridPart.iend(entity);
    for (IntersectionIteratorType intersectionIterator = gridPart.ibegin(entity);
         intersectionIterator != intersectionIteratorEnd;
         ++intersectionIterator) {
      ++neighbours;
    } // walk over all neighbors
    maxNeighbours = std::max(maxNeighbours, neighbours);
  } // walk over all entities
  return maxNeighbours;
} // unsigned int maxNumberOfNeighbors(const GridPartType& gridPart)

} // namespace Information

} // namespace Grid

} // namespace Stuff

} // namespace Dune

#endif // DUNE_STUFF_GRID_INFORMATION_HH
