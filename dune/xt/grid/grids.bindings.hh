#ifndef DUNE_XT_GRID_GRIDS_BINDINGS_HH
#define DUNE_XT_GRID_GRIDS_BINDINGS_HH

#include "grids.hh"

namespace Dune {


// this is used by other headers
typedef YaspGrid<2, EquidistantOffsetCoordinates<double, 2>> YASP_2D_EQUIDISTANT_OFFSET;
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
typedef ALUGrid<2, 2, simplex, conforming> ALU_2D_SIMPLEX_CONFORMING;
#endif


} // namespace Dune

#endif // DUNE_XT_GRID_GRIDS_BINDINGS_HH
