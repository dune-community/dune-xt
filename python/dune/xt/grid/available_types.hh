#ifndef PYTHON_DUNE_XT_GRID_TYPES_HH
#define PYTHON_DUNE_XT_GRID_TYPES_HH

#include <dune/xt/grid/grids.hh>

#include <boost/tuple/tuple.hpp>

namespace Dune {
namespace XT {
namespace Grid {
namespace bindings {

using AvailableTypes = boost::tuple<Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>,
                                    Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>
#if HAVE_DUNE_ALUGRID
                                    ,
                                    Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                    ,
                                    Dune::UGGrid<2>
#endif
#if HAVE_ALBERTA
                                    ,
                                    Dune::AlbertaGrid<2, 2>
#endif
                                    >;
} // namespace bindings
} // namespace Grid
} // namespace XT
} // namespace Dune

#endif // PYTHON_DUNE_XT_GRID_TYPES_HH