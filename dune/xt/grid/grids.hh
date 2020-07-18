// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2009-2020 dune-xt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2019)
//   René Fritze     (2017 - 2019)
//   Tim Keil        (2018)
//   Tobias Leibner  (2016, 2020)

#ifndef DUNE_XT_GRID_GRIDS_HH
#define DUNE_XT_GRID_GRIDS_HH

#include <boost/tuple/tuple.hpp>

#if HAVE_ALBERTA
#  include <dune/xt/common/disable_warnings.hh>
#  include <dune/grid/albertagrid.hh>
#  include <dune/xt/common/reenable_warnings.hh>
#endif

#if HAVE_DUNE_ALUGRID
#  include <dune/alugrid/grid.hh>
#endif

#if HAVE_DUNE_SPGRID
#  include <dune/grid/spgrid.hh>
#  include <dune/grid/spgrid/dgfparser.hh>
#endif

#if HAVE_DUNE_UGGRID || HAVE_UG
#  include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>


// this is used by other headers
typedef Dune::OneDGrid ONED_1D;
typedef Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>> YASP_1D_EQUIDISTANT_OFFSET;
typedef Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>> YASP_2D_EQUIDISTANT_OFFSET;
typedef Dune::YaspGrid<3, Dune::EquidistantOffsetCoordinates<double, 3>> YASP_3D_EQUIDISTANT_OFFSET;
typedef Dune::YaspGrid<4, Dune::EquidistantOffsetCoordinates<double, 4>> YASP_4D_EQUIDISTANT_OFFSET;
#if HAVE_DUNE_ALUGRID
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming> ALU_2D_SIMPLEX_CONFORMING;
typedef Dune::ALUGrid<2, 2, Dune::simplex, Dune::nonconforming> ALU_2D_SIMPLEX_NONCONFORMING;
typedef Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming> ALU_2D_CUBE;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::conforming> ALU_3D_SIMPLEX_CONFORMING;
typedef Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming> ALU_3D_SIMPLEX_NONCONFORMING;
typedef Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming> ALU_3D_CUBE;
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
typedef Dune::UGGrid<2> UG_2D;
typedef Dune::UGGrid<3> UG_3D;
#endif
#if HAVE_ALBERTA
typedef Dune::AlbertaGrid<2, 2> ALBERTA_2D;
typedef Dune::AlbertaGrid<3, 3> ALBERTA_3D;
#endif


namespace Dune {
namespace XT {
namespace Grid {
namespace internal {


// To give better error messages, required below.
template <size_t d>
class ThereIsNoSimplexGridAvailableInDimension
{};


} // namespace internal


/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using Available1dGridTypes = boost::tuple<ONED_1D, YASP_1D_EQUIDISTANT_OFFSET>;

/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using Available2dGridTypes = boost::tuple<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                          ,
                                          ALU_2D_SIMPLEX_CONFORMING,
                                          ALU_2D_SIMPLEX_NONCONFORMING,
                                          ALU_2D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                          ,
                                          UG_2D
#endif
                                          >;

/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 */
using Available3dGridTypes = boost::tuple<YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                          ,
                                          ALU_3D_SIMPLEX_CONFORMING,
                                          ALU_3D_SIMPLEX_NONCONFORMING,
                                          ALU_3D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                          ,
                                          UG_3D
#endif
                                          >;

/**
 * \note Alberta grids are missing here on purpose, these cannot be handled automatically very well.
 * \todo Find a tuple implementation which allows for more than 10 elements!
 */
using AvailableGridTypes = boost::tuple<ONED_1D,
                                        /*YASP_1D_EQUIDISTANT_OFFSET,*/
                                        YASP_2D_EQUIDISTANT_OFFSET,
                                        YASP_3D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                        ,
                                        ALU_2D_SIMPLEX_CONFORMING,
                                        /*ALU_2D_SIMPLEX_NONCONFORMING,*/
                                        /*ALU_2D_CUBE,*/
                                        ALU_3D_SIMPLEX_CONFORMING /*,*/
/*ALU_3D_SIMPLEX_NONCONFORMING,*/
/*ALU_3D_CUBE*/
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                        ,
                                        UG_2D,
                                        UG_3D
#endif
                                        >;


} // namespace Grid
} // namespace XT
} // namespace Dune


using SIMPLEXGRID_1D = ONED_1D;
using SIMPLEXGRID_2D =
#if HAVE_DUNE_ALUGRID
    ALU_2D_SIMPLEX_CONFORMING;
#elif HAVE_DUNE_UGGRID || HAVE_UG
    UG_2D;
#else
    Dune::XT::Grid::internal::ThereIsNoSimplexGridAvailableInDimension<2>;
#endif
using SIMPLEXGRID_3D =
#if HAVE_DUNE_ALUGRID
    ALU_3D_SIMPLEX_CONFORMING;
#elif HAVE_DUNE_UGGRID || HAVE_UG
    UG_3D;
#else
    Dune::XT::Grid::internal::ThereIsNoSimplexGridAvailableInDimension<3>;
#endif


using CUBEGRID_1D = ONED_1D;
using CUBEGRID_2D = YASP_2D_EQUIDISTANT_OFFSET;
using CUBEGRID_3D = YASP_3D_EQUIDISTANT_OFFSET;


#if HAVE_DUNE_ALUGRID || HAVE_DUNE_UGGRID || HAVE_UG
#  define SIMPLEXGRID_2D_AVAILABLE 1
#  define SIMPLEXGRID_3D_AVAILABLE 1
#else
#  define SIMPLEXGRID_2D_AVAILABLE 0
#  define SIMPLEXGRID_3D_AVAILABLE 0
#endif


using GRID_1D = ONED_1D;
using GRID_2D =
#if SIMPLEXGRID_2D_AVAILABLE
    SIMPLEXGRID_2D;
#else
    YASP_2D_EQUIDISTANT_OFFSET;
#endif
using GRID_3D =
#if SIMPLEXGRID_3D_AVAILABLE
    SIMPLEXGRID_3D;
#else
    YASP_3D_EQUIDISTANT_OFFSET;
#endif


#endif // DUNE_XT_GRID_GRIDS_HH
