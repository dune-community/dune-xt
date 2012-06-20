#include "cmake_config.h"

#include <cstdio>
#if defined(USE_PARDG_ODE_SOLVER) && defined(USE_BFG_CG_SCHEME)
#warning("USE_PARDG_ODE_SOLVER enabled, might conflict with custom solvers")
#endif

// the adaption manager might be troublesome with certain gridparts/spaces, so we needed a easy way to disable it
#ifndef ENABLE_ADAPTIVE
#define ENABLE_ADAPTIVE 1
#endif

#if defined(UGGRID) && defined(DEBUG)
#warning("UGGRID in debug mode is likely to produce a segfault")
#endif

#if !defined(POLORDER)
#define POLORDER 0
#warning("using default polorder 0 for all spaces")
#endif // if !defined (POLORDER)

#if !defined(PRESSURE_POLORDER)
#define PRESSURE_POLORDER POLORDER
#endif

#if !defined(VELOCITY_POLORDER)
#define VELOCITY_POLORDER POLORDER
#endif

#if !defined(TESTCASE)
#define TESTCASE TestCase3D
#endif

#define TESTCASE_NAME "TESTCASE"

#if ((defined(SGRID) || defined(ALUGRID_SIMPLEX) || defined(ALUGRID_CUBE)) && (GRIDDIM == 3)) || defined(UGGRID)       \
    || defined(YASPGRID)
// this is no mistake, ALU is indeed only incompatible in 3d
#define OLD_DUNE_GRID_VERSION
#endif // if ( ( defined (SGRID) || defined (ALUGRID_SIMPLEX) || defined (ALUGRID_CUBE ) ) && (GRIDDIM == 3 ) ) ||
// defined (UGGRID) || defined (YASPGRID)

#define USE_GRPAE_VISUALISATION (HAVE_GRAPE && !defined(AORTA_PROBLEM))

#include <vector>
#include <string>

#include <iostream>
#include <cmath>
#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

// !ATTENTION: undef's GRIDDIM

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh> // for the grid

#include <dune/fem/solver/oemsolver/oemsolver.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/dgspace/dgadaptiveleafgridpart.hh>
#include <dune/fem/function/adaptivefunction.hh> // for AdaptiveDiscreteFunction
#include <dune/fem/operator/projection/l2projection.hh>

#include <dune/stuff/printing.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/logging.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/profiler.hh>
#include <dune/stuff/datawriter.hh>
#include <dune/stuff/customprojection.hh>
#include <dune/common/collectivecommunication.hh>

#include <stdio.h>
#include <dune/fem/misc/mpimanager.hh> // An initializer of MPI
#include <dune/fem/space/fvspace/fvspace.hh>
#include <dune/stuff/parametercontainer.hh>
#include <dune/stuff/misc.hh>
#include <dune/stuff/runtimefunction.hh>
#include <dune/stuff/expressions/mathexpr.h>
// do whatever you like to this file to test out simple and small stuff

int main(int argc, char** argv)
{
  Dune::MPIManager::initialize(argc, argv);
  if (!(Parameters().ReadCommandLine(argc, argv))) {
    return 1;
  }
  const bool useLogger = true;
  Logger().Create(Parameters().getParam("loglevel", 62, useLogger),
                  Parameters().getParam("logfile", std::string("dune_stokes"), useLogger),
                  Parameters().getParam("fem.io.logdir", std::string(), useLogger));

  typedef Dune::FunctionSpace<double, double, 2, 2> FSpace;
  FSpace fSpace;
  typedef Stuff::RuntimeFunction<FSpace> RF;
  RF rf("func", fSpace);

  /* ********************************************************************** *
    * initialize the grid                                                    *
    * ********************************************************************** */
  Logging::LogStream& infoStream  = Logger().Info();
  Logging::LogStream& debugStream = Logger().Dbg();
  infoStream << "\n- initialising grid" << std::endl;
  const int gridDim = GridType::dimensionworld;
  Dune::GridPtr<GridType> gridPtr(Parameters().DgfFilename(gridDim));
  int refine_level = (Parameters().getParam("minref", 0)) * Dune::DGFGridInfo<GridType>::refineStepsForHalf();
  gridPtr->globalRefine(refine_level);

  typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;
  GridPartType gridPart_(*gridPtr);

  typedef Dune::DiscontinuousGalerkinSpace<FSpace, GridPartType, 1> DiscreteFunctionSpaceType;
  typedef Dune::AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
  DiscreteFunctionSpaceType disc_space(gridPart_);
  DiscreteFunctionType rf_disc("rf", disc_space);
  typedef Dune::Tuple<const DiscreteFunctionType*> OutputTupleType;
  typedef Dune::DataWriter<GridPartType::GridType, OutputTupleType> DataWriterType;
  Dune::L2Projection<double, double, RF, DiscreteFunctionType>()(rf, rf_disc);
  OutputTupleType out(&rf_disc);
  DataWriterType dt(gridPart_.grid(), "dummy", out, 0, 0);
  dt.write(0.0, 0);

  return 0;
} // main
