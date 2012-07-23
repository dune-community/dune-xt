/**
  \file   examples/grid/provider.cc
  \brief  Demonstrates the capabilities of some Dune::RB::Grid::Providers.
  **/

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/fvspace/fvspace.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/stuff/fem/customprojection.hh>
#endif

// system
#include <iostream>
#include <fstream>

// boost
#include <boost/filesystem.hpp>


#include <dune/common/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/function/expression.hh>

/**
  \brief      Creates a parameter file if it does not exist.

              Nothing is done if the file already exists. If not, a parameter file will be created with all neccessary
              keys and values.
  \param[in]  filename
              (Relative) path to the file.
  **/
void ensureParamFile(std::string filename)
{
  // only write param file if there is none
  if (!boost::filesystem::exists(filename)) {
    std::ofstream file;
    file.open(filename);
    file << "[rb.examples.parameter.function]" << std::endl;
    file << "variable = mu" << std::endl;
    file << "expression.0 = sin(mu[0]) + 2" << std::endl;
    file << "expression.1 = mu[0]^2 + mu[1]" << std::endl;
    file << "expression.2 = 0.5 * mu[2]" << std::endl;
    file << std::endl;
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

/**
  \brief      Fills a Dune::ParameterTree given a parameter file or command line arguments.
  \param[in]  argc
              From \c main()
  \param[in]  argv
              From \c main()
  \param[out] paramTree
              The Dune::ParameterTree that is to be filled.
  **/
void initParamTree(int argc, char** argv, Dune::ParameterTree& paramTree)
{
  if (argc == 1) {
    Dune::ParameterTreeParser::readINITree("function.param", paramTree);
  } else if (argc == 2) {
    Dune::ParameterTreeParser::readINITree(argv[1], paramTree);
  } else {
    Dune::ParameterTreeParser::readOptions(argc, argv, paramTree);
  }
  if (paramTree.hasKey("paramfile")) {
    Dune::ParameterTreeParser::readINITree(paramTree.get<std::string>("paramfile"), paramTree, false);
  }
}

/**
  \brief  Main routine.
  **/
int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);
    // parameter
    ensureParamFile("function.param");
    Dune::ParameterTree paramTree;
    initParamTree(argc, argv, paramTree);
    typedef Dune::GridSelector::GridType GridType;
    // parameter domain
    const int dimDomain = GridType::dimensionworld;
    typedef Dune::FieldVector<double, dimDomain> DomainType;
    DomainType mu(1.0);
    const int dimRange = 2;
    typedef Dune::FieldVector<double, dimRange> RangeType;
    RangeType result(0.0);
    // parameter function
    std::cout << "creating Dune::RB::Paramter::Function< double, " << dimDomain << ", double, " << dimRange << " >... ";
    typedef Dune::Stuff::Function::Expression<double, dimDomain, double, dimRange> ParameterFunctionType;
    ParameterFunctionType parameterFunction(paramTree.sub("rb.examples.parameter.function"));
    std::cout << "done" << std::endl;
    std::cout << "variable is:     " << parameterFunction.variable() << std::endl;
    std::cout << "expressions are: " << parameterFunction.expression(0) << std::endl;
    for (int i = 1; i < dimRange; ++i) {
      std::cout << "                 " << parameterFunction.expression(i) << std::endl;
    }
    parameterFunction.evaluate(mu, result);
    std::cout << "evaluate(" << mu << ") = " << result << std::endl;
    mu = 0.5;
    parameterFunction.evaluate(mu, result);
    std::cout << "evaluate(" << mu << ") = " << result << std::endl;

#if HAVE_DUNE_FEM
    Dune::GridPtr<GridType> gridPtr(paramTree.get("dgf_file", "dummy"));
    typedef Dune::AdaptiveLeafGridPart<GridType> GridPartType;
    GridPartType gridPart_(*gridPtr);

    typedef ParameterFunctionType::FunctionSpaceType FSpace;
    typedef Dune::DiscontinuousGalerkinSpace<FSpace, GridPartType, 1> DiscreteFunctionSpaceType;
    typedef Dune::AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DiscreteFunctionType;
    DiscreteFunctionSpaceType disc_space(gridPart_);
    DiscreteFunctionType rf_disc("rf", disc_space);
    typedef Dune::tuple<const DiscreteFunctionType*> OutputTupleType;
    typedef Dune::DataWriter<GridPartType::GridType, OutputTupleType> DataWriterType;
    Dune::Stuff::Fem::BetterL2Projection::project(parameterFunction, rf_disc);
    OutputTupleType out(&rf_disc);
    DataWriterType dt(gridPart_.grid(), out);
    dt.write();
#endif
  } catch (Dune::Exception& e) {
    std::cout << e.what() << std::endl;
    return 1;
  }
  return 0;
}
