/**
  \file   examples/grid/provider.cc
  \brief  Demonstrates the capabilities of some Dune::RB::Grid::Providers.
  **/

#include "test_common.hh"

#include <dune/common/mpihelper.hh>
#include <dune/common/parametertree.hh>
#include <dune/common/fvector.hh>

//#if HAVE_DUNE_FEM
//  #include <dune/fem/misc/mpimanager.hh>
//  #include <dune/fem/space/fvspace/fvspace.hh>
//  #include <dune/fem/space/dgspace.hh>
//  #include <dune/fem/gridpart/adaptiveleafgridpart.hh>
//  #include <dune/fem/function/adaptivefunction.hh>
//  #include <dune/fem/operator/projection/l2projection.hh>
//  #include <dune/fem/io/file/datawriter.hh>

//  #include <dune/stuff/fem/customprojection.hh>
//#endif
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/function/expression.hh>

//#if HAVE_EIGEN
#include <Eigen/Core>
//#endif // HAVE_EIGEN

int main(int argc, char** argv)
{
  try {
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::Stuff::Function::Expression<double, 2, double, 1> ScalarFunctionType;
    typename ScalarFunctionType::DomainType scalar_x;
    typename ScalarFunctionType::RangeType scalar_value;
    scalar_x[0] = 1.0;
    scalar_x[1] = 2.0;
    ScalarFunctionType scalar_f_from_single_expression("x", "x[0] + x[1]");
    scalar_f_from_single_expression.report("scalar_f_from_single_expression");
    scalar_f_from_single_expression.evaluate(scalar_x, scalar_value);
    std::cout << "scalar_f_from_single_expression(" << scalar_x << ") = " << scalar_value << " (should be 3)"
              << std::endl;
    std::cout << std::endl;

    std::string scalar_variable = "x";
    std::vector<std::string> scalar_expressions;
    scalar_expressions.push_back("2*x[0] + x[1]*x[1]");
    ScalarFunctionType scalar_f_from_expressions(scalar_variable, scalar_expressions);
    scalar_f_from_expressions.report("scalar_f_from_expressions");
    scalar_f_from_expressions.evaluate(scalar_x, scalar_value);
    std::cout << "scalar_f_from_expressions(" << scalar_x << ") = " << scalar_value << " (should be 6)" << std::endl;
    std::cout << std::endl;

    Dune::ParameterTree scalar_paramTree_single_expression;
    scalar_paramTree_single_expression["variable"]   = "x";
    scalar_paramTree_single_expression["expression"] = "x[0] + x[1]*sin(pi/2)";
    ScalarFunctionType scalar_f_from_paramtree_with_single_expression =
        ScalarFunctionType::createFromParamTree(scalar_paramTree_single_expression);
    scalar_f_from_paramtree_with_single_expression.report("scalar_f_from_paramtree_with_single_expression");
    scalar_f_from_paramtree_with_single_expression.evaluate(scalar_x, scalar_value);
    std::cout << "scalar_f_from_paramtree_with_single_expression(" << scalar_x << ") = " << scalar_value
              << " (should be 3)" << std::endl;
    std::cout << std::endl;

    Dune::ParameterTree scalar_paramTree_expressions;
    scalar_paramTree_expressions["variable"]   = "x";
    scalar_paramTree_expressions["expression"] = "[x[0] + x[1]*sin(pi/2)]";
    ScalarFunctionType scalar_f_from_paramtree_with_expressions =
        ScalarFunctionType::createFromParamTree(scalar_paramTree_expressions);
    scalar_f_from_paramtree_with_expressions.report("scalar_f_from_paramtree_with_expressions");
    scalar_f_from_paramtree_with_expressions.evaluate(scalar_x, scalar_value);
    std::cout << "scalar_f_from_paramtree_with_expressions(" << scalar_x << ") = " << scalar_value << " (should be 3)"
              << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    typedef Dune::Stuff::Function::Expression<double, 2, double, 2> VectorFunctionType;
    typename VectorFunctionType::DomainType vector_x;
    typename VectorFunctionType::RangeType vector_value;
    vector_x[0]                 = 1.0;
    vector_x[1]                 = 2.0;
    std::string vector_variable = "x";
    std::vector<std::string> vector_expressions;
    vector_expressions.push_back("2*x[0] + x[1]*x[1]");
    vector_expressions.push_back("-2*x[0] - x[1]*x[1]");
    VectorFunctionType vector_f_from_expressions(vector_variable, vector_expressions);
    vector_f_from_expressions.report("vector_f_from_expressions");
    vector_f_from_expressions.evaluate(vector_x, vector_value);
    std::cout << "vector_f_from_expressions(" << vector_x << ") = " << vector_value << " (should be 6 -6)" << std::endl;
    std::cout << std::endl;

    Dune::ParameterTree vector_paramTree_expressions;
    vector_paramTree_expressions["variable"]   = "x";
    vector_paramTree_expressions["expression"] = "[x[0] + x[1]*sin(pi/2); -1*(x[0] + x[1]*sin(pi/2))]";
    VectorFunctionType vector_f_from_paramtree_with_expressions =
        VectorFunctionType::createFromParamTree(vector_paramTree_expressions);
    vector_f_from_paramtree_with_expressions.report("vector_f_from_paramtree_with_expressions");
    vector_f_from_paramtree_with_expressions.evaluate(vector_x, vector_value);
    std::cout << "vector_f_from_paramtree_with_expressions(" << vector_x << ") = " << vector_value
              << " (should be 3 -3)" << std::endl;

    //#if HAVE_EIGEN
    std::cout << std::endl;
    typedef Dune::Stuff::Function::Expression<double, 50, double, 1> ParameterFunctionType;
    typedef Eigen::VectorXd ParameterType;
    ParameterType scalarParam(1);
    scalarParam(0) = 1.0;
    ParameterType vectorParam(5);
    for (unsigned int i = 0; i < 5; ++i)
      vectorParam(i) = double(i + 2);
    ParameterType tooLargeParam(51);
    for (unsigned int i = 0; i < 51; ++i)
      tooLargeParam(i) = double(i);
    ParameterType paramValue;
    std::cout << std::endl;

    ParameterFunctionType single_paramFunction("mu", "mu[0]");
    single_paramFunction.report("single_paramFunction");
    single_paramFunction.evaluate(scalarParam, paramValue);
    std::cout << "single_paramFunction(" << scalarParam << ") = " << paramValue << " (should be 1)" << std::endl;
    single_paramFunction.evaluate(vectorParam, paramValue);
    std::cout << "single_paramFunction(" << vectorParam << ") = " << paramValue << " (should be 2)" << std::endl;
    // this will fail:
    //    scalar_paramFunction.evaluate(tooLargeParam, paramValue);
    //    std::cout << "scalar_paramFunction(" << tooLargeParam << ") = " << paramValue << " (should fail)"<< std::endl;
    std::cout << std::endl;

    ParameterFunctionType multiple_paramFunction("mu", "mu[0] + mu[1] + mu[2] + mu[3] + mu[4]");
    multiple_paramFunction.report("multiple_paramFunction");
    multiple_paramFunction.evaluate(scalarParam, paramValue);
    std::cout << "multiple_paramFunction(" << scalarParam << ") = " << paramValue
              << " (should be 1, mu[1] .. mu[4] are treated as 0)" << std::endl;
    multiple_paramFunction.evaluate(vectorParam, paramValue);
    std::cout << "multiple_paramFunction(" << vectorParam << ") = " << paramValue << " (should be 20)" << std::endl;
    // this will fail:
    //    multiple_paramFunction.evaluate(tooLargeParam, paramValue);
    //    std::cout << "scalar_paramFunction(" << tooLargeParam << ") = " << paramValue << " ()"<< std::endl;
    //#endif // HAVE_EIGEN

    //#if HAVE_DUNE_FEM
    //    Dune::GridPtr< GridType > gridPtr( paramTree.get("dgf_file", "dummy") );
    //    typedef Dune::AdaptiveLeafGridPart< GridType >
    //        GridPartType;
    //    GridPartType gridPart_(*gridPtr);
    //    typedef ParameterFunctionType::FunctionSpaceType FSpace;
    //    typedef Dune::DiscontinuousGalerkinSpace< FSpace,
    //                                              GridPartType,
    //                                              1 >
    //    DiscreteFunctionSpaceType;
    //    typedef Dune::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
    //    DiscreteFunctionType;
    //    DiscreteFunctionSpaceType disc_space(gridPart_);
    //    DiscreteFunctionType rf_disc("rf", disc_space);
    //    typedef Dune::tuple< const DiscreteFunctionType* >
    //    OutputTupleType;
    //    typedef Dune::DataWriter< GridPartType::GridType,
    //                              OutputTupleType >
    //    DataWriterType;
    //    Dune::Stuff::Fem::BetterL2Projection::project(parameterFunction, rf_disc);
    //    OutputTupleType out(&rf_disc);
    //    DataWriterType dt(gridPart_.grid(), out);
    //    dt.write();
    //#endif
  } catch (Dune::Exception& e) {
    std::cout << e.what() << std::endl;
    return 1;
  }
  return 0;
}
