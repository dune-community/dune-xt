#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include <config.h>
#endif // ifdef HAVE_CMAKE_CONFIG

// system
#include <iostream>
#include <fstream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/timer.hh>

#include <dune/stuff/common/parameter/tree.hh>
//#include <dune/stuff/grid/provider/cornerpoint.hh>
#include <dune/stuff/grid/provider/cube.hh>

using namespace Dune::Stuff;

const std::string id = "grid_provider";

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
    file << "[stuff.grid.provider.cube]" << std::endl;
    file << "level = 2" << std::endl;
    file << std::endl;
    file << "[stuff.grid.provider.cornerpoint]" << std::endl;
    file << "filename = /dune-stuff/data/grid/johansen_formation.grdecl" << std::endl; // has to be an absolute path
    file.close();
  } // only write param file if there is none
} // void ensureParamFile()

template <class GridViewType>
int walkGridView(const GridViewType& gridView)
{
  int numElements = 0;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::template Codim<0>::Iterator IteratorType;
  for (IteratorType iterator = gridView.template begin<0>(); iterator != gridView.template end<0>(); ++iterator) {
    const EntityType& entity = *iterator;
    ++numElements;
  }
  return numElements;
} // int walkGridView(GridViewType& gridView)

/**
  \brief  Main routine.
  **/
int main(int argc, char** argv)
{
  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);
    // parameter
    ensureParamFile(id + ".param");
    Dune::ParameterTree paramTree = Common::Parameter::Tree::init(argc, argv, id + ".param");
    // timer
    Dune::Timer timer;
    // unitcube
    std::cout << "creating Cube... ";
    typedef Grid::Provider::Cube<Dune::GridSelector::GridType> CubeProviderType;
    CubeProviderType cubeProvider(paramTree);
    std::cout << "done (took " << timer.elapsed() << "s)" << std::endl;
    std::cout << "timing... ";
    timer.reset();
    const int numElements = walkGridView(cubeProvider.grid().leafView());
    std::cout << "done (took " << timer.elapsed() << "s, has " << numElements << " elements)" << std::endl;
    std::cout << "visualizing... ";
    timer.reset();
    cubeProvider.visualize(id);
    std::cout << "done (took " << timer.elapsed() << "s, see " << id;
    if (CubeProviderType::dim == 1)
      std::cout << ".vtp";
    else
      std::cout << ".vtu";
    std::cout << ")" << std::endl;
    //    // cornerpoint
    //#ifdef HAVE_DUNE_CORNERPOINT
    //    typedef Grid::Provider::Cornerpoint CornerpointGridProviderType;
    //    CornerpointGridProviderType cornerpointGridProvider(paramTree);
    //    cornerpointGridProvider.visualize(paramTree);
    //#endif
  } catch (Dune::Exception& e) {
    std::cout << e.what() << std::endl;
    return 1;
  }
  return 0;
}
