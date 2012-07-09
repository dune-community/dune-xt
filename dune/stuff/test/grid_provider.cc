/**
  \file   examples/grid/provider.cc
  \brief  Demonstrates the capabilities of some Dune::RB::Grid::Providers.
  **/

#include "config.h"

// system
#include <iostream>
#include <fstream>

// boost
#include <boost/filesystem.hpp>

// dune-common
#include <dune/common/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>
#include <dune/common/timer.hh>

// dune-rb
#include <dune/rb/grid/provider/cornerpoint.hh>
#include <dune/rb/grid/provider/cube.hh>

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
    file << "[rb.grid.provider.cube]" << std::endl;
    file << "level = 2" << std::endl;
    file << "visualize.grid = rb_grid_provider_cube_grid" << std::endl;
    file << "visualize.msGrid = rb_grid_provider_cube_msGrid" << std::endl;
    file << std::endl;
    file << "[rb.grid.provider.cornerpoint]" << std::endl;
    file << "filename = "
            "/home/felix/Projects/dune/dune-lrbms/dune-rb/dune/rb/grid/examples/data/johansen_formation.grdecl"
         << std::endl;
    file << "visualize.grid = rb_grid_provider_cornerpoint_grid" << std::endl;
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
    Dune::ParameterTreeParser::readINITree("provider.param", paramTree);
  } else if (argc == 2) {
    Dune::ParameterTreeParser::readINITree(argv[1], paramTree);
  } else {
    Dune::ParameterTreeParser::readOptions(argc, argv, paramTree);
  }
  if (paramTree.hasKey("paramfile")) {
    Dune::ParameterTreeParser::readINITree(paramTree.get<std::string>("paramfile"), paramTree, false);
  }
}

template <class GridViewType>
int walkGrid(GridViewType& gridView)
{
  int numElements = 0;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::template Codim<0>::Iterator IteratorType;
  for (IteratorType iterator = gridView.template begin<0>(); iterator != gridView.template end<0>(); ++iterator) {
    const EntityType& entity = *iterator;
    ++numElements;
  }
  return numElements;
} // void walkGrid(const GridType& grid)

/*template< class MDGridType >
int walkMDGrid(MDGridType& mdGrid)
{
  int numElements = 0;
  // loop over all subdomains
  for (typename MDGridType::SubDomainIndexType subdomainIndex = 0;
       subdomainIndex <= mdGrid.maxSubDomainIndex;
       ++subdomainIndex) {
    // subdomain grid
    typedef typename MDGridType::SubDomainGrid SDGridType;
    const SDGridType& sdGrid = mdGrid.subDomain(subdomainIndex);
    // subdomain grid view
    typedef typename SDGridType::LeafGridView SDGridViewType;
    const SDGridViewType& sdGridView = sdGrid.leafView();
    // walk the subdomain grid
    typedef typename SDGridViewType::template Codim<0>::Iterator SDElementIteratorType;
    typedef typename SDGridViewType::template Codim<0>::Entity SDElementType;
    for (SDElementIteratorType sdElementIterator = sdGridView.template begin<0>();
         sdElementIterator != sdGridView.template end<0>();
         ++sdElementIterator) {
      const SDElementType& sdElement = *sdElementIterator;
      ++numElements;
    } // walk the subdomain grid
  } // loop over all subdomains
  return numElements;
} // void walkMDGrid(const GridType& grid)*/

template <class GridProviderType>
void measureTiming(GridProviderType& gridProvider)
{
  Dune::Timer timer;
  typename GridProviderType::GridType::LeafGridView gridView = gridProvider.grid().leafView();
  const int numHostGridElements = walkGrid(gridView);
  std::cout << "  host grid:        " << timer.elapsed() << " sec, " << numHostGridElements << " elements" << std::endl;
  timer.reset();
  typename GridProviderType::MsGridType::GridViewType wrappedGridView = gridProvider.msGrid().gridView();
  const int numWrappedGridElements = walkGrid(wrappedGridView);
  std::cout << "  wrapped grid:     " << timer.elapsed() << " sec, " << numWrappedGridElements << " elements"
            << std::endl;
  timer.reset();
  typename GridProviderType::MsGridType::LocalGridViewType localGridView = gridProvider.msGrid().localGridView(0);
  const int numLocalGridElements = walkGrid(localGridView);
  std::cout << "  local grid:       " << timer.elapsed() << " sec, " << numLocalGridElements << " elements"
            << std::endl;
  timer.reset();
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
    ensureParamFile("provider.param");
    Dune::ParameterTree paramTree;
    initParamTree(argc, argv, paramTree);
    // unitcube
    typedef Dune::RB::Grid::Provider::Cube<Dune::GridSelector::GridType,
                                           Dune::RB::Grid::Multiscale::Multidomain<Dune::GridSelector::GridType>>
        CubeProviderType;
    CubeProviderType cubeProvider(paramTree);
// cornerpoint
#ifdef HAVE_DUNE_CORNERPOINT
    typedef Dune::RB::Grid::Provider::Cornerpoint CornerpointGridProviderType;
    CornerpointGridProviderType cornerpointGridProvider(paramTree);
#endif // HAVE_DUNE_CORNERPOINT
    // visualize
    cubeProvider.visualize(paramTree);
#ifdef HAVE_DUNE_CORNERPOINT
    cornerpointGridProvider.visualize(paramTree);
#endif
    // measure timing
    std::cout << std::endl;
    measureTiming(cubeProvider);
  } catch (Dune::Exception& e) {
    std::cout << e.what() << std::endl;
  }
  return 0;
}
