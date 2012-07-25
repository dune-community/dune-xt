#include "test_common.hh"

#include <dune/stuff/grid/information.hh>

int main(int argc, char** argv)
{


  try {
    // mpi
    Dune::MPIHelper::instance(argc, argv);


  } catch (Dune::Exception& e) {
    std::cout << e.what() << std::endl;
    return 1;
  }
  return 0;
}
