
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>

#include <torch/torch.h>

using namespace amrex;


void
NavierStokesBase::test_libtorch()

{

std::cout << "\n WE ARE IN TEST_LIBTORCH ROUTINE \n";


  torch::Tensor tensor = torch::eye(3);
  std::cout << tensor << std::endl;


}

