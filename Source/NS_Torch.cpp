
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>

#include <torch/torch.h>

using namespace amrex;


void
NavierStokesBase::test_libtorch()

{

  std::cout << "\n WE ARE IN TEST_LIBTORCH ROUTINE \n";


  MultiFab& Snew   = get_old_data(State_Type);
  torch::Tensor t1 = torch::zeros({8,8});
  std::cout << "Tensor from array:\n" << t1 << '\n';

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(Snew,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& Snew_vel = Snew.array(mfi,Xvel);

    amrex::ParallelFor(bx, [t1, Snew_vel]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {

      t1[j][i] = Snew_vel(i,j,k);  // Inverting i and j index because of the PyTorch formalism
      amrex::Print() << "i= " << i << " j= " << j << "   array= "  << Snew_vel(i,j,k) << std::endl;
      amrex::Print() << "i= " << i << " j= " << j << "   array Torch= "  << t1[i][j] << std::endl;
    });

  }

std::cout << "Tensor from array:\n" << t1 << '\n';



}

