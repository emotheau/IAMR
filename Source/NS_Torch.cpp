
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>

#include <torch/torch.h>

using namespace amrex;


void
NavierStokesBase::test_libtorch()

{

  std::cout << "\n WE ARE IN TEST_LIBTORCH ROUTINE \n";


// Getting new state data
  MultiFab& Snew   = get_new_data(State_Type);


    // get boxArray
    BoxArray ba = Snew.boxArray();
    // minimalBox() computes a single box to enclose all the boxes
    // enclosedCells() converts it to a cell-centered Box
    Box bx_onegrid = ba.minimalBox().enclosedCells();

    // number of cells in the coarse domain
    Print() << "npts in coarse domain = " << bx_onegrid.numPts() << std::endl;

    // BoxArray, DistributionMapping, and MultiFab with one grid
    BoxArray ba_onegrid(bx_onegrid);
    DistributionMapping dmap_onegrid(ba_onegrid);
    MultiFab Snew_onegrid(ba_onegrid,dmap_onegrid,Snew.nComp(),0);

    // copy data into MultiFab with one grid
    Snew_onegrid.ParallelCopy(Snew,0,0,Snew.nComp(),0,0);


VisMF::Write(Snew,"Snew");
VisMF::Write(Snew_onegrid,"Snew_onegrid");

amrex::Print() << "\n DEBUG amrex::ubound(bx) = " << amrex::ubound(bx_onegrid) << std::endl;
amrex::Print() << "\n DEBUG amrex::lbound(bx) = " << amrex::lbound(bx_onegrid) << std::endl;

const auto bx_onegrid_bounds = amrex::ubound(bx_onegrid);

amrex::Print() << "\n DEBUG box_bounds.x = " << bx_onegrid_bounds.x << std::endl;




  torch::Tensor t1 = torch::zeros({bx_onegrid_bounds.x+1,bx_onegrid_bounds.y+1}); 
  std::cout << "Tensor from array:\n" << t1 << '\n';

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(Snew_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& Snew_vel = Snew_onegrid.array(mfi,Xvel);

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

