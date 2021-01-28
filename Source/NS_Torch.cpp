
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>

#include <torch/torch.h>
#include <torch/script.h> // One-stop header.


using namespace amrex;


void
NavierStokesBase::training_Unet_2d()

{

  std::cout << "\n WE ARE IN training_Unet_2d ROUTINE \n";

////////////////////////////////////////////////////////
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

    const auto bx_onegrid_bounds = amrex::ubound(bx_onegrid);

/////////////////////////////////////////////////////////////

  torch::Tensor t1 = torch::zeros({1,1,bx_onegrid_bounds.x+1,bx_onegrid_bounds.y+1});
//  std::cout << "Tensor from array:\n" << t1 << '\n';

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

      t1[0][0][j][i] = Snew_vel(i,j,k);  // Inverting i and j index because of the PyTorch formalism
//      amrex::Print() << "i= " << i << " j= " << j << "   array= "  << Snew_vel(i,j,k) << std::endl;
//      amrex::Print() << "i= " << i << " j= " << j << "   array Torch= "  << t1[i][j] << std::endl;
    });

  }

//std::cout << "Tensor from array:\n" << t1 << '\n';



amrex::Print() << "DEBUG WE LEAVE training_Unet_2d ROUTINE  \n ";

}




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

// We want to launch Torch on just 1 mpi process
//if (ParallelDescriptor::IOProcessor()){


  torch::Tensor t1 = torch::zeros({1,4,bx_onegrid_bounds.x+1,bx_onegrid_bounds.y+1}); 
//  std::cout << "Tensor from array:\n" << t1 << '\n';

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

      t1[0][0][j][i] = Snew_vel(i,j,k);  // Inverting i and j index because of the PyTorch formalism
      t1[0][1][j][i] = Snew_vel(i,j,k);
      t1[0][2][j][i] = Snew_vel(i,j,k);
      t1[0][3][j][i] = Snew_vel(i,j,k); 
//      amrex::Print() << "i= " << i << " j= " << j << "   array= "  << Snew_vel(i,j,k) << std::endl;
//      amrex::Print() << "i= " << i << " j= " << j << "   array Torch= "  << t1[i][j] << std::endl;
    });

  }

//std::cout << "Tensor from array:\n" << t1 << '\n';



torch::jit::script::Module module;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    module = torch::jit::load("traced_unet_model.pt");
  }
  catch (const c10::Error& e) {
    amrex::Abort( "error loading the model\n");
  }


    std::vector<torch::jit::IValue> inputs{t1};

//   std::vector<torch::jit::IValue> inputs;
//    inputs.push_back(torch::rand({1, 4, 512, 512}));


 at::Tensor output = module.forward(inputs).toTensor();

//} // End of ParallelDescriptor

amrex::Print() << "DEBUG WE ARE AFTER TORCH \n ";

//ParallelDescriptor::Barrier();

}

