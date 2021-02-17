
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>

#include <torch/torch.h>
#include <torch/script.h> // One-stop header.


using namespace amrex;



void
NavierStokesBase::test_libtorch()

{

  std::cout << "\n WE ARE IN TEST_LIBTORCH ROUTINE \n";


// First we convert a multi-box MultiFab to a unique box 

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


  // For debugging
  // VisMF::Write(Snew,"Snew");
  // VisMF::Write(Snew_onegrid,"Snew_onegrid");


// Now we convert the one grid MultiFab to a Torch tensor
  
  // amrex::Print() << "\n DEBUG amrex::ubound(bx) = " << amrex::ubound(bx_onegrid) << std::endl;
  // amrex::Print() << "\n DEBUG amrex::lbound(bx) = " << amrex::lbound(bx_onegrid) << std::endl;

  const auto bx_onegrid_bounds = amrex::ubound(bx_onegrid);

  // amrex::Print() << "\n DEBUG box_bounds.x = " << bx_onegrid_bounds.x << std::endl;

  torch::Device device(torch::kCPU);
  if (torch::cuda::is_available()) {
    std::cout << "CUDA is available! Training on GPU." << std::endl;
  torch::Device   device = torch::Device(torch::kCUDA);
  }
 
  // Here we put 4 components because the trained model has 4
  torch::Tensor t1 = torch::zeros({1,4,bx_onegrid_bounds.x+1,bx_onegrid_bounds.y+1},device); 

  //  std::cout << "Tensor from array:\n" << t1 << '\n';

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(Snew_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& Snew_vel = Snew_onegrid.array(mfi,Xvel);
    FArrayBox  fab_tmp(bx,1);
    auto       fab  = fab_tmp.array(0);

    // Here we retrieve data from device to a fab array on host    
    amrex::ParallelFor(bx, 1, [fab,Snew_vel]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
    fab(i,j,k)=Snew_vel(i,j,k) ;
    });

    // Now we copy the fab array into the Torch Tensor
    {
    for (int jj = 0; jj < bx_onegrid_bounds.y+1; jj++ )
    {
      for (int ii = 0; ii < bx_onegrid_bounds.x+1; ii++ )
      {
      //  std::cout << "\n DEBUG i,j=" << ii << " " << jj << " fab= " << fab(ii,jj,0,0) << std::endl;
      t1[0][0][ii][jj] = fab(ii,jj,0,0);
      t1[0][1][ii][jj] = fab(ii,jj,0,0);
      t1[0][2][ii][jj] = fab(ii,jj,0,0);
      t1[0][3][ii][jj] = fab(ii,jj,0,0);
      }
    }

    amrex::Print() << "\n DEBUG SNEW " << Snew_onegrid[mfi] << "\n";

    }
  }

  std::cout << "Tensor from array:\n" << t1 << '\n';


// Loading the model with TorchScript
  torch::jit::script::Module module;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    module = torch::jit::load("traced_unet_model.pt");
  }
  catch (const c10::Error& e) {
    amrex::Abort( "error loading the model\n");
  }

// Passing the tensor to the model
  std::vector<torch::jit::IValue> inputs{t1};

  // Below is a test with a random tensor 
  //  std::vector<torch::jit::IValue> inputs;
  //  inputs.push_back(torch::rand({1, 4, 512, 512}));

  at::Tensor output = module.forward(inputs).toTensor();


// If we ended up here, it's all good !
  amrex::Print() << "DEBUG WE ARE AFTER TORCH \n ";

}

