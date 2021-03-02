
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <torch/torch.h>
#include <torch/script.h> // One-stop header.
#include <iostream>
#include <memory>

using namespace amrex;



void
NavierStokesBase::test_libtorch()

{

std::cout << "\n WE ARE IN TEST_LIBTORCH ROUTINE \n";

// First we convert a multi-box MultiFab to a unique box 

  // Getting new state data
  MultiFab& Snew   = get_new_data(State_Type);
//  Snew.FillBoundary( geom.periodicity());
  const Real cur_time = state[State_Type].curTime();
  FillPatch(*this,Snew,Snew.nGrow(),cur_time,State_Type,0,NUM_STATE);

/*
    for (MFIter mfi(Snew,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    amrex::Print() << Snew[mfi];
  }
*/

  // get boxArray
  BoxArray ba = Snew.boxArray();
  // minimalBox() computes a single box to enclose all the boxes
  // enclosedCells() converts it to a cell-centered Box
  Box bx_onegrid = ba.minimalBox().enclosedCells();

  // number of cells in the coarse domain
//  Print() << "npts in coarse domain = " << bx_onegrid.numPts() << std::endl;

  // Here we set that we have only 2 components
  // FIX ME : WE HAVE TO DO SOMETHING MORE GENERIC LATER
  int ncomp = 2;

  // BoxArray, DistributionMapping, and MultiFab with one grid
  BoxArray ba_onegrid(bx_onegrid);
  DistributionMapping dmap_onegrid(ba_onegrid);
  MultiFab Snew_onegrid(ba_onegrid,dmap_onegrid,ncomp,1);

  // copy data into MultiFab with one grid
  // FIX ME : SET THAT WE START FROM X_VEL
  Snew_onegrid.ParallelCopy(Snew,0,0,ncomp,1,1);
//    FillPatch(*this,Snew,Snew.nGrow(),cur_time,State_Type,0,NUM_STATE);

  /*
    for (MFIter mfi(Snew_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    amrex::Print() << Snew_onegrid[mfi];
  }
  */

  // For debugging
  VisMF::Write(Snew,"Snew_coarse_b4_ML");
  VisMF::Write(Snew_onegrid,"Snew_onegrid_b4_ML");


// Now we refine Snew_onegrid

  // create a copy of the boxarray and refine it. Use it to create a new Multifab to hold refined data
  // FIX ME : THIS SHOULD BE SET BY USER IN INPUT FILE
  int ratio = 2;

  BoxArray baRefined = ba_onegrid;
  baRefined.refine(ratio);

  DistributionMapping dmapRefined(baRefined);

  // Note that Snew_onegrid should have 0 ghost-cells because the target Torch Tensor will have 0
  MultiFab SnewRefined(baRefined, dmapRefined, ncomp, 0);

  // get geometry object from global scope and create a refined copy
  Box bx_fine_geom = geom.Domain();
  bx_fine_geom.refine(ratio);

  Geometry fine_geom(bx_fine_geom, geom.ProbDomain(), geom.Coord(), geom.isPeriodic());


//amrex::Print() << "\n DEBUG bx_fine_geom = " << bx_fine_geom << std::endl;

//  for (MFIter mfi(Snew_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
//  {
//    amrex::Print() << Snew_onegrid[mfi];
//  }


// Interpolating Snew_onegrid on SnewRefined

  // instantiate interpolator
  Interpolater*  interpolater = &cell_cons_interp;

  // iterate through boxes in multifab and interpolate to refined multifab
  IntVect new_ratio(AMREX_D_DECL(ratio,ratio,ratio));
  const IntVect& rr = new_ratio;

  for (MFIter mfi(SnewRefined); mfi.isValid(); ++mfi)
  {
    FArrayBox& ffab = (SnewRefined)[mfi];
    const FArrayBox& cfab = (Snew_onegrid)[mfi];
    const Box&  bx   = mfi.tilebox();

//    amrex::Print() << "\n DEBUG bx = " << bx << std::endl;


    Vector<BCRec> bx_bcrec(ncomp);

    interpolater->interp(cfab,0,ffab,0,ncomp,bx,rr,
                         geom,fine_geom,bx_bcrec,0,0,RunOn::Host);
  }

  // For debugging
  VisMF::Write(SnewRefined,"SnewRefined_b4_ML");

  /*
    for (MFIter mfi(SnewRefined,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    amrex::Print() << SnewRefined[mfi];
  }
*/

// Now we convert the one grid MultiFab to a Torch tensor

//  amrex::Print() << "\n DEBUG amrex::ubound(bx_fine_geom) = " << amrex::ubound(bx_fine_geom) << std::endl;
//  amrex::Print() << "\n DEBUG amrex::lbound(bx_fine_geom) = " << amrex::lbound(bx_fine_geom) << std::endl;

  const auto bx_onegrid_bounds = amrex::ubound(bx_fine_geom);

  torch::Tensor t1 = torch::zeros({1,ncomp,bx_onegrid_bounds.x+1,bx_onegrid_bounds.y+1}); 

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(SnewRefined,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& SnewRefined_fab = SnewRefined.array(mfi);
    FArrayBox  fab_tmp(bx,ncomp);
    auto       fab  = fab_tmp.array(0);

    // Here we retrieve data from device to a fab array on host    
    amrex::ParallelFor(bx, ncomp, [fab,SnewRefined_fab]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
     fab(i,j,k,n) = SnewRefined_fab(i,j,k,n) ;
    });

    // Now we copy the fab array into the Torch Tensor
    {
    for (int n = 0; n < ncomp; n++ )
    {
      for (int jj = 0; jj < bx_onegrid_bounds.y+1; jj++ )
      {
        for (int ii = 0; ii < bx_onegrid_bounds.x+1; ii++ )
        {
          t1[0][n][ii][jj] = fab(ii,jj,0,n);
        }
      }
    }
//    amrex::Print() << "\n DEBUG SnewRefined " << SnewRefined[mfi] << "\n";

    }
  }

//  amrex::Print() << "Tensor t1 from array :\n" << t1 << '\n';

  /*
    for (MFIter mfi(SnewRefined,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    amrex::Print() << SnewRefined[mfi];
  }
*/

// scale input to match training script TODO: parse the MLPDE-iamr/config/UNet.yaml to directly get scaling factor value.


  t1 = t1*0.1;

  auto bytest1 = torch::jit::pickle_save(t1);
  std::ofstream foutt1("t1.zip", std::ios::out | std::ios::binary);
  foutt1.write(bytest1.data(), bytest1.size());
  foutt1.close();
// Loading the model with TorchScript
  
  torch::jit::script::Module module;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
//    module = torch::jit::load("traced_dummy_module_16.pt");

    module = torch::jit::load("default_sr_large000.pt");

  }
  catch (const c10::Error& e) {
    amrex::Abort( "error loading the model\n");
  }


  std::vector<torch::jit::IValue> inputs{t1};



  at::Tensor output = module.forward(inputs).toTensor();

  auto bytesop = torch::jit::pickle_save(output);
  std::ofstream foutop("output.zip", std::ios::out | std::ios::binary);
  foutop.write(bytesop.data(), bytesop.size());
  foutop.close();

  output += t1;
  auto bytes_op_t1 = torch::jit::pickle_save(output);
  std::ofstream fout_op_t1("output_plus_t1.zip", std::ios::out | std::ios::binary);
  fout_op_t1.write(bytes_op_t1.data(), bytes_op_t1.size());
  fout_op_t1.close();
 // reverse the scaling
  output *= 10;

//  amrex::Print() << "\n DEBUG TENSOR AFTER ML " << output << std::endl;


// Putting back the Torch Tensor to the SnewRefined multifab

  MultiFab CorrectedState(SnewRefined.boxArray(), SnewRefined.DistributionMap(), SnewRefined.nComp(),0);
  CorrectedState.setVal(0.);


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(CorrectedState,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& CorrectedState_fab = CorrectedState.array(mfi);
    FArrayBox  fab_tmp(bx,ncomp);
    auto       fab  = fab_tmp.array(0);

    // Now we copy the Torch Tensor into the temporary fab array
    {
      for (int n = 0; n < ncomp; n++ )
      {
        for (int jj = 0; jj < bx_onegrid_bounds.y+1; jj++ )
        {
          for (int ii = 0; ii < bx_onegrid_bounds.x+1; ii++ )
          {
           fab(ii,jj,0,n) =  output[0][n][ii][jj].item<double>();
          }
        }
      }

//      amrex::Print() << "\n DEBUG fab_tmp_after_ML" << fab_tmp << "\n";

    }

// Here we put back the temporary fab to SnewRefined_fab  
    amrex::ParallelFor(bx, ncomp, [fab,CorrectedState_fab]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
//    amrex::Print() << "\n DEBUG " << i <<  " " << j << " " << n << " " << fab(i,j,k,n) << std::endl; 
      CorrectedState_fab(i,j,k,n) = fab(i,j,k,n);
//       amrex::Print() << "\n DEBUG " << i <<  " " << j << " " << n << " " << fab(i,j,k,n) << " " << CorrectedState_fab(i,j,k,n) << std::endl;
    });
//    amrex::Print() << "\n DEBUG CorrectedState " << CorrectedState[mfi] << "\n";

  }


  Snew_onegrid.setVal(0.);
  
  // TODO fix the time step written into the plotfile. Figure out how to get the current time step?  
  int timestr = cur_time*100000;
  amrex::Print() << timestr << "\n";
  const std::string& hr_pfname = amrex::Concatenate("plotsHR/correctedPltHR",timestr); 
  amrex::WriteSingleLevelPlotfile(hr_pfname, CorrectedState, {"x","y"}, fine_geom, cur_time, 0 );


// We average down CorrectedState to the original Snew_onegrid. Warning, we still have an unused ghost-cell here
  amrex::average_down (CorrectedState, Snew_onegrid, 0,  ncomp, ratio);

//  amrex::Print() << timestr << "\n";
  const std::string& lr_pfname = amrex::Concatenate("plotsLR/correctedPltLR",timestr); 
  amrex::WriteSingleLevelPlotfile(lr_pfname, Snew_onegrid, {"x","y"}, geom, cur_time, 0 );

// Last step is to put back Snew_onegrid in the general State_type data

  Snew.ParallelCopy(Snew_onegrid,0,0,ncomp,0,0);
  FillPatch(*this,Snew,Snew.nGrow(),cur_time,State_Type,0,NUM_STATE);

  // For debugging
  VisMF::Write(Snew,"Snew_coarse_AFTER_ML");
  VisMF::Write(Snew_onegrid,"Snew_onegrid_AFTER_ML");




}

