
#include <NavierStokesBase.H>
#include <AMReX_VisMF.H>

#include <torch/torch.h>
#include <torch/script.h> // One-stop header.


using namespace amrex;



void
NavierStokesBase::test_libtorch()

{

std::cout << "\n WE ARE IN TEST_LIBTORCH ROUTINE \n";


// Getting new state data
MultiFab& Snew   = get_new_data(State_Type);
// get boxarray associated with the new state multifab
BoxArray ba = Snew.boxArray();

// create a copy of the boxarray and refine it. Use it to create a new Multifab to hold refined data
int ratio = 2;

BoxArray baRefined = ba;

baRefined.refine(ratio);

DistributionMapping dmapRefined(baRefined);

MultiFab SnewRefined(baRefined, dmapRefined, Snew.nComp(), Snew.nGrow());

int ncomps = Snew.nComp();

// get geometry object from global scope and create a refined copy
Box bx_geom = geom.Domain();

bx_geom.refine(ratio);

Geometry fine_geom(bx_geom, geom.ProbDomain(), geom.Coord(), geom.isPeriodic());

// instantiate interpolator
Interpolater*  interpolater = &cell_cons_interp;

// iterate through boxes in multifab and interpolate to refined multifab
IntVect new_ratio(AMREX_D_DECL(ratio,ratio,ratio));
const IntVect& rr = new_ratio;

for (MFIter mfi(SnewRefined); mfi.isValid(); ++mfi)
{
  FArrayBox& ffab = (SnewRefined)[mfi];
  const FArrayBox& cfab = (Snew)[mfi];
  const Box&  bx   = mfi.tilebox();
  Vector<BCRec> bx_bcrec(ncomps);

  interpolater->interp(cfab,0,ffab,0,ncomps,bx,rr,
                       geom,fine_geom,bx_bcrec,0,0,RunOn::Host);
}

// minimalBox() computes a single box to enclose all the boxes
// enclosedCells() converts it to a cell-centered Box

Box bx_onegrid = SnewRefined.boxArray().minimalBox().enclosedCells();

// number of cells in the coarse domain
Print() << "npts in coarse domain = " << bx_onegrid.numPts() << std::endl;

// BoxArray, DistributionMapping, and MultiFab with one grid
BoxArray ba_onegrid(bx_onegrid);
DistributionMapping dmap_onegrid(ba_onegrid);
MultiFab Snew_onegrid(ba_onegrid,dmap_onegrid,SnewRefined.nComp(),0);

// copy data into MultiFab with one grid
Snew_onegrid.ParallelCopy(SnewRefined,0,0,SnewRefined.nComp(),0,0);


VisMF::Write(Snew,"Snew");
VisMF::Write(Snew_onegrid,"Snew_onegrid");

amrex::Print() << "\n DEBUG amrex::ubound(bx) = " << amrex::ubound(bx_onegrid) << std::endl;
amrex::Print() << "\n DEBUG amrex::lbound(bx) = " << amrex::lbound(bx_onegrid) << std::endl;

const auto bx_onegrid_bounds = amrex::ubound(bx_onegrid);

amrex::Print() << "\n DEBUG box_bounds.x = " << bx_onegrid_bounds.x << std::endl;

amrex::Print() << "\n DEBUG ncomps = " << ncomps << std::endl;
// We want to launch Torch on just 1 mpi process
//if (ParallelDescriptor::IOProcessor()){


torch::Tensor t1 = torch::zeros({1,2,bx_onegrid_bounds.x+1,bx_onegrid_bounds.y+1}); 


#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(Snew_onegrid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& Snew_vel = Snew_onegrid.array(mfi);

    amrex::ParallelFor(bx, ncomps, [t1, Snew_vel]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      if (n<2) { 
      t1[0][n][j][i] = Snew_vel(i,j,k,n);  // Inverting i and j index because of the PyTorch formalism
      }
    });

  }

torch::jit::script::Module module;
  try {
    // Deserialize the ScriptModule from a file using torch::jit::load().
    module = torch::jit::load("traced_unet_model.pt");
  }
  catch (const c10::Error& e) {
    amrex::Abort( "error loading the model\n");
  }


    std::vector<torch::jit::IValue> inputs{t1};

at::Tensor output = module.forward(inputs).toTensor();
output += t1;
//} // End of ParallelDescriptor

amrex::Print() << "DEBUG WE ARE AFTER TORCH \n ";

MultiFab CorrectedState(Snew_onegrid.boxArray(), Snew_onegrid.DistributionMap(), Snew_onegrid.nComp(), Snew_onegrid.nGrow());

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(CorrectedState,TilingIfNotGPU()); mfi.isValid(); ++mfi)
  {
    const Box& bx = mfi.tilebox();
    auto const& correctedVelocity = CorrectedState.array(mfi);

    amrex::ParallelFor(bx, ncomps, [output, correctedVelocity]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      if (n<2) { 
      correctedVelocity(i,j,k,n) = output[0][n][j][i].item<double>() ;  // Inverting i and j index because of the PyTorch formalism
      }
    });

  }

amrex::average_down (CorrectedState, Snew,
                             0,  ncomps, ratio);



//ParallelDescriptor::Barrier();

}

