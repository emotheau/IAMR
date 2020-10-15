#include <iamr_godunov.H>


using namespace amrex;


void
Godunov::ComputeAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                       MultiFab const& state, const int state_comp,
                       AMREX_D_DECL( MultiFab const& umac,
                                     MultiFab const& vmac,
                                     MultiFab const& wmac),
                       AMREX_D_DECL( MultiFab& xedge,
                                     MultiFab& yedge,
                                     MultiFab& zedge),
                       const int  edge_comp,
                       const bool known_edgestate,
                       AMREX_D_DECL( MultiFab& xfluxes,
                                     MultiFab& yfluxes,
                                     MultiFab& zfluxes),
                       int fluxes_comp,
                       MultiFab const& fq,
                       const int fq_comp,
                       MultiFab const& divu,
                       Vector<BCRec> const& bcs,
                       Geometry const& geom,
                       Gpu::DeviceVector<int>& iconserv,
                       const Real dt,
                       const bool use_ppm,
                       const bool use_forces_in_trans,
                       const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofs()");

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();
        Box const& bxg1 = amrex::grow(bx,1);

        FArrayBox tmpfab(amrex::grow(bx,1),  (4*AMREX_SPACEDIM + 2)*ncomp);

        //
        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp);,
                      const auto& fy = yfluxes.array(mfi,fluxes_comp);,
                      const auto& fz = zfluxes.array(mfi,fluxes_comp););

        AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp);,
                      const auto& yed = yedge.array(mfi,edge_comp);,
                      const auto& zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                      const auto& v = vmac.const_array(mfi);,
                      const auto& w = wmac.const_array(mfi););

        if (not known_edgestate)
        {
            ComputeEdgeState( bx, ncomp,
                              state.array(mfi,state_comp),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              divu.array(mfi),
                              fq.array(mfi,fq_comp),
                              geom, dt, &bcs[0],
                              iconserv.data(),
                              use_ppm,
                              use_forces_in_trans,
                              is_velocity );
        }

        ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
                       AMREX_D_DECL( u, v, w ),
                       AMREX_D_DECL( xed, yed, zed ),
                       geom, ncomp );


        ComputeDivergence( bx,
                           aofs.array(mfi,aofs_comp),
                           AMREX_D_DECL( fx, fy, fz ),
                           AMREX_D_DECL( xed, yed, zed ),
                           AMREX_D_DECL( u, v, w ),
                           ncomp, geom, iconserv.data() );

        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}

#ifdef PLM_USE_EFIELD
void
Godunov::ComputeAofsEF ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                         MultiFab const& state, const int state_comp,
                         AMREX_D_DECL( MultiFab const& umac,
                                       MultiFab const& vmac,
                                       MultiFab const& wmac),
                         AMREX_D_DECL( MultiFab const& udrift,
                                       MultiFab const& vdrift,
                                       MultiFab const& wdrift),
                         AMREX_D_DECL( MultiFab& xedge,
                                       MultiFab& yedge,
                                       MultiFab& zedge),
                         const int  edge_comp,
                         const bool known_edgestate,
                         AMREX_D_DECL( MultiFab& xfluxes,
                                       MultiFab& yfluxes,
                                       MultiFab& zfluxes),
                         int fluxes_comp,
                         MultiFab const& fq,
                         const int fq_comp,
                         MultiFab const& divu,
                         Vector<BCRec> const& bcs,
                         Geometry const& geom,
                         Gpu::DeviceVector<int>& iconserv,
                         const Real dt,
                         const bool use_ppm,
                         const bool use_forces_in_trans,
                         const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofsEF()");

    FArrayBox Ueff[AMREX_SPACEDIM];

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();
        Box const& bxg1 = amrex::grow(bx,1);

        FArrayBox tmpfab(amrex::grow(bx,1),  (4*AMREX_SPACEDIM + 2)*ncomp);

        for ( int n = 0; n < ncomp; n++ ) {
           //
           // Get handlers to Array4
           //
           AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp+n);,
                         const auto& fy = yfluxes.array(mfi,fluxes_comp+n);,
                         const auto& fz = zfluxes.array(mfi,fluxes_comp+n););

           AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp+n);,
                         const auto& yed = yedge.array(mfi,edge_comp+n);,
                         const auto& zed = zedge.array(mfi,edge_comp+n););

           AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                         const auto& v = vmac.const_array(mfi);,
                         const auto& w = wmac.const_array(mfi););

           // Get the effective velocities
           AMREX_D_TERM( Ueff[0].resize(umac[mfi].box(),1);,
                         Ueff[1].resize(vmac[mfi].box(),1);,
                         Ueff[2].resize(wmac[mfi].box(),1););
           AMREX_D_TERM( const auto& udr = udrift.const_array(mfi,n);,
                         const auto& vdr = vdrift.const_array(mfi,n);,
                         const auto& wdr = wdrift.const_array(mfi,n););
           AMREX_D_TERM( const auto& uef = Ueff[0].array();,
                         const auto& vef = Ueff[1].array();,
                         const auto& wef = Ueff[2].array(););
           amrex::ParallelFor(umac[mfi].box(), [uef,udr,u]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
              uef(i,j,k) = u(i,j,k) + udr(i,j,k);
           });
           amrex::ParallelFor(vmac[mfi].box(), [vef,vdr,v]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
              vef(i,j,k) = v(i,j,k) + vdr(i,j,k);
           });
#if (AMREX_SPACEDIM == 3)
           amrex::ParallelFor(wmac[mfi].box(), [wef,wdr,w]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
              wef(i,j,k) = w(i,j,k) + wdr(i,j,k);
           });
#endif

           // Get a local single component version of iconserv
           Gpu::DeviceVector<int> iconserv_comp;
           iconserv_comp.resize(1, 0);
           iconserv_comp[0] = iconserv[n];

           if (not known_edgestate)
           {
               ComputeEdgeState( bx, 1,
                                 state.array(mfi,state_comp+n),
                                 AMREX_D_DECL( xed, yed, zed ),
                                 AMREX_D_DECL( uef, vef, wef ),
                                 divu.array(mfi),
                                 fq.array(mfi,fq_comp+n),
                                 geom, dt, &bcs[n],
                                 iconserv_comp.data(),
                                 use_ppm,
                                 use_forces_in_trans,
                                 is_velocity );
           }

           ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
                          AMREX_D_DECL( uef, vef, wef ),
                          AMREX_D_DECL( xed, yed, zed ),
                          geom, 1 );

           // Scalar used with EF are conserved so the velocity is not used
           ComputeDivergence( bx,
                              aofs.array(mfi,aofs_comp+n),
                              AMREX_D_DECL( fx, fy, fz ),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              1, geom, iconserv_comp.data() );
        } // End loop on ncomp

        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}
#endif


void
Godunov::ComputeSyncAofs ( MultiFab& aofs, const int aofs_comp, const int ncomp,
                           MultiFab const& state, const int state_comp,
                           AMREX_D_DECL( MultiFab const& umac,
                                         MultiFab const& vmac,
                                         MultiFab const& wmac),
                           AMREX_D_DECL( MultiFab const& ucorr,
                                         MultiFab const& vcorr,
                                         MultiFab const& wcorr),
                           AMREX_D_DECL( MultiFab& xedge,
                                         MultiFab& yedge,
                                         MultiFab& zedge),
                           const int  edge_comp,
                           const bool known_edgestate,
                           AMREX_D_DECL( MultiFab& xfluxes,
                                         MultiFab& yfluxes,
                                         MultiFab& zfluxes),
                           int fluxes_comp,
                           MultiFab const& fq,
                           const int fq_comp,
                           MultiFab const& divu,
                           Vector<BCRec> const& bcs,
                           Geometry const& geom,
                           Gpu::DeviceVector<int>& iconserv,
                           const Real dt,
                           const bool use_ppm,
                           const bool use_forces_in_trans,
                           const bool is_velocity  )
{
    BL_PROFILE("Godunov::ComputeAofs()");

    //FIXME - check on adding tiling here
    for (MFIter mfi(aofs); mfi.isValid(); ++mfi)
    {

        const Box& bx   = mfi.tilebox();
        Box const& bxg1 = amrex::grow(bx,1);

        FArrayBox tmpfab(amrex::grow(bx,1),  (4*AMREX_SPACEDIM + 2)*ncomp);

        //
        // Get handlers to Array4
        //
        AMREX_D_TERM( const auto& fx = xfluxes.array(mfi,fluxes_comp);,
                      const auto& fy = yfluxes.array(mfi,fluxes_comp);,
                      const auto& fz = zfluxes.array(mfi,fluxes_comp););

        AMREX_D_TERM( const auto& xed = xedge.array(mfi,edge_comp);,
                      const auto& yed = yedge.array(mfi,edge_comp);,
                      const auto& zed = zedge.array(mfi,edge_comp););

        AMREX_D_TERM( const auto& uc = ucorr.const_array(mfi);,
                      const auto& vc = vcorr.const_array(mfi);,
                      const auto& wc = wcorr.const_array(mfi););

        if (not known_edgestate)
        {

            AMREX_D_TERM( const auto& u = umac.const_array(mfi);,
                          const auto& v = vmac.const_array(mfi);,
                          const auto& w = wmac.const_array(mfi););

            ComputeEdgeState( bx, ncomp,
                              state.array(mfi,state_comp),
                              AMREX_D_DECL( xed, yed, zed ),
                              AMREX_D_DECL( u, v, w ),
                              divu.array(mfi),
                              fq.array(mfi,fq_comp),
                              geom, dt, &bcs[0],
                              iconserv.data(),
                              use_ppm,
                              use_forces_in_trans,
                              is_velocity );
        }

        ComputeFluxes( bx, AMREX_D_DECL( fx, fy, fz ),
                       AMREX_D_DECL( uc, vc, wc ),
                       AMREX_D_DECL( xed, yed, zed ),
                       geom, ncomp );


        ComputeSyncDivergence( bx,
                               aofs.array(mfi,aofs_comp),
                               AMREX_D_DECL( fx, fy, fz ),
                               ncomp, geom );

        Gpu::streamSynchronize();  // otherwise we might be using too much memory
    }

}


void
Godunov::ComputeFluxes ( Box const& bx,
                         AMREX_D_DECL( Array4<Real> const& fx,
                                       Array4<Real> const& fy,
                                       Array4<Real> const& fz),
                         AMREX_D_DECL( Array4<Real const> const& umac,
                                       Array4<Real const> const& vmac,
                                       Array4<Real const> const& wmac),
                         AMREX_D_DECL( Array4<Real const> const& xed,
                                       Array4<Real const> const& yed,
                                       Array4<Real const> const& zed),
                         Geometry const& geom, const int ncomp )
{

    const auto dx = geom.CellSizeArray();

    Real area[AMREX_SPACEDIM];
#if ( AMREX_SPACEDIM == 3 )
    area[0] = dx[1]*dx[2];
    area[1] = dx[0]*dx[2];
    area[2] = dx[0]*dx[1];
#else
    area[0] = dx[0];
    area[1] = dx[1];
#endif

    //
    //  X flux
    //
    const Box& xbx = amrex::surroundingNodes(bx,0);

    amrex::ParallelFor(xbx, ncomp, [fx, umac, xed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fx(i,j,k,n) = xed(i,j,k,n) * umac(i,j,k) * area[0];
    });

    //
    //  y flux
    //
    const Box& ybx = amrex::surroundingNodes(bx,1);

    amrex::ParallelFor(ybx, ncomp, [fy, vmac, yed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fy(i,j,k,n) = yed(i,j,k,n) * vmac(i,j,k) * area[1];
    });

#if (AMREX_SPACEDIM==3)
    //
    //  z flux
    //
    const Box& zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(zbx, ncomp, [fz, wmac, zed, area]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        fz(i,j,k,n) = zed(i,j,k,n) * wmac(i,j,k) * area[2];
    });

#endif

}



void
Godunov::ComputeDivergence ( Box const& bx,
                             Array4<Real> const& div,
                             AMREX_D_DECL( Array4<Real const> const& fx,
                                           Array4<Real const> const& fy,
                                           Array4<Real const> const& fz),
                             AMREX_D_DECL( Array4<Real const> const& xed,
                                           Array4<Real const> const& yed,
                                           Array4<Real const> const& zed),
                             AMREX_D_DECL( Array4<Real const> const& umac,
                                           Array4<Real const> const& vmac,
                                           Array4<Real const> const& wmac),
                             const int ncomp, Geometry const& geom,
                             int const* iconserv )
{

    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM==3)
    Real qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
    Real qvol = dxinv[0] * dxinv[1];
#endif

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (iconserv[n])
        {
            div(i,j,k,n) =  qvol *
                (
                         fx(i+1,j,k,n) -  fx(i,j,k,n)
                       + fy(i,j+1,k,n) -  fy(i,j,k,n)
#if (AMREX_SPACEDIM==3)
                       + fz(i,j,k+1,n) -  fz(i,j,k,n)
#endif
                );
        }
        else
        {
            div(i,j,k,n) = 0.5*dxinv[0]*( umac(i+1,j,k  ) +  umac(i,j,k  ))
                *                       (  xed(i+1,j,k,n) -   xed(i,j,k,n))
                +          0.5*dxinv[1]*( vmac(i,j+1,k  ) +  vmac(i,j,k  ))
                *                       (  yed(i,j+1,k,n) -   yed(i,j,k,n))
#if (AMREX_SPACEDIM==3)
                +          0.5*dxinv[2]*( wmac(i,j,k+1  ) +  wmac(i,j,k  ))
                *                       (  zed(i,j,k+1,n) -   zed(i,j,k,n))
#endif
                ;
       }

    });
}


void
Godunov::ComputeSyncDivergence ( Box const& bx,
                                 Array4<Real> const& div,
                                 AMREX_D_DECL( Array4<Real const> const& fx,
                                               Array4<Real const> const& fy,
                                               Array4<Real const> const& fz),
                                 const int ncomp, Geometry const& geom )
{

    const auto dxinv = geom.InvCellSizeArray();

#if (AMREX_SPACEDIM==3)
    Real qvol = dxinv[0] * dxinv[1] * dxinv[2];
#else
    Real qvol = dxinv[0] * dxinv[1];
#endif

    amrex::ParallelFor(bx, ncomp,[=]
    AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        div(i,j,k,n) +=  qvol * (
              fx(i+1,j,k,n) -  fx(i,j,k,n)
            + fy(i,j+1,k,n) -  fy(i,j,k,n)
#if (AMREX_SPACEDIM==3)
            + fz(i,j,k+1,n) -  fz(i,j,k,n)
#endif
            );
    });
}
