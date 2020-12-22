/* We evolve forward in time a set of functions called the
 *  "conservative variables", and any time the conserv's
 *  are updated, we must solve for the primitive variables
 *  (rho, pressure, velocities) using a Newton-Raphson
 *  technique, before reconstructing & evaluating the RHSs
 *  of the MHD equations again.
 *
 * This file contains the driver routine for this Newton-
 *  Raphson solver. Truncation errors in conservative
 *  variables can lead to no physical solutions in
 *  primitive variables. We correct for these errors here
 *  through a number of tricks described in the appendices
 *  of http://arxiv.org/pdf/1112.0568.pdf.
 *
 * This is a wrapper for the 2d solver of Noble et al. See
 *  harm_utoprim_2d.c for references and copyright notice
 *  for that solver. This wrapper was primarily written by
 *  Zachariah Etienne & Yuk Tung Liu, in 2011-2013.
 *
 * For optimal compatibility, this wrapper is licensed under
 *  the GPL v2 or any later version.
 *
 * Note that this code assumes a simple gamma law for the
 *  moment, though it would be easy to extend to a piecewise
 *  polytrope. */

// Standard #include's
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>

#if defined(ENABLE_STANDALONE_IGM_C2P_SOLVER)
#include "standalone_conserv_to_prims_main_function.h"
#elif defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#include "standalone_conserv_to_prims_main_function_random_data.h"
#else
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "IllinoisGRMHD_headers.h"
#include "harm_primitives_headers.h"
#include "harm_u2p_util.c"
#include "inlined_functions.C"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"

extern "C" void IllinoisGRMHD_conserv_to_prims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // We use proper C++ here, for file I/O later.
  using namespace std;
#endif

  /**********************************
   * Piecewise Polytropic EOS Patch *
   *   Setting up the EOS struct    *
   **********************************/
  /*
   * The short piece of code below takes care
   * of initializing the EOS parameters.
   * Please refer to the "inlined_functions.C"
   * source file for the documentation on the
   * function.
   */
  eos_struct eos;
  initialize_EOS_struct_from_input(eos);


  // These BSSN-based variables are not evolved, and so are not defined anywhere that the grid has moved.
  // Here we convert ADM variables (from ADMBase) to the BSSN-based variables expected by this routine.
  IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij(cctkGH,cctk_lsh,  gxx,gxy,gxz,gyy,gyz,gzz,alp,
                                                                              gtxx,gtxy,gtxz,gtyy,gtyz,gtzz,
                                                                              gtupxx,gtupxy,gtupxz,gtupyy,gtupyz,gtupzz,
                                                                              phi_bssn,psi_bssn,lapm1);


#if defined(ENABLE_STANDALONE_IGM_C2P_SOLVER) || defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#else
  if(CCTK_EQUALS(Symmetry,"equatorial")) {
    // SET SYMMETRY GHOSTZONES ON ALL CONSERVATIVE VARIABLES!
    int ierr=0;
    ierr+=CartSymGN(cctkGH,"IllinoisGRMHD::grmhd_conservatives");
    // FIXME: UGLY. Filling metric ghostzones is needed for, e.g., Cowling runs.
    ierr+=CartSymGN(cctkGH,"lapse::lapse_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_vars");
    ierr+=CartSymGN(cctkGH,"bssn::BSSN_AH");
    ierr+=CartSymGN(cctkGH,"shift::shift_vars");
    if(ierr!=0) CCTK_VError(VERR_DEF_PARAMS,"IllinoisGRMHD ERROR (grep for it, foo!)  :(");
  }
#endif


  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  //FIXME: Replace this timing benchmark with something more meaningful, like the avg # of Newton-Raphson iterations per gridpoint!
  /*
    struct timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);
  */

  int failures=0,font_fixes=0,vel_limited_ptcount=0;
  int pointcount=0;
  int failures_inhoriz=0;
  int pointcount_inhoriz=0;
  int tau_fixes_applied=0;

  // Moving these up. Will flag major failures
  // if NAN is found
  int ff_avg_P                = 0;
  int ff_not_enough_neighbors = 0;
  int major_averages          = 0;
  int major_failures          = 0;

  // Set a mask to track C2P failures. Initialize to zero.
  // int local_number_of_pts = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
  // int c2p_mask[local_number_of_pts];
  // for(int pt=0;pt<local_number_of_pts;pt++) c2p_mask[pt]=0;

  int pressure_cap_hit=0;

  CCTK_REAL error_int_numer=0,error_int_denom=0;
#if( USE_ENTROPY_EQUATION )
  CCTK_REAL __attribute__((unused)) error_int_numer_1d         = 0;
  CCTK_REAL __attribute__((unused)) error_int_numer_1dee       = 0;
  CCTK_REAL __attribute__((unused)) error_int_numer_1dee2      = 0;
  CCTK_REAL smallest_numer_err = 0;
  int fails_2d=0,fails_1d=0,fails_1dee=0,fails_1dee2=0;
  int  used_2d=0, used_1d=0, used_1dee=0, used_1dee2=0;
#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA
  CCTK_REAL prim_error_int_numer_2d    = 0;
  CCTK_REAL prim_error_int_numer_1d    = 0;
  CCTK_REAL prim_error_int_numer_1dee  = 0;
  CCTK_REAL prim_error_int_numer_1dee2 = 0;
  CCTK_REAL prim_error_int_denom       = 0;
#endif
#endif

  int imin=0,jmin=0,kmin=0;
  int imax=cctk_lsh[0],jmax=cctk_lsh[1],kmax=cctk_lsh[2];

  int rho_star_fix_applied=0;
  long n_iter=0;

#if( USE_ENTROPY_EQUATION )
#if defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,fails_2d,fails_1d,fails_1dee,fails_1dee2,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_denom,error_int_numer_1d,error_int_numer_1dee,error_int_numer_1dee2,pressure_cap_hit,rho_star_fix_applied,n_iter,prim_error_int_numer_2d,prim_error_int_numer_1d,prim_error_int_numer_1dee,prim_error_int_numer_1dee2,major_failures,tau_fixes_applied) schedule(static)
#else
#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,fails_2d,fails_1d,fails_1dee,fails_1dee2,used_2d,used_1d,used_1dee,used_1dee2,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_numer_1d,error_int_numer_1dee,error_int_numer_1dee2,smallest_numer_err,error_int_denom,pressure_cap_hit,rho_star_fix_applied,n_iter,major_failures,tau_fixes_applied) schedule(static)
#endif // defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#else
#pragma omp parallel for reduction(+:failures,vel_limited_ptcount,font_fixes,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_denom,pressure_cap_hit,rho_star_fix_applied,n_iter,major_failures,tau_fixes_applied) schedule(static)
#endif // USE_ENTROPY_EQUATION
  for(int k=kmin;k<kmax;k++)
    for(int j=jmin;j<jmax;j++)
      for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

        int ww;
        CCTK_REAL METRIC[NUMVARS_FOR_METRIC],dummy=0;
        ww=0;
        // FIXME: NECESSARY?
        //psi_bssn[index] = exp(phi[index]);
        METRIC[ww] = phi_bssn[index];ww++;
        METRIC[ww] = dummy;          ww++; // Don't need to set psi.
        METRIC[ww] = gtxx[index];    ww++;
        METRIC[ww] = gtxy[index];    ww++;
        METRIC[ww] = gtxz[index];    ww++;
        METRIC[ww] = gtyy[index];    ww++;
        METRIC[ww] = gtyz[index];    ww++;
        METRIC[ww] = gtzz[index];    ww++;
        METRIC[ww] = lapm1[index];   ww++;
        METRIC[ww] = betax[index];   ww++;
        METRIC[ww] = betay[index];   ww++;
        METRIC[ww] = betaz[index];   ww++;
        METRIC[ww] = gtupxx[index];  ww++;
        METRIC[ww] = gtupyy[index];  ww++;
        METRIC[ww] = gtupzz[index];  ww++;
        METRIC[ww] = gtupxy[index];  ww++;
        METRIC[ww] = gtupxz[index];  ww++;
        METRIC[ww] = gtupyz[index];  ww++;


        CCTK_REAL PRIMS[MAXNUMVARS];
        ww=0;
        PRIMS[ww] = rho_b[index]; ww++;
        PRIMS[ww] = P[index];     ww++;
        PRIMS[ww] = vx[index];    ww++;
        PRIMS[ww] = vy[index];    ww++;
        PRIMS[ww] = vz[index];    ww++;
        PRIMS[ww] = Bx[index];    ww++;
        PRIMS[ww] = By[index];    ww++;
        PRIMS[ww] = Bz[index];    ww++;


#ifndef ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA
        CCTK_REAL CONSERVS[NUM_CONSERVS] = {rho_star[index], mhd_st_x[index],mhd_st_y[index],mhd_st_z[index],tau[index]};
#if( USE_ENTROPY_EQUATION )
        CCTK_REAL S_entropy_local = S_entropy[index];
        CCTK_REAL S_star_local    = S_star[index];
#endif
#endif


        CCTK_REAL METRIC_LAP_PSI4[NUMVARS_METRIC_AUX];
        SET_LAPSE_PSI4(METRIC_LAP_PSI4,METRIC);


        CCTK_REAL METRIC_PHYS[NUMVARS_FOR_METRIC];
        METRIC_PHYS[GXX]   = METRIC[GXX]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GXY]   = METRIC[GXY]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GXZ]   = METRIC[GXZ]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GYY]   = METRIC[GYY]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GYZ]   = METRIC[GYZ]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GZZ]   = METRIC[GZZ]*METRIC_LAP_PSI4[PSI4];
        METRIC_PHYS[GUPXX] = METRIC[GUPXX]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXY] = METRIC[GUPXY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPXZ] = METRIC[GUPXZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYY] = METRIC[GUPYY]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPYZ] = METRIC[GUPYZ]*METRIC_LAP_PSI4[PSIM4];
        METRIC_PHYS[GUPZZ] = METRIC[GUPZZ]*METRIC_LAP_PSI4[PSIM4];


        CCTK_REAL TUPMUNU[10]      ,TDNMUNU[10];
#if( USE_ENTROPY_EQUATION )
        CCTK_REAL TUPMUNU_1D[10]   ,TDNMUNU_1D[10];
        CCTK_REAL TUPMUNU_1DEE[10] ,TDNMUNU_1DEE[10];
        CCTK_REAL TUPMUNU_1DEE2[10],TDNMUNU_1DEE2[10];
#endif

        CCTK_REAL shift_xL = METRIC_PHYS[GXX]*METRIC[SHIFTX] + METRIC_PHYS[GXY]*METRIC[SHIFTY] + METRIC_PHYS[GXZ]*METRIC[SHIFTZ];
        CCTK_REAL shift_yL = METRIC_PHYS[GXY]*METRIC[SHIFTX] + METRIC_PHYS[GYY]*METRIC[SHIFTY] + METRIC_PHYS[GYZ]*METRIC[SHIFTZ];
        CCTK_REAL shift_zL = METRIC_PHYS[GXZ]*METRIC[SHIFTX] + METRIC_PHYS[GYZ]*METRIC[SHIFTY] + METRIC_PHYS[GZZ]*METRIC[SHIFTZ];
        CCTK_REAL beta2L   = shift_xL*METRIC[SHIFTX] + shift_yL*METRIC[SHIFTY] + shift_zL*METRIC[SHIFTZ];


        // Compute 4-metric, both g_{\mu \nu} and g^{\mu \nu}.
        // This is for computing T_{\mu \nu} and T^{\mu \nu}. Also the HARM con2prim lowlevel function requires them.
        CCTK_REAL g4dn[4][4],g4up[4][4];
        g4dn[0][0] = -SQR(METRIC_LAP_PSI4[LAPSE]) + beta2L;
        g4dn[0][1] = g4dn[1][0] = shift_xL;
        g4dn[0][2] = g4dn[2][0] = shift_yL;
        g4dn[0][3] = g4dn[3][0] = shift_zL;
        g4dn[1][1]              = METRIC_PHYS[GXX];
        g4dn[1][2] = g4dn[2][1] = METRIC_PHYS[GXY];
        g4dn[1][3] = g4dn[3][1] = METRIC_PHYS[GXZ];
        g4dn[2][2]              = METRIC_PHYS[GYY];
        g4dn[2][3] = g4dn[3][2] = METRIC_PHYS[GYZ];
        g4dn[3][3]              = METRIC_PHYS[GZZ];

        CCTK_REAL alpha_inv_squared=SQR(METRIC_LAP_PSI4[LAPSEINV]);
        g4up[0][0] = -1.0*alpha_inv_squared;
        g4up[0][1] = g4up[1][0] = METRIC[SHIFTX]*alpha_inv_squared;
        g4up[0][2] = g4up[2][0] = METRIC[SHIFTY]*alpha_inv_squared;
        g4up[0][3] = g4up[3][0] = METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[1][1]              = METRIC_PHYS[GUPXX] - METRIC[SHIFTX]*METRIC[SHIFTX]*alpha_inv_squared;
        g4up[1][2] = g4up[2][1] = METRIC_PHYS[GUPXY] - METRIC[SHIFTX]*METRIC[SHIFTY]*alpha_inv_squared;
        g4up[1][3] = g4up[3][1] = METRIC_PHYS[GUPXZ] - METRIC[SHIFTX]*METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[2][2]              = METRIC_PHYS[GUPYY] - METRIC[SHIFTY]*METRIC[SHIFTY]*alpha_inv_squared;
        g4up[2][3] = g4up[3][2] = METRIC_PHYS[GUPYZ] - METRIC[SHIFTY]*METRIC[SHIFTZ]*alpha_inv_squared;
        g4up[3][3]              = METRIC_PHYS[GUPZZ] - METRIC[SHIFTZ]*METRIC[SHIFTZ]*alpha_inv_squared;

#if( defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA) || defined(ENABLE_STANDALONE_IGM_C2P_SOLVER) )
        struct output_stats stats_dummy;
#if( !USE_ENTROPY_EQUATION )
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(1,
                                                                          PRIMS,stats_dummy,eos,METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS);
#else
        CCTK_REAL S_entropy_local = S_entropy[index];
        CCTK_REAL S_star_local, CONSERVS[NUM_CONSERVS];
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(1,
                                                                          PRIMS,stats_dummy,eos,METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS,
                                                                          &S_entropy_local, &S_star_local );
#endif
#endif


        //FIXME: might slow down the code.
#if( !USE_ENTROPY_EQUATION )
        if(CCTK_isnan(CONSERVS[RHOSTAR]*CONSERVS[STILDEX]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]*CONSERVS[TAUENERGY]*PRIMS[BX_CENTER]*PRIMS[BY_CENTER]*PRIMS[BZ_CENTER])) {
          CCTK_VInfo(CCTK_THORNSTRING,"NAN FOUND: i,j,k = %d %d %d, x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e",
                     i,j,k,x[index],y[index],z[index],index,
                     CONSERVS[STILDEX],CONSERVS[STILDEY],CONSERVS[STILDEZ],CONSERVS[RHOSTAR],CONSERVS[TAUENERGY],
                     PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6]);
        }
#else
        if(CCTK_isnan(CONSERVS[RHOSTAR]*CONSERVS[STILDEX]*CONSERVS[STILDEY]*CONSERVS[STILDEZ]*CONSERVS[TAUENERGY]*PRIMS[BX_CENTER]*PRIMS[BY_CENTER]*PRIMS[BZ_CENTER]*S_star_local)) {
          CCTK_VInfo(CCTK_THORNSTRING,"NAN FOUND: i,j,k = %d %d %d, x,y,z = %e %e %e , index=%d st_i = %e %e %e, rhostar = %e, tau = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e, S_star = %e",
                     i,j,k,x[index],y[index],z[index],index,
                     CONSERVS[STILDEX],CONSERVS[STILDEY],CONSERVS[STILDEZ],CONSERVS[RHOSTAR],CONSERVS[TAUENERGY],
                     PRIMS[BX_CENTER],PRIMS[BY_CENTER],PRIMS[BZ_CENTER],METRIC_PHYS[GXX],METRIC_PHYS[GXY],METRIC_PHYS[GXZ],METRIC_PHYS[GYY],METRIC_PHYS[GYZ],METRIC_PHYS[GZZ],METRIC_LAP_PSI4[PSI6],
                     S_star_local);
            major_failures++;
          /* // Mark this as a failure because we will need to fix this point later. */
          /* failures++; */
          /* c2p_mask[index]=1; */
          /* // Skip the rest of the function */
          /* continue; */
        }
#endif

        // Here we use _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
        rho_star_flux[index] = CONSERVS[RHOSTAR];
        st_x_flux[index]     = CONSERVS[STILDEX];
        st_y_flux[index]     = CONSERVS[STILDEY];
        st_z_flux[index]     = CONSERVS[STILDEZ];
        tau_flux[index]      = CONSERVS[TAUENERGY];
#if( USE_ENTROPY_EQUATION )
        S_star_flux[index]   = S_star_local;
#endif

        CCTK_REAL rho_star_orig = CONSERVS[RHOSTAR];
        CCTK_REAL mhd_st_x_orig = CONSERVS[STILDEX];
        CCTK_REAL mhd_st_y_orig = CONSERVS[STILDEY];
        CCTK_REAL mhd_st_z_orig = CONSERVS[STILDEZ];
        CCTK_REAL tau_orig      = CONSERVS[TAUENERGY];
#if( USE_ENTROPY_EQUATION )
        // CCTK_REAL S_star_orig   = S_star_local;
#endif

#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA
        ww=0;
        CCTK_REAL rhob_orig = PRIMS[ww]; ww++;
        CCTK_REAL P_orig    = PRIMS[ww]; ww++;
        CCTK_REAL vx_orig   = PRIMS[ww]; ww++;
        CCTK_REAL vy_orig   = PRIMS[ww]; ww++;
        CCTK_REAL vz_orig   = PRIMS[ww]; ww++;
        CCTK_REAL Bx_orig   = PRIMS[ww]; ww++;
        CCTK_REAL By_orig   = PRIMS[ww]; ww++;
        CCTK_REAL Bz_orig   = PRIMS[ww]; ww++;
#endif


        int check=0;
        struct output_stats stats;
        stats.n_iter=0;
        stats.vel_limited=0;
        stats.failure_checker=0;

#if( USE_ENTROPY_EQUATION )
        stats.check_2d = stats.check_1d = stats.check_1dee = stats.check_1dee2 = 0;
        CCTK_REAL S_entropy_local_1d    = S_entropy_local;
        CCTK_REAL S_entropy_local_1dee  = S_entropy_local;
        CCTK_REAL S_entropy_local_1dee2 = S_entropy_local;
        CCTK_REAL S_star_local_1d       = S_star_local;
        CCTK_REAL S_star_local_1d_ee    = S_star_local;
        CCTK_REAL S_star_local_1d_ee2   = S_star_local;
        CCTK_REAL PRIMS_1D[MAXNUMVARS]     , PRIMS_1DEE[MAXNUMVARS]     , PRIMS_1DEE2[MAXNUMVARS];
        CCTK_REAL CONSERVS_1D[NUM_CONSERVS], CONSERVS_1DEE[NUM_CONSERVS], CONSERVS_1DEE2[NUM_CONSERVS];
#endif


        stats.font_fixed=0;
#if( USE_ENTROPY_EQUATION )
        stats.fails_2d=stats.fails_1d=stats.fails_1dee=stats.fails_1dee2=0;
#endif
        CCTK_REAL rho_star_atm = METRIC_LAP_PSI4[PSI6]*rho_b_atm;
        CCTK_REAL atmprefactor = 0.0;
        if( CONSERVS[RHOSTAR] > atmprefactor*rho_star_atm ) {
          // Apply the tau floor
            stats.tau_fix_applied=0;
            apply_tau_floor(index,tau_atm,rho_b_atm,Psi6threshold,PRIMS,METRIC,METRIC_PHYS,METRIC_LAP_PSI4,stats,eos,  CONSERVS);
            tau_fixes_applied += stats.tau_fix_applied;

          for(int ii=0;ii<3;ii++) {
#if( !USE_ENTROPY_EQUATION )
            check = harm_primitives_gammalaw_lowlevel(index,i,j,k,x,y,z,METRIC,METRIC_PHYS,METRIC_LAP_PSI4,
                                                      CONSERVS,PRIMS,  g4dn,g4up,   stats,eos);
#else
            check = harm_primitives_gammalaw_lowlevel(index,i,j,k,x,y,z,METRIC,METRIC_PHYS,METRIC_LAP_PSI4,
                                                      CONSERVS, CONSERVS_1D, CONSERVS_1DEE, CONSERVS_1DEE2,
                                                      PRIMS   , PRIMS_1D   , PRIMS_1DEE   , PRIMS_1DEE2   ,
                                                      S_star_local,S_entropy_local,  g4dn,g4up,   stats,eos);
            fails_2d    += stats.fails_2d;
            fails_1d    += stats.fails_1d;
            fails_1dee  += stats.fails_1dee;
            fails_1dee2 += stats.fails_1dee2;
#endif
            if(check==0) ii=4;
            else stats.failure_checker+=100000;

          }
        }
        else {
            // Otherwise, reset to atmosphere
            // rho
            PRIMS[RHOB]     = rho_b_atm;
            int pp_idx      = find_polytropic_K_and_Gamma_index(eos,rho_b_atm);
            CCTK_REAL K     = eos.K_ppoly_tab[pp_idx];
            CCTK_REAL Gamma = eos.Gamma_ppoly_tab[pp_idx];
            // P = P_cold
            PRIMS[PRESSURE] = K * pow(rho_b_atm,Gamma);
            // v^i
            PRIMS[VX] = - METRIC[SHIFTX];
            PRIMS[VY] = - METRIC[SHIFTY];
            PRIMS[VZ] = - METRIC[SHIFTZ];

#if( USE_ENTROPY_EQUATION )
            PRIMS_1D[RHOB] = PRIMS_1DEE[RHOB] = PRIMS_1DEE2[RHOB] = PRIMS[RHOB];

            PRIMS_1D[PRESSURE] = PRIMS_1DEE[PRESSURE] = PRIMS_1DEE2[PRESSURE] = PRIMS[PRESSURE];

            PRIMS_1D[VX] = PRIMS_1DEE[VX] = PRIMS_1DEE2[VX] = PRIMS[VX];
            PRIMS_1D[VY] = PRIMS_1DEE[VY] = PRIMS_1DEE2[VY] = PRIMS[VY];
            PRIMS_1D[VZ] = PRIMS_1DEE[VZ] = PRIMS_1DEE2[VZ] = PRIMS[VZ];
#endif
            rho_star_fix_applied++;
        }


        CCTK_REAL error_int_numerl=0.0;
#if( USE_ENTROPY_EQUATION )
        CCTK_REAL error_int_numer_1dl=0.0,error_int_numer_1deel=0.0,error_int_numer_1dee2l=0.0;
#endif

        // Enforce limits on primitive variables and recompute conservatives.
        static const int already_computed_physical_metric_and_inverse=1;
#if( !USE_ENTROPY_EQUATION )
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,
                                                                          PRIMS,stats,eos,METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS);
#else
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,
                                                                          PRIMS,stats,eos,METRIC,g4dn,g4up, TUPMUNU,TDNMUNU,CONSERVS,
                                                                          &S_entropy_local, &S_star_local );
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,
                                                                          PRIMS_1D,stats,eos,METRIC,g4dn,g4up, TUPMUNU_1D,TDNMUNU_1D,CONSERVS_1D,
                                                                          &S_entropy_local_1d, &S_star_local_1d );
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,
                                                                          PRIMS_1DEE,stats,eos,METRIC,g4dn,g4up, TUPMUNU_1DEE,TDNMUNU_1DEE,CONSERVS_1DEE,
                                                                          &S_entropy_local_1dee, &S_star_local_1d_ee );
        IllinoisGRMHD_enforce_limits_on_primitives_and_recompute_conservs(already_computed_physical_metric_and_inverse,
                                                                          PRIMS_1DEE2,stats,eos,METRIC,g4dn,g4up, TUPMUNU_1DEE2,TDNMUNU_1DEE2,CONSERVS_1DEE2,
                                                                          &S_entropy_local_1dee2, &S_star_local_1d_ee2 );
#endif


        // Now we compute the difference between original & new conservatives, for diagnostic purposes:
        error_int_numerl = fabs(CONSERVS[TAUENERGY] - tau_orig) + fabs(CONSERVS[RHOSTAR] - rho_star_orig) +
          fabs( CONSERVS[STILDEX] - mhd_st_x_orig ) + fabs( CONSERVS[STILDEY] - mhd_st_y_orig ) + fabs( CONSERVS[STILDEZ] - mhd_st_z_orig );
        error_int_denom += tau_orig + rho_star_orig + fabs(mhd_st_x_orig) + fabs(mhd_st_y_orig) + fabs(mhd_st_z_orig);

        /* error_int_numer += error_int_numerl; */

#if( USE_ENTROPY_EQUATION )
        error_int_numer_1dl    = fabs( CONSERVS_1D[TAUENERGY] - tau_orig ) + fabs( CONSERVS_1D[RHOSTAR] - rho_star_orig ) +
          fabs( CONSERVS_1D[STILDEX] - mhd_st_x_orig ) + fabs( CONSERVS_1D[STILDEY] - mhd_st_y_orig ) + fabs( CONSERVS_1D[STILDEZ] - mhd_st_z_orig );
        error_int_numer_1deel  = fabs( CONSERVS_1DEE[TAUENERGY] - tau_orig ) + fabs( CONSERVS_1DEE[RHOSTAR] - rho_star_orig ) +
          fabs( CONSERVS_1DEE[STILDEX] - mhd_st_x_orig ) + fabs( CONSERVS_1DEE[STILDEY] - mhd_st_y_orig ) + fabs( CONSERVS_1DEE[STILDEZ] - mhd_st_z_orig );
        error_int_numer_1dee2l = fabs( CONSERVS_1DEE2[TAUENERGY] - tau_orig ) + fabs( CONSERVS_1DEE2[RHOSTAR] - rho_star_orig ) +
          fabs( CONSERVS_1DEE2[STILDEX] - mhd_st_x_orig ) + fabs( CONSERVS_1DEE2[STILDEY] - mhd_st_y_orig ) + fabs( CONSERVS_1DEE2[STILDEZ] - mhd_st_z_orig );

        /* error_int_numer_1d    += error_int_numer_1dl; */
        /* error_int_numer_1dee  += error_int_numer_1deel; */
        /* error_int_numer_1dee2 += error_int_numer_1dee2l; */
#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA
        prim_error_int_numer_2d    += fabs( PRIMS[RHOB] - rhob_orig )          + fabs( PRIMS[PRESSURE] - P_orig )
          + fabs( PRIMS[VX] - vx_orig )              + fabs( PRIMS[VY] - vy_orig )              + fabs( PRIMS[VZ] - vz_orig )
          + fabs( PRIMS[BX_CENTER] - Bx_orig )       + fabs( PRIMS[BY_CENTER] - By_orig )       + fabs( PRIMS[BZ_CENTER] - Bz_orig );

        prim_error_int_numer_1d    += fabs( PRIMS_1D[RHOB] - rhob_orig )       + fabs( PRIMS_1D[PRESSURE] - P_orig )
          + fabs( PRIMS_1D[VX] - vx_orig )           + fabs( PRIMS_1D[VY] - vy_orig )           + fabs( PRIMS_1D[VZ] - vz_orig )
          + fabs( PRIMS_1D[BX_CENTER] - Bx_orig )    + fabs( PRIMS_1D[BY_CENTER] - By_orig )    + fabs( PRIMS_1D[BZ_CENTER] - Bz_orig );

        prim_error_int_numer_1dee  += fabs( PRIMS_1DEE[RHOB] - rhob_orig )     + fabs( PRIMS_1DEE[PRESSURE] - P_orig )
          + fabs( PRIMS_1DEE[VX] - vx_orig )         + fabs( PRIMS_1DEE[VY] - vy_orig )         + fabs( PRIMS_1DEE[VZ] - vz_orig )
          + fabs( PRIMS_1DEE[BX_CENTER] - Bx_orig )  + fabs( PRIMS_1DEE[BY_CENTER] - By_orig )  + fabs( PRIMS_1DEE[BZ_CENTER] - Bz_orig );

        prim_error_int_numer_1dee2 += fabs( PRIMS_1DEE2[RHOB] - rhob_orig )    + fabs( PRIMS_1DEE2[PRESSURE] - P_orig )
          + fabs( PRIMS_1DEE2[VX] - vx_orig )        + fabs( PRIMS_1DEE2[VY] - vy_orig )        + fabs( PRIMS_1DEE2[VZ] - vz_orig )
          + fabs( PRIMS_1DEE2[BX_CENTER] - Bx_orig ) + fabs( PRIMS_1DEE2[BY_CENTER] - By_orig ) + fabs( PRIMS_1DEE2[BZ_CENTER] - Bz_orig );

        prim_error_int_denom       += rhob_orig + P_orig + fabs(vx_orig) + fabs(vy_orig) + fabs(vz_orig) + fabs(Bx_orig) + fabs(By_orig) + fabs(Bz_orig);
#endif
#endif

#if( USE_ENTROPY_EQUATION )

        if( CCTK_isnan(error_int_numerl)       ) error_int_numerl       = 1e300;
        if( CCTK_isnan(error_int_numer_1dl)    ) error_int_numer_1dl    = 1e300;
        if( CCTK_isnan(error_int_numer_1deel)  ) error_int_numer_1deel  = 1e300;
        if( CCTK_isnan(error_int_numer_1dee2l) ) error_int_numer_1dee2l = 1e300;

        // The denominator is the same for all, so it suffices to compare the numerators
        int which_routine_to_use;
        if( error_int_numerl < error_int_numer_1dl ) {
          if( error_int_numerl < error_int_numer_1deel ) {
            if( error_int_numerl < error_int_numer_1dee2l ) {
              which_routine_to_use = 0;
            }
            else {
              which_routine_to_use = 3;
            }
          }
          else if( error_int_numer_1deel < error_int_numer_1dee2l ) {
            which_routine_to_use = 2;
          }
          else {
            which_routine_to_use = 3;
          }
        }
        else {
          if( error_int_numer_1dl < error_int_numer_1deel ) {
            if( error_int_numer_1dl < error_int_numer_1dee2l ) {
              which_routine_to_use = 1;
            }
            else {
              which_routine_to_use = 3;
            }
          }
          else if( error_int_numer_1deel < error_int_numer_1dee2l ) {
            which_routine_to_use = 2;
          }
          else {
            which_routine_to_use = 3;
          }
        }

        CCTK_REAL *BESTCONSERVS, *BESTPRIMS, *BEST_TDNMUNU, S_entropy_best, S_star_best;
        if( which_routine_to_use == 0 ) {
          BESTCONSERVS       = CONSERVS;
          BESTPRIMS          = PRIMS;
          BEST_TDNMUNU       = TDNMUNU;
          S_entropy_best     = S_entropy_local;
          S_star_best        = S_star_local;
          smallest_numer_err += error_int_numerl;
          used_2d++;
        }
        else if( which_routine_to_use == 1 ) {
          BESTCONSERVS       = CONSERVS_1D;
          BESTPRIMS          = PRIMS_1D;
          BEST_TDNMUNU       = TDNMUNU_1D;
          S_entropy_best     = S_entropy_local_1d;
          S_star_best        = S_star_local_1d;
          smallest_numer_err += error_int_numer_1dl;
          used_1d++;
        }
        else if( which_routine_to_use == 2 ) {
          BESTCONSERVS       = CONSERVS_1DEE;
          BESTPRIMS          = PRIMS_1DEE;
          BEST_TDNMUNU       = TDNMUNU_1DEE;
          S_entropy_best     = S_entropy_local_1dee;
          S_star_best        = S_star_local_1d_ee;
          smallest_numer_err += error_int_numer_1deel;
          used_1dee++;
        }
        else {
          BESTCONSERVS       = CONSERVS_1DEE2;
          BESTPRIMS          = PRIMS_1DEE2;
          BEST_TDNMUNU       = TDNMUNU_1DEE2;
          S_entropy_best     = S_entropy_local_1dee2;
          S_star_best        = S_star_local_1d_ee2;
          smallest_numer_err += error_int_numer_1dee2l;
          used_1dee2++;
        }
#else
        CCTK_REAL *BESTCONSERVS = CONSERVS;
        CCTK_REAL *BESTPRIMS    = PRIMS;
        CCTK_REAL *BEST_TDNMUNU = TDNMUNU;
#endif

        // Update conservatives, in case they needed to be modified.
        rho_star[index] = BESTCONSERVS[RHOSTAR];
        mhd_st_x[index] = BESTCONSERVS[STILDEX];
        mhd_st_y[index] = BESTCONSERVS[STILDEY];
        mhd_st_z[index] = BESTCONSERVS[STILDEZ];
        tau[index]      = BESTCONSERVS[TAUENERGY];
#if( USE_ENTROPY_EQUATION )
        S_star[index]   = S_star_best;
#endif

        rho_b[index] = BESTPRIMS[RHOB];
        P[index]     = BESTPRIMS[PRESSURE];
        vx[index]    = BESTPRIMS[VX];
        vy[index]    = BESTPRIMS[VY];
        vz[index]    = BESTPRIMS[VZ];
#if( USE_ENTROPY_EQUATION )
        S_entropy[index] = S_entropy_best;
#endif


        if(update_Tmunu) {
          ww=0;
          eTtt[index] = BEST_TDNMUNU[ww]; ww++;
          eTtx[index] = BEST_TDNMUNU[ww]; ww++;
          eTty[index] = BEST_TDNMUNU[ww]; ww++;
          eTtz[index] = BEST_TDNMUNU[ww]; ww++;
          eTxx[index] = BEST_TDNMUNU[ww]; ww++;
          eTxy[index] = BEST_TDNMUNU[ww]; ww++;
          eTxz[index] = BEST_TDNMUNU[ww]; ww++;
          eTyy[index] = BEST_TDNMUNU[ww]; ww++;
          eTyz[index] = BEST_TDNMUNU[ww]; ww++;
          eTzz[index] = BEST_TDNMUNU[ww];
        }

        vel_limited_ptcount+=stats.vel_limited;
        pointcount++;
        /***************************************************************************************************************************/
        failure_checker[index] = stats.failure_checker;
        n_iter += stats.n_iter;
      }

  if(CCTK_Equals(verbose, "essential") || CCTK_Equals(verbose, "essential+iteration output")) {
#if( !USE_ENTROPY_EQUATION )
    CCTK_VInfo(CCTK_THORNSTRING,"C2P: Lev: %d NumPts= %d | Fixes: Font= %d VL= %d rho*= %d | Failures: %d InHoriz= %d / %d | Error: %.3e, ErrDenom: %.3e | %.2f iters/gridpt",
               (int)GetRefinementLevel(cctkGH),
               pointcount,font_fixes,vel_limited_ptcount,rho_star_fix_applied,
               failures,
               failures_inhoriz,pointcount_inhoriz,
               error_int_numer/error_int_denom,error_int_denom,
               (double)n_iter/( (double)(cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]) ));
#else
    CCTK_VInfo(CCTK_THORNSTRING,"C2P: Lev: %d Npts=%d | 2d=%d/%d 1d=%d/%d ee=%d/%d e2=%d/%d | Fixes: VL=%d rho*=%d FF=%d+%d=%d TF=%d PC=%d | Fails: %d+%d=%d(I:%d/%d) | Error: %.3e, ErrDenom: %.3e",
               (int)GetRefinementLevel(cctkGH),
               pointcount,
               fails_2d,used_2d,
               fails_1d,used_1d,
               fails_1dee,used_1dee,
               fails_1dee2,used_1dee2,
               vel_limited_ptcount,
               rho_star_fix_applied,
               ff_not_enough_neighbors,ff_avg_P,font_fixes,tau_fixes_applied,pressure_cap_hit,
               major_failures,major_averages,failures,
	       failures_inhoriz,pointcount_inhoriz,
               smallest_numer_err/error_int_denom,error_int_denom);
#endif
  }

  //if(pressure_cap_hit!=0) {
    //CCTK_VInfo(CCTK_THORNSTRING,"PRESSURE CAP HIT %d TIMES!  Outputting debug file!",pressure_cap_hit);
  //}

  // Very useful con2prim debugger. If the primitives (con2prim) solver fails, this will output all data needed to
  //     debug where and why the solver failed. Strongly suggested for experimenting with new fixes.
  if( (conserv_to_prims_debug==1 && error_int_numer/error_int_denom > 0.05) || failures || major_failures ) {

    ofstream myfile;
    char filename[100];
    srand(time(NULL));
    if( failures || major_failures ) {
      sprintf(filename,"major_failure-%d_%d-%e.dat",failures,major_failures,smallest_numer_err/error_int_denom);
    }
    else {
      sprintf(filename,"primitives_debug-%e.dat",smallest_numer_err/error_int_denom);
    }
    //Alternative, for debugging purposes as well:
    //srand(time(NULL));
    //sprintf(filename,"primitives_debug-%d.dat",rand());
    myfile.open (filename, ios::out | ios::binary);
    //myfile.open ("data.bin", ios::out | ios::binary);
    myfile.write((char*)cctk_lsh, 3*sizeof(int));

    myfile.write((char*)&GAMMA_SPEED_LIMIT, 1*sizeof(CCTK_REAL));

    myfile.write((char*)&rho_b_max, 1*sizeof(CCTK_REAL));
    myfile.write((char*)&rho_b_atm, 1*sizeof(CCTK_REAL));
    myfile.write((char*)&tau_atm, 1*sizeof(CCTK_REAL));

    myfile.write((char*)&Psi6threshold, 1*sizeof(CCTK_REAL));

    myfile.write((char*)&update_Tmunu, 1*sizeof(int));

    myfile.write((char*)&neos,                   1*sizeof(int));
    myfile.write((char*)&Gamma_th,               1*sizeof(CCTK_REAL));
    myfile.write((char*)&K_ppoly_tab0,            1*sizeof(CCTK_REAL));
    myfile.write((char*)Gamma_ppoly_tab_in,   neos*sizeof(CCTK_REAL));
    myfile.write((char*)rho_ppoly_tab_in, (neos-1)*sizeof(CCTK_REAL));

    int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
    myfile.write((char*)x,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)y,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)z,   (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char *)failure_checker, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTtt, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTtx, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTty, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTtz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTxx, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTxy, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTxz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTyy, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTyz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)eTzz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)alp, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)gxx, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)gxy, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)gxz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)gyy, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)gyz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)gzz, fullsize*sizeof(CCTK_REAL));
    myfile.write((char *)psi_bssn, fullsize*sizeof(CCTK_REAL));

    myfile.write((char*)phi_bssn, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxx, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtxz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtyy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtyz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtzz, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)gtupxx, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupxz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupyy, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupyz, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)gtupzz, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)betax, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betay, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)betaz, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)lapm1, (fullsize)*sizeof(CCTK_REAL));

    // HERE WE USE _flux variables as temp storage for original values of conservative variables.. This is used for debugging purposes only.
    myfile.write((char*)tau_flux,      (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_x_flux, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_y_flux, (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)st_z_flux, (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)rho_star_flux, (fullsize)*sizeof(CCTK_REAL));
#if( USE_ENTROPY_EQUATION )
    myfile.write((char*)S_star_flux, (fullsize)*sizeof(CCTK_REAL));
#endif

    myfile.write((char*)Bx,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)By,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)Bz,   (fullsize)*sizeof(CCTK_REAL));

    myfile.write((char*)vx,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vy,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)vz,   (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)P,    (fullsize)*sizeof(CCTK_REAL));
    myfile.write((char*)rho_b,(fullsize)*sizeof(CCTK_REAL));
#if( USE_ENTROPY_EQUATION )
    myfile.write((char*)S_entropy, (fullsize)*sizeof(CCTK_REAL));
#endif

    int checker=1063; myfile.write((char*)&checker,sizeof(int));

    myfile.close();
    CCTK_VInfo(CCTK_THORNSTRING,"Finished writing %s",filename);
  }
  if( major_failures ) {
#if( !defined(ENABLE_STANDALONE_IGM_C2P_SOLVER) && !defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA) )
    CCTK_VError(VERR_DEF_PARAMS,"Major failure! None of the fixes worked or NAN found! Terminating the run...\n");
#else
    fprintf(stderr,"Major failure! None of the fixes worked or NAN found! Terminating the run...\n");
    std::exit(1);
#endif
  }

#ifdef ENABLE_STANDALONE_IGM_C2P_SOLVER
  return 0; // int main() requires an integer be returned
#endif

}

#include "harm_primitives_lowlevel.C"
