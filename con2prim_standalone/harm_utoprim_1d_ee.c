#ifndef __HARM_UTOPRIM_1D_EE__C__
#define __HARM_UTOPRIM_1D_EE__C__

/**********************************************************************************
 * This is a modified version of the original HARM code which is compatible with
 * IllinoisGRMHD. The modifications below are documented in pedagogical Jupyter
 * notebooks, which are available at: www.nrpyplus.net
 *
 * Author: Leo Werneck (leow155@gmail.com)
 * Date  : July/August 2020
 **********************************************************************************/

/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble,
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to
    solve the relativistic magnetohydrodynamic equations of motion on a
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model.

    You are morally obligated to cite the following two papers in his/her
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003,
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006,
        Astrophysical Journal, 641, 626.


    Further, we strongly encourage you to obtain the latest version of
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1d_ee.c:
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS, which
        is

             P = Sc rho^(GAMMA-1) / gamma

    Uses the 1D_W method:
       -- solves for one independent variable (W) via a 1D
          Newton-Raphson method
       -- solves for rho using Newton-Raphson using the definition of W :
          W = Dc ( Dc + GAMMA Sc rho^(GAMMA-1) / (GAMMA-1) ) / rho

       -- can be used (in principle) with a general equation of state.

  -- Currently returns with an error state (>0) if a negative rest-mass
      density or internal energy density is calculated.  You may want
      to change this aspect of the code so that it still calculates the
      velocity and so that you can floor the densities.  If you want to
      change this aspect of the code please comment out the "return(retval)"
      statement after "retval = 5;" statement in Utoprim_new_body();

******************************************************************************/

// If NEWT_DIM is already defined, make sure it has
// the correct value (1, in this case); otherwise,
// redefine it.
#ifdef  NEWT_DIM
#if(    NEWT_DIM != (1) )
#undef  NEWT_DIM
#define NEWT_DIM (1)
#endif
#endif

// Declarations:
static int Utoprim_new_body_ee(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                               CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL gamma_times_S );
static void func_W(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                   CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                   eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL rho_in);

static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                     CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                     eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W_in);

/**********************************************************************/
/******************************************************************

  Utoprim_1d_ee():

  -- Driver for new prim. var. solver.  The driver just translates
     between the two sets of definitions for U and P.  The user may
     wish to alter the translation as they see fit.

     It assumes that on input/output:


              /  rho u^t           \
         U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))
              |  T^t_\mu           |
              \   B^i              /


             /    rho        \
	 P = |    uu         |
             | \tilde{u}^i   |
             \   B^i         /


         S = sqrt(-det(g_{a b})) u^t P / rho^(GAMMA-1)

     ala HARM.

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
         S       = entropy density  = sqrt(-det(g_{a b})) u^t P / rho^(GAMMA-1)
       gcov[NDIM][NDIM] = covariant form of the metric;
       gcon[NDIM][NDIM] = contravariant form of the metric;
       gdet             = sqrt( - determinant of the metric);
       prim[NPR] = primitive variables (guess on input, calculated values on
                                        output if there are no problems);

   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.;

******************************************************************/

int Utoprim_1d_ee(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL gcon[NDIM][NDIM],
                  CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL S_star )
{

  if( U[0] <= 0. ) {
    fprintf(stderr,"Negative U[0] found!!  We encourage you to figure out how this weird thing happened!! \n");
    return(-100);
  }

  /* First update the primitive B-fields */
  CCTK_REAL inv_gdet = 1.0 / gdet;
  for(int i = BCON1; i <= BCON3; i++) prim[i] = U[i] * inv_gdet;

  /* Set the geometry variables: */
  CCTK_REAL alpha      = 1.0/sqrt(-gcon[0][0]);
  CCTK_REAL geomfactor = alpha * inv_gdet;

  /* Transform the CONSERVED variables into the new system
   *
   * The input CONSERVED entropy variable is
   *
   * S_star := sqrt(-g) * S * u0
   *
   * Using u0 = gamma/alpha, with gamma the Lorentz factor,
   * and sqrt(-g), we get
   *
   * alpha / sqrt(-g) * S_star = alpha / sqrt(-g) * ( sqrt(-g) * S * u0 )
   *                           = ( gamma/alpha ) * alpha * S
   *                           = gamma * S,
   * and thus
   * .---------------------------------.
   * | gamma * S = geomfactor * S_star |
   * .---------------------------------.
   */
  CCTK_REAL gamma_times_S = geomfactor * S_star;
  CCTK_REAL U_tmp[NPR];
  U_tmp[RHO]              = geomfactor * U[RHO];
  U_tmp[UU]               = geomfactor * (U[UU] - U[RHO]);
  U_tmp[UTCON1]           = geomfactor * U[UTCON1];
  U_tmp[UTCON2]           = geomfactor * U[UTCON2];
  U_tmp[UTCON3]           = geomfactor * U[UTCON3];
  U_tmp[BCON1 ]           = geomfactor * U[BCON1];
  U_tmp[BCON2 ]           = geomfactor * U[BCON2];
  U_tmp[BCON3 ]           = geomfactor * U[BCON3];

  /* Transform the PRIMITIVE variables into the new system */
  CCTK_REAL prim_tmp[NPR];
  prim_tmp[RHO   ] = prim[RHO   ];
  prim_tmp[UU    ] = prim[UU    ];
  prim_tmp[UTCON1] = prim[UTCON1];
  prim_tmp[UTCON2] = prim[UTCON2];
  prim_tmp[UTCON3] = prim[UTCON3];
  prim_tmp[BCON1]  = U_tmp[BCON1];
  prim_tmp[BCON2]  = U_tmp[BCON2];
  prim_tmp[BCON3]  = U_tmp[BCON3];

  /* Perform the C2P inversion */
  int ret = Utoprim_new_body_ee(eos, U_tmp, gcov, gcon, gdet, prim_tmp, n_iter, gamma_times_S);

  /* Transform new primitive variables back if there was no problem : */
  if( ret == 0 ) {
    prim[RHO   ] = prim_tmp[RHO   ];
    prim[UU    ] = prim_tmp[UU    ];
    prim[UTCON1] = prim_tmp[UTCON1];
    prim[UTCON2] = prim_tmp[UTCON2];
    prim[UTCON3] = prim_tmp[UTCON3];
  }

  return( ret );

}

/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body_ee():

     -- Attempt an inversion from U to prim using the initial guess prim.

     -- This is the main routine that calculates auxiliary quantities for the
        Newton-Raphson routine.

  -- assumes that
             /  rho gamma        \
         U = |  alpha T^t_\mu    |
             \  alpha B^i        /



               /    rho        \
        prim = |    uu         |
               | \tilde{u}^i   |
               \  alpha B^i   /


return:  (i*100 + j)  where
         i = 0 ->  Newton-Raphson solver either was not called (yet or not used)
                   or returned successfully;
             1 ->  Newton-Raphson solver did not converge to a solution with the
                   given tolerances;
             2 ->  Newton-Raphson procedure encountered a numerical divergence
                   (occurrence of "nan" or "+/-inf";

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: utsq<0 w/ initial p[] guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0;

**********************************************************************************/

static int Utoprim_new_body_ee(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                               CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL gamma_times_S )
{

  // Assume ok initially:
  int retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  // B^{\mu}
  CCTK_REAL Bcon[NDIM];
  Bcon[0] = 0.;
  Bcon[1] = U[BCON1];
  Bcon[2] = U[BCON2];
  Bcon[3] = U[BCON3];

  // Leo says: declare harm auxiliary variables struct
  harm_aux_vars_struct harm_aux;
  // gamma_times_S
  harm_aux.gamma_times_S = gamma_times_S;
  // B^2 = g_{\mu\nu}B^{\mu}B^{\nu} (note that B^{0} = 0)
  harm_aux.Bsq = gcov[1][1]*Bcon[1]*Bcon[1]
               + gcov[2][2]*Bcon[2]*Bcon[2]
               + gcov[3][3]*Bcon[3]*Bcon[3]
            + 2*(gcov[1][2]*Bcon[1]*Bcon[2]
               + gcov[1][3]*Bcon[1]*Bcon[3]
               + gcov[2][3]*Bcon[2]*Bcon[3]);

  // Q_{\mu}
  CCTK_REAL Qcov[NDIM];
  Qcov[0] = U[QCOV0];
  Qcov[1] = U[QCOV1];
  Qcov[2] = U[QCOV2];
  Qcov[3] = U[QCOV3];

  // Q^{\mu}
  CCTK_REAL Qcon[NDIM];
  raise_g(Qcov,gcon,Qcon);

  // Q.B = Q_{\mu}B^{\mu} (again, note that B^{0} = 0)
  harm_aux.QdotB   = Qcov[1]*Bcon[1] + Qcov[2]*Bcon[2] + Qcov[3]*Bcon[3];

  // (Q.B)^2 = (Q.B) * (Q.B)
  harm_aux.QdotBsq = harm_aux.QdotB*harm_aux.QdotB;

  // g^{00} = alpha^{-2} => alpha = (g^{00})^{-1/2} = 1 / sqrt(-g^{00})
  CCTK_REAL alpha = 1.0/sqrt(-gcon[0][0]);

  // Q.n = n_{\mu}Q^{\mu} = -alpha Q^{0}, since n_{\mu} = (-alpha,0,0,0)
  harm_aux.Qdotn = -alpha * Qcon[0];

  // Q^2 = Q.Q = Q_{\mu}Q^{\mu}
  harm_aux.Qsq = 0.0;
  for(int mu=0;mu<4;mu++) harm_aux.Qsq += Qcov[mu]*Qcon[mu];

  // \tilde{Q}^{2} = Q^{2} + (Q.n)^{2}
  harm_aux.Qtsq = harm_aux.Qsq + harm_aux.Qdotn*harm_aux.Qdotn;

  harm_aux.D = U[RHO];

  /* calculate W from last timestep and use for guess */
  CCTK_REAL utsq = 0.0;
  for(int i=1;i<4;i++)
    for(int j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1];

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  }
  if(utsq < 0.0 || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval);
  }

  CCTK_REAL gammasq = 1.0 + utsq;   // Lorentz factor squared
  harm_aux.gamma    = sqrt(gammasq); // Lorentz factor, not to be confused with Gamma
  CCTK_REAL rho0    = harm_aux.D / harm_aux.gamma;

  /* Set basic cold polytrope-based hybrid EOS quantities */
  int poly_idx         = find_polytropic_K_and_Gamma_index(eos,rho0); // Which piece of the PPEOS to use
  CCTK_REAL Gamma_poly = eos.Gamma_ppoly_tab[poly_idx];               // Set Gamma
  CCTK_REAL Gm1        = Gamma_poly - 1.0;                            // HARM auxiliary variable
  CCTK_REAL Gm1_inv    = 1.0/Gm1;                                     // HARM auxiliary variable
  CCTK_REAL rho_Gm1    = pow(rho0,Gm1);                               // HARM auxiliary variable

  /* The definition of the entropy density, S, is
   *
   * S = P / rho^(Gamma - 1)
   *
   * Thus we have
   * .-------------------------.
   * | P = rho^(Gamma - 1) * S |
   * .-------------------------.
   */
  CCTK_REAL p = harm_aux.gamma_times_S * rho_Gm1 / harm_aux.gamma;
  CCTK_REAL u = p * Gm1_inv;
  CCTK_REAL w = rho0 + u + p;

  CCTK_REAL W_last = w*gammasq;

  // Make sure that W is large enough so that v^2 < 1 :
  int i_increase = 0;
  while( (( W_last*W_last*W_last * ( W_last + 2.*harm_aux.Bsq )
	    - harm_aux.QdotBsq*(2.*W_last + harm_aux.Bsq) ) <= W_last*W_last*(harm_aux.Qtsq-harm_aux.Bsq*harm_aux.Bsq))
	 && (i_increase < 10) ) {
    W_last *= 10.;
    i_increase++;
  }

  // Calculate W:
  CCTK_REAL x_1d[1] = {W_last};

  retval = general_newton_raphson( x_1d, 1, eos, harm_aux, rho0, func_W );

  CCTK_REAL W = x_1d[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (W == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{
    if(W <= 0. || W > W_TOO_BIG) {
      retval = 3;
      return(retval);
    }
  }

  CCTK_REAL rho_g    =  rho0;
  CCTK_REAL x_rho[1] = {rho0};

  int ntries = 0;
  while (  (retval = general_newton_raphson( x_rho, 1, eos, harm_aux, W, func_rho )) &&  ( ntries++ < 10 )  ) {
    rho_g *= 10.;
    x_rho[0] = rho_g;
  }

  if( (retval != 0) ) {
    retval = 10;
    return(retval);
  }

  // Calculate v^2 :
  rho0       = x_rho[0];
  poly_idx   = find_polytropic_K_and_Gamma_index(eos,rho0);
  Gamma_poly = eos.Gamma_ppoly_tab[poly_idx];
  Gm1        = Gamma_poly - 1.0;
  rho_Gm1    = pow(rho0,Gm1);

  CCTK_REAL rel_err = (harm_aux.D != 0.0) ? fabs((harm_aux.D-rho0)/harm_aux.D) : ( (rho0 != 0.0) ? fabs((harm_aux.D-rho0)/rho0) : 0.0 );
  utsq = ( rel_err > 1e-15 ) ? (harm_aux.D-rho0)*(harm_aux.D+rho0)/(rho0*rho0) : 0.0;

  gammasq        = 1.+utsq;
  harm_aux.gamma = sqrt(gammasq);

  if( utsq < 0. ) {
    retval = 4;
    return(retval);
  }

  // Recover the primitive variables from the scalars and conserved variables:

  w = W / gammasq;
  p = harm_aux.gamma_times_S * rho_Gm1 / harm_aux.gamma;
  u = p * Gm1_inv;

  // User may want to handle this case differently, e.g. do NOT return upon
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  int treat_floor_as_failure = 0;
  if( treat_floor_as_failure && ((rho0 <= 0.) || (u <= 0.)) ) {
    retval = 5;
    return(retval);
  }

  prim[RHO] = rho0;
  prim[UU]  = u;

  CCTK_REAL g_o_WBsq = harm_aux.gamma/(W+harm_aux.Bsq);
  CCTK_REAL QdB_o_W  = harm_aux.QdotB / W;

  // n_{mu} = (-alpha,0,0,0)
  // n^{mu} = g^{munu}n_{nu} = g^{mu 0}n_{0} = -alpha * g^{mu 0)
  CCTK_REAL ncon[4];
  for(int i=1;i<4;i++)
    ncon[i] = - alpha * gcon[i][0]; // No need to set the 0th component

  // Update \tilde{u}^{i}
  prim[UTCON1] = g_o_WBsq * ( Qcon[1] + ncon[1] * harm_aux.Qdotn + QdB_o_W*Bcon[1] );
  prim[UTCON2] = g_o_WBsq * ( Qcon[2] + ncon[2] * harm_aux.Qdotn + QdB_o_W*Bcon[2] );
  prim[UTCON3] = g_o_WBsq * ( Qcon[3] + ncon[3] * harm_aux.Qdotn + QdB_o_W*Bcon[3] );

  /* Done! */
  return(retval);

}

/**********************************************************************/
/*********************************************************************************
   func_W()

        -- calculates the residuals, and Newton step for general_newton_raphson();
        -- for this method, x=W here;

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void func_W(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                   CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                   eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL rho_in)
{

  // Set W from input
  CCTK_REAL W = x[0];

  // Set input rho for NR
  CCTK_REAL rho      =  rho_in;
  CCTK_REAL x_rho[1] = {rho_in};
  // Find rho from W
  int retval, ntries = 0;
  while (  (retval = general_newton_raphson( x_rho, 1, eos, harm_aux, W, func_rho)) &&  ( ntries++ < 10 )  ) {
    rho     *= 10.;
    x_rho[0] = rho;
  }
  // Set rho to the NR output
  rho = x_rho[0];

  // Auxiliary variable declarations
  // Gamma
  CCTK_REAL Gamma = eos.Gamma_ppoly_tab[find_polytropic_K_and_Gamma_index(eos,rho)];
  // Gamma-1
  CCTK_REAL Gm1 = Gamma - 1.0;
  // B^2 + B^2
  CCTK_REAL two_Bsq = harm_aux.Bsq + harm_aux.Bsq;
  // D^2
  CCTK_REAL t2 = harm_aux.D*harm_aux.D;
  // (Q.B)^2 * D^2
  CCTK_REAL t4 = harm_aux.QdotBsq*t2;
  // rho^2
  CCTK_REAL t6 = rho*rho; // Leo says: added, since will be used twice
  // B^4
  CCTK_REAL t7 = harm_aux.Bsq*harm_aux.Bsq;
  // rho^2 - D^2
  CCTK_REAL t15 = t6 - t2; // Leo says: changed from -(D-rho)*(D+rho)
  // D^(-2)
  CCTK_REAL t24 = 1/t2;
  // W + 2B^2
  CCTK_REAL t200 = W + two_Bsq;
  // (Q.B)^2 * B^2 * D^2
  CCTK_REAL t300 = harm_aux.Bsq*t4; // Leo says: changed from QdotBsq*Bsq*t2
  // \tilde{Q}^2 * D^2
  CCTK_REAL t400 = harm_aux.Qtsq*t2;
  // D * Gamma * gamma * S
  CCTK_REAL s200 = harm_aux.D * Gamma * harm_aux.gamma_times_S;
  // rho^(Gamma - 1)
  CCTK_REAL rho_Gm1 = pow(rho,Gm1);
  // -rho^2 / ( W * rho - rho^(Gamma-1) * D * Gamma * gamma * S )
  CCTK_REAL drho_dW = -t6/( -rho_Gm1*s200 + W*rho); // Leo says: changed rho*rho -> t6
  // rho^(Gamma-1) * drho_dW (see variable above)
  CCTK_REAL t1000 = rho*drho_dW;

  // Compute the residual and the needed Jacobian component
  resid[0]  = (t300+(t4+t4+(t400+t15*(t7+(t200)*W))*W)*W)*t24;
  jac[0][0] = 2*(t4+(t400+t15*t7+(3.0*t15*harm_aux.Bsq+t7*t1000+(t15+t15+t1000*(t200))*W)*W)*W)*t24;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = - 0.5*(*df);

  return;

}

/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho from W via the polytrope:

        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
      IGM update:
        eos   = struct containing polytropic eos quantities
gamma_times_S = HARM's Sc variable, set in the Utoprim_1d_ee() function
        W_in  = HARM's W_for_gnr2 variable, set in the Utoprim_new_body() function
          D   = gamma * rho, set in the Utoprim_new_body() functoin
 *********************************************************************************/
// for the isentropic version: eq. (27)
static void func_rho(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                     CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                     eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W_in)
{

  // Set rho and W
  CCTK_REAL rho = x[0];
  CCTK_REAL W   = W_in;

  // Auxiliary variable declarations
  // Gamma
  CCTK_REAL Gamma = eos.Gamma_ppoly_tab[find_polytropic_K_and_Gamma_index(eos,rho)];
  // Gamma-1
  CCTK_REAL Gm1 = Gamma - 1.0;
  // rho^(Gamma-1)
  CCTK_REAL t40 = pow(rho,Gm1);
  // rho^(Gamma-2)
  CCTK_REAL t14 = t40/rho;
  // Gamma/(Gamma-1) * gamma * S
  CCTK_REAL s100 = Gamma/Gm1 * harm_aux.gamma_times_S;
  // D * Gamma * gamma * S
  CCTK_REAL s200 = harm_aux.D * Gamma * harm_aux.gamma_times_S;

  // Compute the residual and the needed Jacobian component
  resid[0] = (rho*W+(-t40*s100-harm_aux.D)*harm_aux.D);
  jac[0][0] = -t14*s200 + W;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = - 0.5*(*df);

  return;

}

/******************************************************************************
             END   OF   UTOPRIM_1D_EE.C
 ******************************************************************************/
#endif
