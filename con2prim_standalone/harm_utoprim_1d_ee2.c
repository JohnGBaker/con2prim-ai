#ifndef __HARM_UTOPRIM_1D_EE2__C__
#define __HARM_UTOPRIM_1D_EE2__C__

/**********************************************************************************
 * This is a modified version of the original HARM code which is compatible with
 * IllinoisGRMHD. The modifications below are documented in pedagogical Jupyter
 * notebooks, which are available at: www.nrpyplus.net
 *
 * Author: Leo R. Werneck (wernecklr@gmail.com)
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

utoprim_1d_ee2.c:
---------------

  -- uses eq. (27) of Noble  et al. or the "momentum equation" and ignores
        the energy equation (29) in order to use the additional EOS, which
        is

             P = Sc rho^(GAMMA-1) / gamma

    Uses a method similiar to  1D_W method:
       -- solves for one independent variable (rho) via a 1D
          Newton-Raphson method
       -- by substituting
          W = Dc ( Dc + GAMMA Sc rho^(GAMMA-1) / (GAMMA-1) ) / rho
           into Qtsq equation, one can get one equation for
           one unknown (rho)

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
static int Utoprim_new_body2(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                             CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet,  CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL gamma_times_S );

static void func_rho2(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                      CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                      eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W_in);

/**********************************************************************/
/******************************************************************

  Utoprim_1d_ee2():

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
       gcov[NDIM][NDIM] = covariant form of the metric ;
       gcon[NDIM][NDIM] = contravariant form of the metric ;
       gdet             = sqrt( - determinant of the metric) ;
       prim[NPR] = primitive variables (guess on input, calculated values on
                                        output if there are no problems);

   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.  ;

******************************************************************/

int Utoprim_1d_ee2( eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL gcon[NDIM][NDIM],
                    CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL S_star )
{

  if( U[0] <= 0. ) {
    fprintf(stderr,"Negative U[0] found!!  We encourage you to figure out how this weird thing happened!! \n");
    return(-100);
  }

  /* First update the primitive B-fields */
  CCTK_REAL inv_gdet = 1.0 / gdet;
  for(int i = BCON1; i <= BCON3; i++) prim[i] = U[i] * inv_gdet ;

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
  U_tmp[BCON1 ]           = geomfactor * U[BCON1] ;
  U_tmp[BCON2 ]           = geomfactor * U[BCON2] ;
  U_tmp[BCON3 ]           = geomfactor * U[BCON3] ;


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
  int ret = Utoprim_new_body2(eos, U_tmp, gcov, gcon, gdet, prim_tmp, n_iter, gamma_times_S);

  /* Transform new primitive variables back if there was no problem : */
  if( ret == 0 ) {
    prim[RHO   ] = prim_tmp[RHO   ];
    prim[UU    ] = prim_tmp[UU    ];
    prim[UTCON1] = prim_tmp[UTCON1];
    prim[UTCON2] = prim_tmp[UTCON2];
    prim[UTCON3] = prim_tmp[UTCON3];
  }

  return( ret ) ;

}


/**********************************************************************/
/**********************************************************************************

  Utoprim_new_body():

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
                   (occurrence of "nan" or "+/-inf" ;

         j = 0 -> success
             1 -> failure: some sort of failure in Newton-Raphson;
             2 -> failure: utsq<0 w/ initial p[] guess;
             3 -> failure: W<0 or W>W_TOO_BIG
             4 -> failure: v^2 > 1
             5 -> failure: rho,uu <= 0 ;

**********************************************************************************/

static int Utoprim_new_body2(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                             CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet,  CCTK_REAL prim[NPR], long &n_iter, CCTK_REAL gamma_times_S )
{

  // Assume ok initially:
  int retval = 0;

  // Calculate various scalars (Q.B, Q^2, etc)  from the conserved variables:
  // B^{\mu}
  CCTK_REAL Bcon[NDIM];
  Bcon[0] = 0. ;
  Bcon[1] = U[BCON1];
  Bcon[2] = U[BCON2];
  Bcon[3] = U[BCON3];

  // Leo says: declare harm auxiliary variables struct
  // This struct is defined in the harm_utoprim_1d_ee.c file
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
  /* calculate W from last timestep and use for guess */
  CCTK_REAL utsq = 0.0;
  for(int i=1;i<4;i++)
    for(int j=1;j<4;j++) utsq += gcov[i][j]*prim[UTCON1+i-1]*prim[UTCON1+j-1];

  if( (utsq < 0.) && (fabs(utsq) < 1.0e-13) ) {
    utsq = fabs(utsq);
  }
  if(utsq < 0.0 || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval) ;
  }

  CCTK_REAL gammasq = 1.0 + utsq;    // Lorentz factor squared
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

  CCTK_REAL x_1d[1] = {rho0};
  retval = general_newton_raphson( x_1d, 1, eos, harm_aux, W_last, func_rho2 );
  rho0 = x_1d[0];

  /* Problem with solver, so return denoting error before doing anything further */
  if( (retval != 0) || (rho0 == FAIL_VAL) ) {
    retval = retval*100+1;
    return(retval);
  }
  else{
    if(rho0 > W_TOO_BIG) {
      retval = 3;
      return(retval) ;
    }
  }

  // Calculate v^2 :
  poly_idx   = find_polytropic_K_and_Gamma_index(eos,rho0);
  Gamma_poly = eos.Gamma_ppoly_tab[poly_idx];
  Gm1        = Gamma_poly - 1.0;
  rho_Gm1    = pow(rho0,Gm1);

  CCTK_REAL rel_err = (harm_aux.D != 0.0) ? fabs((harm_aux.D-rho0)/harm_aux.D) : ( (rho0 != 0.0) ? fabs((harm_aux.D-rho0)/rho0) : 0.0 );
  utsq = ( rel_err > 1e-15 ) ? (harm_aux.D-rho0)*(harm_aux.D+rho0)/(rho0*rho0) : 0.0;

  gammasq        = 1.+utsq;
  harm_aux.gamma = sqrt(gammasq);
  //*gamma_out = gamma;

  if( utsq < 0. ) {
    retval = 4;
    return(retval);
  }

  // Recover the primitive variables from the scalars and conserved variables:
  p = harm_aux.gamma_times_S * rho_Gm1 / harm_aux.gamma;
  u = p / Gm1;
  w = rho0 + u + p;
  CCTK_REAL W = w * gammasq;


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
  return(retval) ;

}

/**********************************************************************/
/*********************************************************************************
   func_rho():

        -- residual/jacobian routine to calculate rho Qtsq equation with
            the definition of W
        W  =  ( 1 + GAMMA * K_atm * rho^(GAMMA-1)/(GAMMA-1) ) D^2 / rho
              substituted in.

     Arguments:
          x   = current value of independent var's (on input & output);
         dx   = Newton-Raphson step (on output);
        resid = residuals based on x (on output);
         jac  = Jacobian matrix based on x (on output);
         f    =  resid.resid/2  (on output)
        df    = -2*f;  (on output)
         n    = dimension of x[];
 *********************************************************************************/
static void func_rho2(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                      CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                      eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W_in)
{



  // Set rho
  CCTK_REAL rho = x[0];

  // Auxiliary variable declarations
  // Gamma
  CCTK_REAL Gamma = eos.Gamma_ppoly_tab[find_polytropic_K_and_Gamma_index(eos,rho)];
  // Gamma-1
  CCTK_REAL Gm1 = Gamma - 1.0;
  // rho^2
  CCTK_REAL rhosq = rho*rho;
  // Gamma/(Gamma-1) * gamma * S * rho^(Gamma-1)
  CCTK_REAL t100 = (Gamma/Gm1) * harm_aux.gamma_times_S * pow(rho,Gm1);
  // D^(-2)
  CCTK_REAL t200 = 1.0/(harm_aux.D*harm_aux.D);
  // W
  CCTK_REAL W = harm_aux.D * ( harm_aux.D + t100 ) / rho;
  // dWdrho
  CCTK_REAL dWdrho = harm_aux.D * ( -harm_aux.D + t100*(Gamma-2.0) ) / rhosq;
  // W^2
  CCTK_REAL t1 = W*W;
  // B^2 + W
  CCTK_REAL t2 = harm_aux.Bsq+W;
  // (B^2 + W)^2
  CCTK_REAL t3 = t2*t2;
  // v^2
  CCTK_REAL rel_err = (harm_aux.D != 0.0) ? fabs((harm_aux.D-rho)/harm_aux.D) : ( ( rho != 0.0 ) ? fabs((harm_aux.D-rho)/rho) : 0.0 );
  CCTK_REAL vsq = ( rel_err > 1e-15 ) ? ((harm_aux.D-rho)*(harm_aux.D+rho)*t200) : 0.0;
  // d(v^2)/drho
  CCTK_REAL dvsqdrho = -2*rho*t200;
  // B^4
  CCTK_REAL t12 = harm_aux.Bsq*harm_aux.Bsq;
  // v^2 * dwrho
  CCTK_REAL t17 = dWdrho*vsq;

  // Compute the residual and the needed Jacobian component
  resid[0]  = t1*(harm_aux.Qtsq-vsq*t3)+harm_aux.QdotBsq*(t2+W);
  jac[0][0] = 2*harm_aux.QdotBsq*dWdrho
    +((harm_aux.Qtsq-vsq*t12)*2*dWdrho+(-6*t17*harm_aux.Bsq-dvsqdrho*t12
                      +(-2*dvsqdrho*harm_aux.Bsq-4*t17-dvsqdrho*W)*W)*W)*W;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = - 0.5*(*df);

  return;

}

/******************************************************************************
             END   OF   UTOPRIM_1D_EE2.C
 ******************************************************************************/
#endif // __HARM_UTOPRIM_1D_EE2__C__
