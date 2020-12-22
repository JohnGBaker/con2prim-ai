#ifndef __HARM_UTOPRIM_1D__C__
#define __HARM_UTOPRIM_1D__C__

/**********************************************************************************
 * This is a modified version of the original HARM code which is compatible with
 * IllinoisGRMHD. The modifications below are documented in pedagogical Jupyter
 * notebooks, which are available at: www.nrpyplus.net
 *
 * Author: Leo Werneck (leow155@gmail.com)
 * Date  : July/August 2020
 **********************************************************************************/
/*
  -------------------------------------------------------------------------------
  Copyright 2005 Scott C. Noble, Charles F. Gammie,
  Jonathan C. McKinney, and Luca Del Zanna


  This file is part of PVS-GRMHD.

  PVS-GRMHD is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  PVS-GRMHD is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with PVS-GRMHD; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------------------
*/

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************

utoprim_1d.c:
---------------

    Uses the 1D_W method:
       -- solves for one independent variable (W) via a 1D
          Newton-Raphson method
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

typedef struct _harm_auxiliary_vars_ {
  CCTK_REAL Bsq,QdotBsq,Qsq,Qtsq,Qdotn,QdotB,D,gamma,gamma_times_S;
} harm_aux_vars_struct;


// Declarations:
static int Utoprim_new_body_1d( eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                                CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter );

static void func_1d_orig(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                         CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                         eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL dummy);

static int general_newton_raphson( CCTK_REAL x[], int n, eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL indep_var_in,
                                   void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                                                  CCTK_REAL [][NEWT_DIM], CCTK_REAL *,
                                                  CCTK_REAL *, int,
                                                  eos_struct, harm_aux_vars_struct, CCTK_REAL) );

static CCTK_REAL vsq_calc_1d(harm_aux_vars_struct harm_aux, CCTK_REAL W);
static CCTK_REAL pressure_W_vsq_1d(eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W, CCTK_REAL vsq);
static CCTK_REAL dpdW_calc_vsq_1d(CCTK_REAL W, CCTK_REAL vsq);
static CCTK_REAL dpdvsq_calc_1d(eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W, CCTK_REAL vsq);
static CCTK_REAL dvsq_dW_1d(harm_aux_vars_struct harm_aux, CCTK_REAL W);

/**********************************************************************/
/******************************************************************

  Utoprim_1d():

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


     ala HARM.

   Arguments:
       U[NPR]    = conserved variables (current values on input/output);
       gcov[NDIM][NDIM] = covariant form of the metric;
       gcon[NDIM][NDIM] = contravariant form of the metric;
       gdet             = sqrt( - determinant of the metric);
       prim[NPR] = primitive variables (guess on input, calculated values on
                                        output if there are no problems);

   -- NOTE: for those using this routine for special relativistic MHD and are
            unfamiliar with metrics, merely set
              gcov = gcon = diag(-1,1,1,1)  and gdet = 1.;

******************************************************************/

int Utoprim_1d(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM], CCTK_REAL gcon[NDIM][NDIM],
               CCTK_REAL gdet, CCTK_REAL prim[NPR], long &n_iter)
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

  /* Transform the CONSERVED variables into the new system */
  CCTK_REAL U_tmp[NPR];
  U_tmp[RHO]    = geomfactor * U[RHO];
  U_tmp[UU]     = geomfactor * (U[UU] - U[RHO]);
  U_tmp[UTCON1] = geomfactor * U[UTCON1];
  U_tmp[UTCON2] = geomfactor * U[UTCON2];
  U_tmp[UTCON3] = geomfactor * U[UTCON3];
  U_tmp[BCON1 ] = geomfactor * U[BCON1];
  U_tmp[BCON2 ] = geomfactor * U[BCON2];
  U_tmp[BCON3 ] = geomfactor * U[BCON3];


  /* Transform the PRIMITIVE variables into the new system */
  CCTK_REAL prim_tmp[NPR];
  prim_tmp[RHO   ] = prim[RHO   ];
  prim_tmp[UU    ] = prim[UU    ];
  prim_tmp[UTCON1] = prim[UTCON1];
  prim_tmp[UTCON2] = prim[UTCON2];
  prim_tmp[UTCON3] = prim[UTCON3];
  prim_tmp[BCON1]  =  U_tmp[BCON1 ];
  prim_tmp[BCON2]  =  U_tmp[BCON2 ];
  prim_tmp[BCON3]  =  U_tmp[BCON3 ];

  /* Perform the C2P inversion */
  int ret = Utoprim_new_body_1d(eos, U_tmp, gcov, gcon, gdet, prim_tmp, n_iter);

  /* Use new primitive variables if there was no problem : */
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

  Utoprim_new_body_1d():

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

static int Utoprim_new_body_1d(eos_struct eos, CCTK_REAL U[NPR], CCTK_REAL gcov[NDIM][NDIM],
                               CCTK_REAL gcon[NDIM][NDIM], CCTK_REAL gdet,  CCTK_REAL prim[NPR], long &n_iter)
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
  if(utsq < 0. || utsq > UTSQ_TOO_BIG) {
    retval = 2;
    return(retval);
  }

  CCTK_REAL gammasq = 1.0 + utsq;   // Lorentz factor squared
  harm_aux.gamma    = sqrt(gammasq); // Lorentz factor, not to be confused with Gamma

  // Always calculate rho from D and gamma so that using D in EOS remains consistent
  //   i.e. you don't get positive values for dP/d(vsq) .
  CCTK_REAL rho0    = harm_aux.D / harm_aux.gamma;

  // EOS quantities: internal energy, pressure, enthalpy
  CCTK_REAL u = prim[UU];
  CCTK_REAL p = pressure_rho0_u(eos,rho0,u);
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

  // We need a dummy variable to keep the function call in this
  // file consistent with the ones in the harm_utoprim_1d_ee.c
  // and harm_utoprim_1d_ee2.c files.
  CCTK_REAL dummy = 0.0;
  retval = general_newton_raphson( x_1d, 1, eos, harm_aux, dummy, func_1d_orig );

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

  // Calculate v^2 :
  CCTK_REAL vsq = vsq_calc_1d(harm_aux,W);
  if( vsq >= 1. ) {
    retval = 4;
    return(retval);
  }

  // Recover the primitive variables from the scalars and conserved variables:
  CCTK_REAL gtmp = sqrt(1. - vsq);
  harm_aux.gamma = 1./gtmp;
  rho0           = harm_aux.D * gtmp;

  w = W * (1. - vsq);
  p = pressure_rho0_w(eos,rho0,w);
  u = w - (rho0 + p);

  // User may want to handle this case differently, e.g. do NOT return upon
  // a negative rho/u, calculate v^i so that rho/u can be floored by other routine:
  int treat_floor_as_failure = 0;
  if( treat_floor_as_failure && ((rho0 <= 0.) || (u <= 0.)) ) {
    retval = 5;
    return(retval);
  }

  prim[RHO] = rho0;
  prim[UU] = u;

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

  /* set field components */
  for(int i = BCON1; i <= BCON3; i++) prim[i] = U[i];

  /* Done! */
  return(retval);

}

/**********************************************************************/
/************************************************************

  general_newton_raphson():

    -- performs Newton-Rapshon method on an arbitrary system.

    -- inspired in part by Num. Rec.'s routine newt();

    Arguements:

       -- x[]   = set of independent variables to solve for;
       -- n     = number of independent variables and residuals;
       -- funcd = name of function that calculates residuals, etc.;

*****************************************************************/
static int general_newton_raphson( CCTK_REAL x[], int n, eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL indep_var_in,
                                   void (*funcd) (CCTK_REAL [], CCTK_REAL [], CCTK_REAL [],
                                   CCTK_REAL [][NEWT_DIM], CCTK_REAL *,
                                   CCTK_REAL *, int,
                                   eos_struct, harm_aux_vars_struct, CCTK_REAL) )
{

  // Initialize various parameters and variables:
  CCTK_REAL errx  = 1.0;
  CCTK_REAL f     = 1.0;
  CCTK_REAL df    = 1.0;
  int i_extra     = 0;
  int doing_extra = 0;

  int n_iter = 0;

  /* Start the Newton-Raphson iterations : */
  int keep_iterating = 1;
  CCTK_REAL dx[NEWT_DIM], resid[NEWT_DIM], jac[NEWT_DIM][NEWT_DIM];
  while( keep_iterating ) {

    /* returns with new dx, f, df */
    (*funcd) (x, dx, resid, jac, &f, &df, n, eos, harm_aux, indep_var_in );

    /* don't use line search : */
    x[0] += dx[0];

    /****************************************/
    /* Calculate the convergence criterion */
    /****************************************/

    /* For the new criterion, always look at error in "W" : */
    // METHOD specific:
    errx  = (x[0]==0.0) ? fabs(dx[0]) : fabs(dx[0]/x[0]);

    /********************************************/
    /* Make sure that the new x[] is physical : */
    /********************************************/
    x[0] = fabs(x[0]);

    /*****************************************************************************/
    /* If we've reached the tolerance level, then just do a few extra iterations */
    /*  before stopping                                                          */
    /*****************************************************************************/
    if( (fabs(errx) <= NEWT_TOL) && (doing_extra == 0) && (EXTRA_NEWT_ITER > 0) ) {
      doing_extra = 1;
    }

    if( doing_extra == 1 ) i_extra++;

    if( ((fabs(errx) <= NEWT_TOL)&&(doing_extra == 0)) ||
        (i_extra > EXTRA_NEWT_ITER) || (n_iter >= (MAX_NEWT_ITER-1)) ) {
      keep_iterating = 0;
    }

    n_iter++;

  }   // END of while(keep_iterating)


  /*  Check for bad untrapped divergences : */
  if( (!std::isfinite(f)) || (!std::isfinite(df)) || (!std::isfinite(x[0]))  ) {
    return(2);
  }

  if( fabs(errx) <= NEWT_TOL ){
    return(0);
  }
  else if( (fabs(errx) <= MIN_NEWT_TOL) && (fabs(errx) > NEWT_TOL) ){
    return(0);
  }
  else {
    return(1);
  }

  return(0);

}

/**********************************************************************/
/*********************************************************************************
   func_1d_orig():

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

static void func_1d_orig(CCTK_REAL x[], CCTK_REAL dx[], CCTK_REAL resid[],
                         CCTK_REAL jac[][NEWT_DIM], CCTK_REAL *f, CCTK_REAL *df, int n,
                         eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL dummy)
{

  // Set W from input
  CCTK_REAL W = x[0];

  // W^2
  CCTK_REAL Wsq = W*W;

  // v^2
  CCTK_REAL vsq = vsq_calc_1d(harm_aux,W);

  // Make sure that v^2 is physically reasonable, and if not make it so:
  const CCTK_REAL dv = 1.0e-15;
  vsq = ( vsq < -dv ) ?  0.      : fabs(vsq);
  vsq = ( vsq > 1. )  ?  (1.-dv) : vsq;

  // Compute P from W and v^2
  CCTK_REAL p_tmp = pressure_W_vsq_1d(eos,harm_aux,W,vsq);

  // Jacobian is calculated using full differentiation w.r.t. W
  CCTK_REAL dvsq = dvsq_dW_1d( harm_aux,W );
  CCTK_REAL dp1 = dpdW_calc_vsq_1d( W, vsq );
  CCTK_REAL dp2 = dpdvsq_calc_1d( eos, harm_aux, W, vsq );
  CCTK_REAL dpdW = dp1  + dp2*dvsq;

  // Compute the residual and the needed Jacobian component
  resid[0]  = W + 0.5 * harm_aux.Bsq * ( 1. + vsq ) - 0.5*harm_aux.QdotBsq/Wsq + harm_aux.Qdotn - p_tmp;
  jac[0][0] = 1. - dpdW + harm_aux.QdotBsq/(Wsq*W) + 0.5*harm_aux.Bsq*dvsq;
  // Set dx (NR step), f, and df (see function description above)
  dx[0] = - resid[0]/jac[0][0];
  *df   = - resid[0]*resid[0];
  *f    = - 0.5*(*df);

}

/**********************************************************************/
/****************************************************************************
   vsq_calc_1d():

      -- evaluate v^2 (spatial, normalized velocity) from
            W = \gamma^2 w

****************************************************************************/
static CCTK_REAL vsq_calc_1d(harm_aux_vars_struct harm_aux, CCTK_REAL W)
{

  CCTK_REAL Wsq = W*W;
  CCTK_REAL Xsq = (harm_aux.Bsq + W) * (harm_aux.Bsq + W);

  return( ( Wsq * harm_aux.Qtsq + harm_aux.QdotBsq * (harm_aux.Bsq + 2.*W)) / (Wsq*Xsq) );
}

/**********************************************************************/
/****************************************************************************
   dvsq_dW_1d():

      -- evaluate the partial derivative of v^2 w.r.t. W

****************************************************************************/
static CCTK_REAL dvsq_dW_1d(harm_aux_vars_struct harm_aux, CCTK_REAL W)
{

  CCTK_REAL X  = harm_aux.Bsq + W;
  CCTK_REAL W3 = W*W*W;
  CCTK_REAL X3 = X*X*X;

  return( -2.*( harm_aux.Qtsq/X3 + harm_aux.QdotBsq * (3*W*X + harm_aux.Bsq*harm_aux.Bsq) / ( W3 * X3 ) )  );
}

/**********************************************************************
 **********************************************************************

 The following routines specify the equation of state.  All routines
  above here should be indpendent of EOS.  If the user wishes
  to use another equation of state, the below functions must be replaced
  by equivalent routines based upon the new EOS.

  **********************************************************************
  **********************************************************************/

/**********************************************************************/
/**********************************************************************
  pressure_W_vsq():

        -- Gamma-law equation of state;
        -- pressure as a function of W, vsq, and D:
**********************************************************************/

static CCTK_REAL pressure_W_vsq_1d(eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W, CCTK_REAL vsq)
{

#if defined(ENABLE_STANDALONE_IGM_C2P_SOLVER) || defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#else
  DECLARE_CCTK_PARAMETERS;
#endif

  // Compute gamma^{-2} = 1 - v^{2} and gamma^{-1}
  CCTK_REAL inv_gammasq = 1.0 - vsq;
  CCTK_REAL inv_gamma   = sqrt(inv_gammasq);

  // Set D
  CCTK_REAL D = harm_aux.D;

  // Compute rho_b = D / gamma
  CCTK_REAL rho_b = D*inv_gamma;

  // Compute P_cold and eps_cold
  CCTK_REAL P_cold, eps_cold;
  compute_P_cold__eps_cold(eos,rho_b, P_cold,eps_cold);

  // Compute p = P_{cold} + P_{th}
  return( ( P_cold + (Gamma_th - 1.0)*( W*inv_gammasq - D*inv_gamma*( 1.0 + eps_cold ) ) )/Gamma_th );

}


/**********************************************************************/
/**********************************************************************
  dpdW_calc_vsq_1d():

      -- partial derivative of pressure with respect to W;
**********************************************************************/
static CCTK_REAL dpdW_calc_vsq_1d(CCTK_REAL W, CCTK_REAL vsq)
{

#if defined(ENABLE_STANDALONE_IGM_C2P_SOLVER) || defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#else
  DECLARE_CCTK_PARAMETERS;
#endif

  return( (Gamma_th - 1.0) * (1.0 - vsq) /  Gamma_th  ) ;

}

/**********************************************************************/
/**********************************************************************
  dpdvsq_calc_1d():

      -- partial derivative of pressure with respect to vsq
**********************************************************************/
static CCTK_REAL dpdvsq_calc_1d(eos_struct eos, harm_aux_vars_struct harm_aux, CCTK_REAL W, CCTK_REAL vsq)
{

  // This sets Gamma_th
#if defined(ENABLE_STANDALONE_IGM_C2P_SOLVER) || defined(ENABLE_STANDALONE_IGM_C2P_SOLVER_RANDOM_DATA)
#else
  DECLARE_CCTK_PARAMETERS;
#endif

  // Set D
  CCTK_REAL D = harm_aux.D;

  // Set gamma and rho
  CCTK_REAL gamma = 1.0/sqrt(1.0 - vsq);
  CCTK_REAL rho_b = D/gamma;

  // Compute P_cold and eps_cold
  CCTK_REAL P_cold, eps_cold;
  compute_P_cold__eps_cold(eos,rho_b, P_cold,eps_cold);

  // Set basic polytropic quantities
  int polytropic_index = find_polytropic_K_and_Gamma_index(eos,rho_b);
  CCTK_REAL Gamma_ppoly_tab = eos.Gamma_ppoly_tab[polytropic_index];

  /* Now we implement the derivative of P_cold with respect
   * to v^{2}, given by
   *  ----------------------------------------------------
   * | dP_cold/dvsq = gamma^{2 + Gamma_{poly}/2} P_{cold} |
   *  ----------------------------------------------------
   */
  CCTK_REAL dPcold_dvsq = P_cold * pow(gamma,2.0 + 0.5*Gamma_ppoly_tab);


  /* Now we implement the derivative of eps_cold with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------------
   * | deps_cold/dvsq = gamma/(D*(Gamma_ppoly_tab-1)) * (dP_cold/dvsq + gamma^{2} P_cold / 2) |
   *  -----------------------------------------------------------------------------------
   */
  CCTK_REAL depscold_dvsq = ( gamma/(D*(Gamma_ppoly_tab-1.0)) ) * ( dPcold_dvsq + 0.5*gamma*gamma*P_cold );

  /* Now we implement the derivative of p_hybrid with respect
   * to v^{2}, given by
   *  -----------------------------------------------------------------------------
   * | dp/dvsq = Gamma_th^{-1}( dP_cold/dvsq                                       |
   * |                          + (Gamma_{th}-1)*(-W                               |
   * |                                            + D gamma (1 + eps_cold)/2       |
   * |                                            - (D/gamma) * deps_cold/dvsq) )  |
   *  -----------------------------------------------------------------------------
   */
  return( ( dPcold_dvsq + (Gamma_th-1.0)*( -W + D*gamma*(1+eps_cold)/2.0 - D*depscold_dvsq/gamma ) )/Gamma_th );
}

/******************************************************************************
             END   OF   UTOPRIM_1D.C
******************************************************************************/
#endif // __HARM_UTOPRIM_1D__C__
