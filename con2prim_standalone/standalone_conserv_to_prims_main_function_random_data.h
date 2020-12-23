

#define CCTK_REAL double
#define CCTK_EQUALS(a,b) (strcmp((a),(b))==0)
#define CCTK_Equals(a,b) (strcmp((a),(b))==0)
#define cGH int

#define RANDOM_NUMBER_BET_m1_and_p1 (2.0*(double)rand()/(double)RAND_MAX - 1.0)

#define CCTK_isnan std::isnan

int conserv_to_prims_debug = 1;
char verbose[100];

CCTK_REAL GAMMA_SPEED_LIMIT,rho_b_atm,tau_atm, rho_b_max, Psi6threshold;

CCTK_REAL Gamma_th, K_ppoly_tab0;
CCTK_REAL rho_ppoly_tab_in[10],Gamma_ppoly_tab_in[11];

int neos;
int update_Tmunu;

#define CCTK_THORNSTRING ""
#define CCTK_WARN_ALERT  ""
#define Symmetry "none"
#define CCTK_GFINDEX3D(IGNORE,i,j,k) ((i) + cctk_lsh[0] * ((j) + cctk_lsh[1] * (k)))
#define GetRefinementLevel(cctkGH) 0

#define LOOP_ALL_GRIDPOINTS for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++)

#include <stdarg.h>
#include <string.h>
#include "IllinoisGRMHD_headers.h"
#include "harm_primitives_headers.h"
#include "harm_u2p_util.c"
#include "inlined_functions.C"
#include "apply_tau_floor__enforce_limits_on_primitives_and_recompute_conservs.C"
#include "convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij.C"

int CCTK_VInfo(const char *thorn, const char *format, ...) {
  va_list ap;
  fprintf (stdout, "INFO (NOTHORN): ");
  va_start (ap, format);
  vfprintf (stdout, format, ap);
  va_end (ap);
  fprintf (stdout, "\n");
  return 0;
}

int *cctkGH;

int main(int argc, const char *argv[]) {

#ifdef WRITE_TEST_DATA 
  int write_data=1;
  //Select 1 or 0 for variations in the output
  int write_tmunudn=1;
  int write_tmunuup=1;
  int write_bssn=1;
  int write_metric4dn=1;
  int write_metric4up=1;
  //int write_ref_prims=1;
#else
  int write_data=0;
#endif
  
  if( (argc < 2+write_data) || (argc > 3+write_data) ) {
    std::cerr<<"Error: Correct usage: ./C2P_standalone_random_data number_of_test_points"<<(write_data?" output_file_name":"")<<" random_seed(optional, default=0)\n"<<std::endl;
    exit(1);
  }

  sprintf(verbose,"essential+iteration output");
  // We use proper C++ here, for file I/O later.
  using namespace std;

  // Set parameters needed by IGM
  // These mimic the values found in the parameter
  // file TOV-NRPyPlusTOVID.par
  GAMMA_SPEED_LIMIT = 10.0;
  tau_atm = 4.876083025795607e-12;
  rho_b_atm = 1.292852735094440e-10;
  rho_b_max = 1e300;
  Psi6threshold = 1e100;

  // These are simple polytrope EOS parameters
  // that can also be found in TOV-NRPyPlusTOVID.par
  neos = 1;
  Gamma_th = 2.0;
  K_ppoly_tab0 = 1.0;
  rho_ppoly_tab_in[0] = 0.0;
  Gamma_ppoly_tab_in[0] = 2.0;

  // Not completely sure about this one
  update_Tmunu = 1;

  // Set the local shape (i.e. Nx,Ny,Nz)
  int cctk_lsh[3] = {(int)(atof(argv[1])+0.5),1,1};

  // Set fullsize of arrays
  int fullsize=cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];

  // Read in the domain size, which is arbitrary
  CCTK_REAL domain_size = 10.0;

  // Allocate memory for (x,y,z)
  CCTK_REAL *x = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *y = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *z = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);

  // Compute (dx,dy,dz)
  CCTK_REAL dx = 2.0*domain_size/cctk_lsh[0];
  CCTK_REAL dy = 2.0*domain_size/cctk_lsh[1];
  CCTK_REAL dz = 2.0*domain_size/cctk_lsh[2];

  // Set (x,y,z)
  for(int i=0;i<cctk_lsh[0];i++) x[i] = -domain_size + dx*i;
  for(int j=0;j<cctk_lsh[1];j++) y[j] = -domain_size + dy*j;
  for(int k=0;k<cctk_lsh[2];k++) z[k] = -domain_size + dz*k;

  // Not sure what to do with these yet
  CCTK_REAL *failure_checker = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTtt = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTtx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTty = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTtz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTxx = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTxy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTxz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTyy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTyz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *eTzz = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);

  // Allocate memory for the ADM metric. We will initialize the
  // ADM metric to a perturbation around flat space, for which
  // alpha = 1, beta^{i} = 0, and gamma_{ij} = diag(1,1,1).
  CCTK_REAL *alp       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gxx       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gxy       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gxz       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gyy       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gyz       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gzz       = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *betax     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *betay     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *betaz     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);

  // Primitive variables: We will initialize these to random values as well
  CCTK_REAL *P         = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *rho_b     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *Bx        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *By        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *Bz        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *vx        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *vy        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *vz        = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
#if( USE_ENTROPY_EQUATION )
  CCTK_REAL *S_entropy = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
#endif

#ifdef WRITE_TEST_DATA
  //If set up for writing data, set the output file name
  ofstream outfile;
  if( write_data ){
    string outfilename(argv[2]);
    outfile.open(outfilename);
  }
  //Write header line
  outfile<<"#CONSERVS[5],FLUIDPRIMS[5],ADMMETRIC[10],BFIELD[3]";
  if( write_tmunudn ) outfile<<",TMNUDN[10]";
  if( write_tmunuup ) outfile<<",TMNUUP[10]";
  if( write_bssn ) outfile<<",METRICBSSN[18]";
  if( write_metric4dn ) outfile<<",METRIC4DN[10]";
  if( write_metric4up ) outfile<<",METRIC4UP[10]";
  //if( write_ref_prims ) outfile<<",MODEL_PRIMS[10]";
  outfile<<endl;
#endif

  // Set the random seed to ensure consistent
  // random number generation between different
  // runs of the code
  int random_seed = 0;
  if( argc == 3+write_data ) random_seed = atoi(argv[2+write_data]);
  srand(random_seed);

  

  // Set the pertubation. Note that we will perturb
  // the metric quantities by a random value in the
  // range [-perturbation,+perturbation]
  CCTK_REAL perturbation     = 1e-2;
  CCTK_REAL avg_absBx        = 0.0;
  CCTK_REAL avg_absBy        = 0.0;
  CCTK_REAL avg_absBz        = 0.0;
  CCTK_REAL avg_Bsq          = 0.0;
  CCTK_REAL avg_Bsq_over_rho = 0.0;
  CCTK_REAL avg_Bsq_over_P   = 0.0;
  CCTK_REAL avg_vsq          = 0.0;
  CCTK_REAL avg_absvx        = 0.0;
  CCTK_REAL avg_absvy        = 0.0;
  CCTK_REAL avg_absvz        = 0.0;
  CCTK_REAL avg_gamma        = 0.0;
  // Initialize the metric and primitives
  LOOP_ALL_GRIDPOINTS {
    int index    = CCTK_GFINDEX3D(cctkGH,i,j,k); // cctkGH is ignored here
    alp[index]   = 1.0 + 50.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    gxx[index]   = 1.0 + 50.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    gxy[index]   = 0.0 + 10.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    gxz[index]   = 0.0 + 10.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    gyy[index]   = 1.0 + 50.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    gyz[index]   = 0.0 + 10.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    gzz[index]   = 1.0 + 10.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    betax[index] = 0.0 + 50.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    betay[index] = 0.0 + 50.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    betaz[index] = 0.0 + 50.0*perturbation*RANDOM_NUMBER_BET_m1_and_p1;

    rho_b[index] = 0.9 + perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    // P[index]     = K_ppoly_tab0 * pow(rho_b[index],Gamma_ppoly_tab_in[0]);
    // Different test: Let us randomly set epsilon, the internal energy
    // density, and compute the pressure using:
    // P = (Gamma-1)*rho*epsilon
    CCTK_REAL eps = 0.9 + perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    P[index]     = (Gamma_ppoly_tab_in[0]-1.0) * rho_b[index] * eps;

    Bx[index]    = 2.5*RANDOM_NUMBER_BET_m1_and_p1;
    By[index]    = 2.5*RANDOM_NUMBER_BET_m1_and_p1;
    Bz[index]    = 2.5*RANDOM_NUMBER_BET_m1_and_p1;
    vx[index]    = 0.25*RANDOM_NUMBER_BET_m1_and_p1;
    vy[index]    = 0.25*RANDOM_NUMBER_BET_m1_and_p1;
    vz[index]    = 0.25*RANDOM_NUMBER_BET_m1_and_p1;
    // Higher velocity test (aiming for v^2 ~ 0.65)
    // vx[index]    = 0.47 + perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    // vy[index]    = 0.47 + perturbation*RANDOM_NUMBER_BET_m1_and_p1;
    // vz[index]    = 0.47 + perturbation*RANDOM_NUMBER_BET_m1_and_p1;

    avg_absBx        += fabs(Bx[index]);
    avg_absBy        += fabs(By[index]);
    avg_absBz        += fabs(Bz[index]);
    CCTK_REAL Bsq     = Bx[index]*Bx[index] + By[index]*By[index] + Bz[index]*Bz[index];
    avg_Bsq          += Bsq;
    avg_Bsq_over_P   += Bsq/P[index];
    avg_Bsq_over_rho += Bsq/rho_b[index];

    avg_absvx        += fabs(vx[index]);
    avg_absvy        += fabs(vy[index]);
    avg_absvz        += fabs(vz[index]);
    CCTK_REAL vsq     = vx[index]*vx[index] + vy[index]*vy[index] + vz[index]*vz[index];
    avg_vsq          += vsq;
    avg_gamma        += 1.0/sqrt(1.0 - vsq);

#if( USE_ENTROPY_EQUATION )
    S_entropy[index] = P[index] * pow(rho_b[index],1.0-Gamma_ppoly_tab_in[0]);
#endif
  }
  avg_absBx        /= fullsize;
  avg_absBy        /= fullsize;
  avg_absBz        /= fullsize;
  avg_Bsq          /= fullsize;
  avg_Bsq_over_P   /= fullsize;
  avg_Bsq_over_rho /= fullsize;

  avg_absvx        /= fullsize;
  avg_absvy        /= fullsize;
  avg_absvz        /= fullsize;
  avg_vsq          /= fullsize;
  avg_gamma        /= fullsize;

  printf("Averages:\n");
  printf("|Bx| = %.3e |By| = %.3e |Bz| = %.3e  Bsq = %.3e Bsq/P = %.3e Bsq/rho = %.3e\n",avg_absBx,avg_absBy,avg_absBz,avg_Bsq,avg_Bsq_over_P,avg_Bsq_over_rho);
  printf("|vx| = %.3e |vy| = %.3e |vz| = %.3e  vsq = %.3e gamma = %.3e\n",avg_absvx,avg_absvy,avg_absvz,avg_vsq,avg_gamma);

  // All BSSN quantities below, including lapm1, are initialized
  // by the C2P driver by calling the function:
  // IllinoisGRMHD_convert_ADM_to_BSSN__enforce_detgtij_eq_1__and_compute_gtupij()
  CCTK_REAL *psi_bssn = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *phi_bssn = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtxx     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtxy     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtxz     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtyy     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtyz     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtzz     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupxx   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupxy   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupxz   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupyy   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupyz   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *gtupzz   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *lapm1    = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);

  // Conservative variables: These we compute from the primitives,
  // but we do it inside the C2P driver since we need Psi6 to compute
  // sqrt(-g).
  CCTK_REAL *rho_star = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *tau      = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *mhd_st_x = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *mhd_st_y = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *mhd_st_z = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
#if( USE_ENTROPY_EQUATION )
  CCTK_REAL *S_star   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
#endif

  // HERE WE USE _flux variables as temp storage for original values of conservative variables..
  // This is used for debugging purposes only.
  CCTK_REAL *rho_star_flux = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *st_x_flux     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *st_y_flux     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *st_z_flux     = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
  CCTK_REAL *tau_flux      = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
#if( USE_ENTROPY_EQUATION )
  CCTK_REAL *S_star_flux   = (CCTK_REAL *)malloc(sizeof(CCTK_REAL)*fullsize);
#endif
