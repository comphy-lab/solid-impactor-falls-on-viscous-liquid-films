/** Title: Elastic impactor hitting a vicous liquid film.
# Version: 1.0
# Updated: Dec 5, 2024

# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

# Main feature: 
- 2D+axi viscoelastic scalar implementation of the spherical impactor falling on a viscous liquid film.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "../src-local/log-conform-viscoelastic-scalar-2D.h"
#define FILTERED
#include "../src-local/two-phase-TF-VE.h"
#include "tension.h"

// error tolerances
#define fErr (1e-3)
#define VelErr (1e-3)
#define KErr (1e-3)
#define AErr (1e-3)

// Calculations!
#define Xdist (1.040)
#define R2Drop(x,y,z) (sq(x - Xdist) + sq(y))
// domain
#define Ldomain 4                                // Dimension of the domain

// boundary conditions
u.t[left] = dirichlet(0.0);
f1[left] = 0.0;
f2[left] = 1.0;

int MAXlevel;
double tmax, We, Ohd, Bo, Ohf, hf, Ec, De;
#define MINlevel 4                                            // minimum level
#define tsnap (0.01)

char comm[80], restartFile[80], logFile[80];
int main(int argc, char const *argv[]) {

  MAXlevel = 10;

  tmax = 2e0;

  We = 1e1; // We is 1 for 0.167 m/s <816*0.167^2*0.00075/0.017>
  Ohd = 1e-2; // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  Ec = 1e0;
  De = 1e30;
  Ohf = 1e0;
  hf = 0.10;

  if (hf == 0){
    fprintf(ferr, "We have a problem. Wrong code. Change code or film height.\n");
    return 1;
  }
  fprintf(ferr, "Level %d tmax %g. We %g, Ohd %g, Ec %g, De 1e30, Ohf %3.2e, hf %3.2f\n", MAXlevel, tmax, We, Ohd, Ec, Ohf, hf);

  L0=Ldomain;
  X0=-hf; Y0=0.;
  init_grid (1 << (MINlevel));

  rho1 = 1.0; mu1 = Ohd/sqrt(We); 
  G1 = Ec/We; lambda1 = De*sqrt(We);
  
  rho2 = 1.0; mu2 = Ohf/sqrt(We);
  G2 = 0.0; lambda2 = 0.0;

  rho3 = 1e-3; mu3 = 1e-4/sqrt(We); // to run things fast for now, we are not using the right density and viscosity ratio. It is advisable to use the navier-stokes/conserving.h for the correct air-water density ratio -- #TODO: fixme for production runs. 
  G3 = 0.0; lambda3 = 0.0;

  f1.sigma = 1e2/We; f2.sigma = 1.0/We;

  /*
  In the folder called "intermediate", we will store the snapshot files as the simulation progresses. See event: writingFiles.
  */
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  /*
  We save the restart file every t = tsnap such that if the simulation is interrupted, we can recover from the last saved state.
  */
  sprintf (restartFile, "restart");
  /*
  We also save the log file every few time steps.
  */
  sprintf (logFile, "logFile.dat");

  run();

}

event init(t = 0){
  if(!restore (file = "dump")){
    refine((R2Drop(x,y,z) < 1.44) && (level < MAXlevel));
    fraction (f1, 1. - R2Drop(x,y,z));
    fraction (f2, -(x-1e-2));
    foreach () {
      u.x[] = -1e0*f1[];
      u.y[] = 0.0;
    }
  }
}


event adapt(i++) {
  /*
  We adapt the grid to resolve the interface and the flow. For infinite Weber number, these should be fine. 
  For finite Weber number, perhaps adaptation based on curvature is also needed. 
  */
 scalar KAPPA1[], KAPPA2[];
 curvature (f1, KAPPA1); curvature (f2, KAPPA2);
  adapt_wavelet ((scalar *){f1, f2, u.x, u.y, A11, A12, A22, KAPPA1, KAPPA2},
    (double[]){fErr, fErr, VelErr, VelErr, AErr, AErr, AErr, KErr, KErr},
    MAXlevel, 4);
}

// Outputs
event writingFiles (t = 0; t += tsnap; t <= tmax + tsnap) {
  dump (file = restartFile);
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

event logWriting (i++) {

  double ke = 0.;
  foreach (reduction(+:ke)){
    ke += 0.5*(2*pi*y)*(f1[]+f2[])*(sq(u.y[])+sq(u.x[]))*sq(Delta);
  }

  if (pid() == 0){
    static FILE * fp;
    if (i == 0) {
      fprintf (ferr, "i dt t ke rmax\n");
      fp = fopen (logFile, "w");
      fprintf(fp, "Level %d tmax %g. We %g, Ohd %g, Ec %g, De 1e30, Ohf %3.2e, hf %3.2f\n", MAXlevel, tmax, We, Ohd, Ec, Ohf, hf);
      fprintf (fp, "i dt t ke\n");
    } else {
      fp = fopen (logFile, "a");
    }

    fprintf (fp, "%d %g %g %g\n", i, dt, t, ke);
    fclose(fp);
    fprintf (ferr, "%d %g %g %g\n", i, dt, t, ke);
    
    if (ke < -1e-10){
      dump(file="FailedDumpKeNegative");
      fprintf(ferr, "Something very bad is happening! Stopping!\n");
      return 1;
    }
    if (ke < 1e-10 && i > 100){
      dump(file="StoppedDump");
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
    if (ke > 1e4){
      dump(file="FailedDump");
      fprintf(ferr, "kinetic energy blew up! Stopping!\n");
      fp = fopen (logFile, "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
  }

}

event end (t = tmax + tsnap) {
  dump(file="FinalDump");

  static FILE * fp;
  fp = fopen (logFile, "a");
  fprintf(fp, "Level %d tmax %g. We %g, Ohd %g, Ec %g, De 1e30, Ohf %3.2e, hf %3.2f\n", MAXlevel, tmax, We, Ohd, Ec, Ohf, hf);
  fprintf(fp, "End of simulation.\n");
  fclose(fp);
}