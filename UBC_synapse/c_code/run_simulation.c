#include <math.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////////////////////////
// UBC cleft glutamate diffusion model coupled to AMPA kinetic model (Raman and Trussell) //                                                                //
////////////////////////////////////////////////////////////////////////////////////////////

// Note: the element [i,j] of a matrix in R (i=row, j=col) corresponds to element [(j-1)*nrows + (i-1)] in C

void check_import()
{
    printf("c_code successfully imported\n");
}

double sum(double v[], int n)
{
  double s=0;
  for(int i=0;i<n;i++)
    s+=v[i];
  return s;
}


void sim(int *ig,           // identity matrix for intraglomerular space (1 in glomerulus, 0 outside)
             int *sy,         // identity matrix for synaptic spaces (0s and 1s)
             int nr,      // number of rows in c1 and c2
             int nc,      // number of column in c1 and c2
             int nit,        // number of iterations, i.e. length(spat)
             int *spat,       // synaptic release pattern (0s and 1s)
             double tstep,   // time step of stimulation in seconds
             double q,   // quantum of glutamate released at each synaptic site (in mM)
             double qi,     // equivalent of coefficient of diffusion inside glomerular space
             double qo,    // equivalent of coefficient of diffusion at boundary with extraglomerular space
             int *coords,      // coordinates of points of interest for continuous mesures
             int n,      // number of points of interest for continuous mesures
             double n_AMPA,  // number of AMPA channels per pixel
             double e_AMPA,  // inversion potential of AMPA channel
             double g_AMPA,  // conductance of AMPA channel
             double rm,       // membrane resistance
             double cm,       // membrane capacitance
             double v0,      // membrane resting potential
             double *resglu,  // variable to store the results of concentrations to the interet points
             double *resAMPAtot,  // variable to store the results of the overall open probability at all pixels
             double *resV     // vector to store the membrane potential variations in time
)
{
  ///////////////////////////////
  // DECLARATIONS OF VARIABLES //
  ///////////////////////////////
  // Declare integers and dereference pointers given as arguments
  int i,j,k,l,m ;
  double vm, vmw, i_tot ;
  double Ientr[nr*nc] ;
  // Create doubles for glutamate concentration matrices
  double c[nr*nc],cw[nr*nc],ca[nr*nc],cb[nr*nc],cc[nr*nc],cd[nr*nc];
  // Create doubles for identity matrices
  int iga[nr*nc],igb[nr*nc],igc[nr*nc],igd[nr*nc];
  // Create doubles for AMPA state matrices
  double C0[nr*nc],C1[nr*nc],C2[nr*nc],O2f[nr*nc],O2s[nr*nc],D1[nr*nc],D2[nr*nc],C3[nr*nc],O3[nr*nc];
  // Working copies
  double C0w[nr*nc],C1w[nr*nc],C2w[nr*nc],O2fw[nr*nc],O2sw[nr*nc],D1w[nr*nc],D2w[nr*nc],C3w[nr*nc],O3w[nr*nc];
  // Rate constants
  double kC0C1,kC1C0,kC1C2,kC2C1,kC2O2f,kO2fC2,kC2O2s,kO2sC2,kC1D1,kD1C1,kC2D2,kD2C2,kD1D2,kD2D1,kD2C3,kC3D2,kC3O3,kO3C3;
  double ts=tstep*1000;   // convert to ms
  kC0C1=ts*30;    // ! depends on glutamate concentration
  kC1C0=ts*0.3;
  kC1C2=ts*20;    // ! depends on glutamate concentration
  kC2C1=ts*600;
  kC2O2f=ts*60;
  kO2fC2=ts*3;
  kC2O2s=ts*3;
  kO2sC2=ts*0.35;
  kC1D1=ts*1;
  kD1C1=ts*0.3;
  kC2D2=ts*27;
  kD2C2=ts*0.014;
  kD1D2=ts*20;    // ! depends on glutamate concentration
  kD2D1=ts*1.038;
  kD2C3=ts*3.33;  // ! depends on glutamate concentration
  kC3D2=ts*0.22;
  kC3O3=ts*0.006;
  kO3C3=ts*2;


  /////////////////////
  // INITIALISATIONS //
  /////////////////////
  // Initialize Vm
  vm=v0 ; vmw=v0 ;
  // Initialize Ientr
  for(i=0;i<nr*nc;i++)
    Ientr[i]=0;
  // Initialize c and cw (working copy of c)
  for(i=0;i<nr*nc;i++)
    c[i]=0; cw[i]=0;
    // Initialize working copies of ig
    // iga: version of ig in which first row is deleted (elements moved one row up)
  for(l=2;l<(nr+1);l++){
    for(m=1;m<(nc+1);m++){
      iga[(m-1)*nr+(l-2)]=ig[(m-1)*nr+(l-1)];
    }
  }
  for(m=1;m<(nc+1);m++)
    iga[(m-1)*nr+(nr-1)]=0;
  // igb: version of ig in which last row is deleted (elements moved one row down)
  for(l=1;l<nr;l++){
    for(m=1;m<(nc+1);m++){
      igb[(m-1)*nr+l]=ig[(m-1)*nr+(l-1)];
    }
  }
  for(m=1;m<(nc+1);m++)
    igb[(m-1)*nr]=0;
  // igc: version of ig in which first column is deleted (elements moved one row left)
  for(l=1;l<(nr+1);l++){
    for(m=2;m<(nc+1);m++){
      igc[(m-2)*nr+(l-1)]=ig[(m-1)*nr+(l-1)];
    }
  }
  for(l=1;l<(nr+1);l++)
    igc[(nc-1)*nr+(l-1)]=0;
  // igd: version of ig in which last column is deleted (elements moved one row right)
  for(l=1;l<(nr+1);l++){
    for(m=1;m<nc;m++){
      igd[m*nr+(l-1)]=ig[(m-1)*nr+(l-1)];
    }
  }
  for(l=1;l<(nr+1);l++)
    igd[(l-1)]=0;

  // Initialize AMPA state matrices
  for(i=0;i<nr*nc;i++){
    C0[i]=1; C1[i]=0; C2[i]=0; O2f[i]=0; O2s[i]=0; D1[i]=0; D2[i]=0; C3[i]=0; O3[i]=0;
  }


  ////////////////
  // SIMULATION //
  ////////////////
  for(k=0;k<nit;k++){
    // Calculate working copies of c
    // ca: version of c in which first row is deleted (elements moved one row up)
    for(l=2;l<(nr+1);l++){
      for(m=1;m<(nc+1);m++){
        ca[(m-1)*nr+(l-2)]=c[(m-1)*nr+(l-1)];
      }
    }
    for(m=1;m<(nc+1);m++){
      ca[(m-1)*nr+(nr-1)]=0;
    }
    // cb: version of c in which last row is deleted (elements moved one row down)
    for(l=1;l<nr;l++){
      for(m=1;m<(nc+1);m++){
        cb[(m-1)*nr+l]=c[(m-1)*nr+(l-1)];
      }
    }
    for(m=1;m<(nc+1);m++){
      cb[(m-1)*nr]=0;
    }
    // cc: version of c in which first column is deleted (elements moved one row left)
    for(l=1;l<(nr+1);l++){
      for(m=2;m<(nc+1);m++){
        cc[(m-2)*nr+(l-1)]=c[(m-1)*nr+(l-1)];
      }
    }
    for(l=1;l<(nr+1);l++){
      cc[(nc-1)*nr+(l-1)]=0;
    }
    // cd: version of c in which last column is deleted (elements moved one row right)
    for(l=1;l<(nr+1);l++){
      for(m=1;m<nc;m++){
        cd[m*nr+(l-1)]=c[(m-1)*nr+(l-1)];
      }
    }
    for(l=1;l<(nr+1);l++){
      cd[(l-1)]=0;
    }
    // Core computation
    for(i=0;i<nr;i++){
      for(j=0;j<nc;j++){
        // Diffusion
        cw[j*nr+i] = c[j*nr+i]+ // glutamate concentration in the (i,j) "pixel" at previous time step
          q*spat[k]*sy[j*nr+i]+ // if synaptic site, add synaptically-released glutamate
          ig[j*nr+i]*(qi*(iga[j*nr+i]*(ca[j*nr+i]-c[j*nr+i])+ // diffusion with pixel below
          igb[j*nr+i]*(cb[j*nr+i]-c[j*nr+i])+ // diffusion with pixel above
          igc[j*nr+i]*(cc[j*nr+i]-c[j*nr+i])+ // diffusion with pixel to the right
          igd[j*nr+i]*(cd[j*nr+i]-c[j*nr+i]))- // diffusion with pixel to the left
          qo*(4-iga[j*nr+i]-igb[j*nr+i]-igc[j*nr+i]-igd[j*nr+i])*c[j*nr+i]); // diffusion outside glomerulus
        // Calculate AMPA state probabilities
        C0w[j*nr+i]=C0[j*nr+i]*(1-kC0C1*cw[j*nr+i])+C1[j*nr+i]*kC1C0;
        C1w[j*nr+i]=C1[j*nr+i]*(1-kC1C0-kC1D1-kC1C2*cw[j*nr+i])+C0[j*nr+i]*kC0C1*cw[j*nr+i]+D1[j*nr+i]*kD1C1+C2[j*nr+i]*kC2C1;
        C2w[j*nr+i]=C2[j*nr+i]*(1-kC2C1-kC2D2-kC2O2s-kC2O2f)+C1[j*nr+i]*kC1C2*cw[j*nr+i]+
          D2[j*nr+i]*kD2C2+O2s[j*nr+i]*kO2sC2+O2f[j*nr+i]*kO2fC2;
        O2fw[j*nr+i]=O2f[j*nr+i]*(1-kO2fC2)+C2[j*nr+i]*kC2O2f;
        O2sw[j*nr+i]=O2s[j*nr+i]*(1-kO2sC2)+C2[j*nr+i]*kC2O2s;
        D1w[j*nr+i]=D1[j*nr+i]*(1-kD1C1-kD1D2*cw[j*nr+i])+C1[j*nr+i]*kC1D1+D2[j*nr+i]*kD2D1;
        D2w[j*nr+i]=D2[j*nr+i]*(1-kD2D1-kD2C2-kD2C3*cw[j*nr+i])+D1[j*nr+i]*kD1D2*cw[j*nr+i]+C2[j*nr+i]*kC2D2+C3[j*nr+i]*kC3D2;
        C3w[j*nr+i]=C3[j*nr+i]*(1-kC3D2-kC3O3)+D2[j*nr+i]*kD2C3*cw[j*nr+i]+O3[j*nr+i]*kO3C3;
        O3w[j*nr+i]=O3[j*nr+i]*(1-kO3C3)+C3[j*nr+i]*kC3O3;
        // Compute the entering current
        // only at synaptic junctions, where there are AMPA receptors
        Ientr[j*nr+i] = (O2fw[j*nr+i]+O2sw[j*nr+i]+O3w[j*nr+i])*sy[j*nr+i]*n_AMPA*g_AMPA*(vm - e_AMPA) ;
      }
    }
    // Copy cw into c for next iteration
    for(i=0;i<nr*nc;i++)
      c[i]=cw[i];
    // Copy state matrices for next iteration
    for(i=0;i<nr*nc;i++){
      C0[i]=C0w[i]; C1[i]=C1w[i]; C2[i]=C2w[i]; O2f[i]=O2fw[i]; O2s[i]=O2sw[i]; D1[i]=D1w[i]; D2[i]=D2w[i]; C3[i]=C3w[i]; O3[i]=O3w[i];
    }
    // Compute the membrane voltage
    i_tot = sum(Ientr,nr*nc) ;
    vmw = vm - (i_tot + (vm - v0)/rm)*(ts/cm)/1000  ;
    vm = vmw ;
    // Return membrane voltage
    resV[k]=vm ;
    // Output the sum of all open state probabilities (over all "pixels")
    resAMPAtot[k]=sum(O2s,nr*nc)+sum(O2f,nr*nc)+sum(O3,nr*nc);
    // Return sum of glutamate concentration values at set of points of interest - continuous measures
    for(i=0;i<n;i++){
      l=coords[i];
      m=coords[i+n];
      resglu[k*n+i] = c[m*nr+l];
    }
  }
}

  