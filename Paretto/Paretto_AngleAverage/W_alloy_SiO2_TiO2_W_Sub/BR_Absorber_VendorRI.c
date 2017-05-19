#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<malloc.h>
#include"memory.h"
#include<complex.h>
#include<time.h>
// Function Prototypes
int *VEC_INT(int dim);
double *VEC_DOUBLE(int dim);
char *VEC_CHAR(int dim);
double complex *VEC_CDOUBLE(int dim);
void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C);
void TransferMatrix(double thetaI, double k0, double complex *rind, double *d,
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21);
void Bruggenman(double f, double epsD, double complex epsM, double *eta, double *kappa);
void MaxwellGarnett(double f, double epsD, double complex epsM, double *eta, double *kappa);
void Lorentz(double we, double de, double w, double *epsr, double *epsi);
int ReadDielectric(char *file, double *lambda, double complex *epsM);
int IsDominated(int idx, int LENGTH, double *O1,double *O2);
void ReadBRRind(int numBR, double lambda, double *BRlambda, double complex *BRind, double *n, double *k); 
// Gauss-Legendre headers
void legendre_compute_glr ( int n, double x[], double w[] );
void legendre_compute_glr0 ( int n, double *p, double *pp );
void legendre_compute_glr1 ( int n, double *roots, double *ders );
void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
void legendre_handle ( int n, double a, double b );
void r8mat_write ( char *output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
double rk2_leg ( double t, double tn, double x, int n );
void timestamp ( void );
double ts_mult ( double *u, double h, int n );
double wtime ( );



// Global variables
int Nlayer;
int polflag;
double c = 299792458;
double pi=3.141592653589793;

int main(int argc, char* argv[]) {
  // integer variables
  int i, j, k, F1, F2, TK, NI, NV;
  // complex double precision variables
  double complex m11P, m21P, rP, tP, cosLP;
  double complex m11S, m21S, rS, tS, cosLS;

  double complex *rind, nlow, nhi;
  double complex *sio2_rind, *tio2_rind, *w_sub, *w_ald, *alumina_ald;

  double h=6.626e-34;
  double kb = 1.38064852e-23;
  double rho;

  double SE, PU;

  // real double precision variables
  double RP, TP, AP, RS, TS, AS;
  double *d, thetaI, lambda, k0, alphaP, betaP, alphaS, betaS;
  double sti, n1, n2, thetaT, rp, Rp, Tangle;
  double eta, kappa;
  double we, de, w;
  double dalloy, d1, d2, d3, d4, vf1, vf2, epsbg, fac1, fac2;

  // Lists for spectral efficiency
  double *LamList, *EmissS, *EmissP, *clam;
  // Variables for Spectral Efficiency
  double Temp, lbg;
  // This is intentionally larger than number of wavelength values in data files
  int NumLam=10000;

  // For Legendre quadrature:
  int leg_N = 10;  // 10 points between 0 and 2pi
  double leg_a = 0.;
  double leg_b = pi/2.; 
  double *leg_x, *leg_w;
  
  // Allocate memory for array of Legendre roots and weights
  leg_x = (double *)malloc(leg_N*sizeof(double));
  leg_w = (double *)malloc(leg_N*sizeof(double));

  legendre_handle( leg_N, leg_a, leg_b );

  legendre_compute_glr(leg_N, leg_x, leg_w);
  
  rescale( leg_a, leg_b, leg_N, leg_x, leg_w);

  // test the grid
  for (int leg_i=0; leg_i<leg_N; leg_i++) {

    printf("  %20.16f  %20.16f\n",leg_x[leg_i],leg_w[leg_i]);

  }

  //  Allocate arrays for spectral efficiency
  LamList = VEC_DOUBLE(NumLam);
  EmissS  = VEC_DOUBLE(NumLam);
  EmissP  = VEC_DOUBLE(NumLam);
  clam    = VEC_DOUBLE(NumLam);

  //   Metal epsilon array
  sio2_rind   = VEC_CDOUBLE(NumLam);
  tio2_rind   = VEC_CDOUBLE(NumLam);
  w_sub       = VEC_CDOUBLE(NumLam);
  w_ald       = VEC_CDOUBLE(NumLam);
  alumina_ald = VEC_CDOUBLE(NumLam);

  FILE *fp;

  // Character string(s)
  char *write, *line, *sio2file, *tio2file, *w_sub_file, *w_ald_file, *alumina_ald_file;

  write   = VEC_CHAR(1000);
  line    = VEC_CHAR(1000);
  sio2file = VEC_CHAR(1000);
  tio2file = VEC_CHAR(1000);
  w_sub_file = VEC_CHAR(1000);
  w_ald_file = VEC_CHAR(1000);
  alumina_ald_file = VEC_CHAR(1000);
  strcpy(sio2file,"DIEL/sio2_cspline.txt");
  strcpy(tio2file,"DIEL/tio2_cspline.txt");
  strcpy(w_sub_file,"DIEL/W_Palik.txt");
  strcpy(w_ald_file,"DIEL/w_ald_cspline.txt");
  strcpy(alumina_ald_file,"DIEL/ald_al2o3_cspline.txt");
  int CheckNum;
  // How many data points are in the file W_Palik.txt?  This function  will tell us
  // Substrate W data - dielectric function
  NumLam =   ReadDielectric(w_sub_file, LamList, w_sub);
  // BR data - refractive index data from vendor
  CheckNum = ReadDielectric(sio2file, clam, sio2_rind);
  CheckNum = ReadDielectric(tio2file, clam, tio2_rind); 
  // W data - refractive index from ALD sample 
  CheckNum = ReadDielectric(w_ald_file, clam, w_ald);
  // Alumina ata - refractive index from ALD sample
  CheckNum = ReadDielectric(alumina_ald_file, clam, alumina_ald);  

  //  Did we pass a filename to the program?
  if (argc==1) {
    exit(0);
  }

  strcpy(write,argv[1]);

  // initialize variables to be read
  Nlayer=0;
  d1=0.0;
  d2=0.0;
  nlow=0.;
  nhi=0.;
  vf1=0.0;
  vf2=0.0;
  epsbg=0.0;
  d3=0.0;
  d4=0.0;
  lbg=0.0;
  Temp=0.0;
  // Open the file for writing!
  fp = fopen(write,"r");
  printf("  going to read from file %s\n",write);
  fflush(stdout);
 
  
    fscanf(fp,"%s",line);  // Nlayer
    fscanf(fp,"%i",&Nlayer);
    fscanf(fp,"%s",line);   // dalloy
    fscanf(fp,"%lf",&dalloy);
    fscanf(fp,"%s",line);   //  d1
    fscanf(fp,"%lf",&d1);
    fscanf(fp,"%s",line);  //  nlow
    fscanf(fp,"%lf",&nlow);
    fscanf(fp,"%s",line);  // d2
    fscanf(fp,"%lf",&d2);
    fscanf(fp,"%s",line);  // nhi
    fscanf(fp,"%lf",&nhi);
    fscanf(fp,"%s",line);  // d3 (spacer)
    fscanf(fp,"%lf",&d3);
    fscanf(fp,"%s",line);  // d4 (W/alloy underlayer)
    fscanf(fp,"%lf",&d4);
    fscanf(fp,"%s",line);  // volume fraction 1
    fscanf(fp,"%lf",&vf1); 
    fscanf(fp,"%s",line);  // volume fraction 2
    fscanf(fp,"%lf",&vf2); 
    fscanf(fp,"%s",line);  // epsbg
    fscanf(fp,"%lf",&epsbg);
    fscanf(fp,"%s",line);  // Temp
    fscanf(fp,"%lf",&Temp);  
    fscanf(fp,"%s",line);  //Lambda bg
    fscanf(fp,"%lf",&lbg);

  polflag=1;

  int numVf, numNlayers, numFac, *NLa, *PF, numT;
  double *VFa, *SEA, *SFAC, *SDA, *Tem;
  
  numNlayers=20;
  numFac=20;
  numVf = 20;

  NLa = (int*)malloc((numNlayers*sizeof(int)));
  SEA = (double*)malloc((numFac*numVf*numNlayers*numFac*sizeof(double)));
  SDA = (double*)malloc((numFac*numVf*numNlayers*numFac*sizeof(double)));
  PF  = (int*)malloc((numFac*numVf*numNlayers*numFac*sizeof(int)));
  SFAC= (double*)malloc((numFac*sizeof(double)));
  VFa = (double*)malloc((numVf*sizeof(double)));

  d = VEC_DOUBLE(1000);
  rind = VEC_CDOUBLE(1000);
 


  for (TK=0; TK<numVf; TK++) {

    vf1 = 0. + (1./20)*TK;
    VFa[TK] = vf1;
  
    // Loop over different factors to multiply d1 by
    for (F1=0; F1<numFac; F1++) {

      fac1 = 0.6 + (1./20)*F1;   
      SFAC[F1] = fac1;

      // Loop over different factors to multiply d2 by
      for (F2=0; F2<numFac; F2++) {

        fac2 = 0.6 + (1./20)*F2;

        //  Loop over different number of layers
        for (NI=0; NI<numNlayers; NI++) { 

          Nlayer = 5 + NI;

          // NLa vector stores the value of Nlayer_i for layter use
          NLa[NI] = Nlayer;

          // Vector to store the thicknesses in micrometers of each layer
          d[0] = 0.;

          d[1] = dalloy;

          // Refractive index of air
          rind[0] = 1.00 + 0.*I;

          rind[1] = 1.00 + 0.*I;

          // Now start the Bragg Reflector
          for (i=2; i<Nlayer-2; i++) {

            if (i%2==0) {
              d[i] = d1*fac1;
              rind[i] = nlow + 0.*I;
            }
            else {
              d[i] = d2*fac2;
              rind[i] = nhi + 0.*I;
            }
          }

          d[Nlayer-3] = 0.01;
          rind[Nlayer-3] = sqrt(epsbg) + 0.*I;
          // W layer that is the substrate for the Bragg Reflector
          d[Nlayer-2] = 0.9;
          // Temporary - will replace with Tungsten!
          rind[Nlayer-2] = 1.0 + 0.*I;
 
          // Air underneath
          d[Nlayer-1] = 0.;
          rind[Nlayer-1] = 1.0 + 0.*I;
 
 
         //  Top/Bottom layer RI for Transmission calculation
         n1 = creal(rind[0]);
         n2 = creal(rind[Nlayer-1]);
   
         // Normal incidence
         double PUP, PUS, SEP, SES;
         double SE_Den_sum, SE_Num_sum, BBsum;
         double dtheta = pi/90.;
         double dlambda;
         SE_Num_sum = 0.;
         SE_Den_sum = 0.;
         BBsum = 0.;

         for (k=0; k<leg_N; k++) {

           // thetaI = k*dtheta; -> regular grid
           // Gauss-Legendre grid
           thetaI = leg_x[k];
           dtheta = leg_w[k];

           double BB_theta=0;
           double BB_int = 0;
           double P_SE_Num = 0;
           double P_SE_Den = 0;
           double S_SE_Num = 0;
           double S_SE_Den = 0;

      
           for (i=0; i<NumLam; i++) {
 
             lambda = LamList[i];    // Lambda in meters
             if (i==0) {
               dlambda = fabs(LamList[1]-LamList[0]);
             }
             else if (i>=(NumLam-1)) {
               dlambda = fabs(LamList[NumLam-1] - LamList[NumLam-2]);
             }
             else {
               dlambda = fabs(LamList[i+1]-LamList[i]);
             }

             k0 = 2*pi*1e-6/lambda;  // k0 in inverse microns - verified
             w=2*pi*c/lambda;        // angular frequency 
 
             epsbg = creal(alumina_ald[i]*alumina_ald[i]);
             double complex epsald  = w_ald[i]*w_ald[i];
             // Alloy superstrate Layer (Layer 1 in the structure [Layer 0 is air!])
             //MaxwellGarnett(vf1, epsbg, epsald, &eta, &kappa);
             Bruggenman(vf1,epsbg, epsald, &eta, &kappa);
             rind[1] = eta + I*kappa; 

             // Now start the PC
             for (j=2; j<Nlayer-3; j++) {

               if (j%2==0) {
                 rind[j] = sio2_rind[i];
               }
               else {
                 rind[j] = tio2_rind[i];
               }
             }

             // Alumina layer
             rind[Nlayer-3] = alumina_ald[i];

             // W substrate layer (Layer N-2 in the structure [layer N-1 is air!])
             MaxwellGarnett(1.0, epsbg, w_sub[i], &eta, &kappa);
             rind[Nlayer-2] = eta + I*kappa;

             // p-polarized
             polflag = 1; 
             // Solve the Transfer Matrix Equations
             TransferMatrix(thetaI, k0, rind, d, &cosLP, &betaP, &alphaP, &m11P, &m21P);
             // s-polarized
             polflag = 2;
             TransferMatrix(thetaI, k0, rind, d, &cosLS, &betaS, &alphaS, &m11S, &m21S);
             // power density per unit solid angle, etc
             rho = (2*h*c*c/pow(lambda,5))*(1/(exp(h*c/(lambda*kb*Temp))-1));
 
             // Fresnel reflection and trans coefficient - p pol
             rP = m21P/m11P; 
             tP = 1./m11P;

             // Fresnel reflection and trans coefficient - s pol
             rS = m21S/m11S;
             tS = 1./m11S;

             // Reflectance, which is a real quantity between 0 and 1
             // p-pol
             RP = creal(rP*conj(rP));
             // s-pol
             RS = creal(rS*conj(rS));
             
             // p-pol T and A
             Tangle =  n2*creal(cosLP)/(n1*cos(thetaI));
             TP = creal(tP*conj(tP))*Tangle;
             AP = 1 - RP - TP;
 
             // s-pol T and A
             Tangle = n2*creal(cosLS)/(n1*cos(thetaI));
             TS = creal(tS*conj(tS))*Tangle;
             AS = 1 - RS - TS;

             // Store absorbance/emissivity in array Emiss
             EmissP[i] = AP;
             EmissS[i] = AS;
      
             BB_theta += rho*cos(thetaI)*dlambda;
             BB_int   += pi*rho*dlambda;

             if (lambda<=lbg) {
               P_SE_Num += (lambda/lbg)*AP*rho*cos(thetaI)*dlambda;
               S_SE_Num += (lambda/lbg)*AS*rho*cos(thetaI)*dlambda;
             }
             // Denomenator of Spectral Efficiency
             P_SE_Den += AP*rho*cos(thetaI)*dlambda;
             S_SE_Den += AS*rho*cos(thetaI)*dlambda;
             //printf(" ( %12.10f) ( %12.10f) %12.10f  %12.10f\n",BB_theta-P_SE_Den,BB_theta-S_SE_Num,AP,AS);

           }
           // Sum of sin(theta) dtheta - the factor of 2pi comes from the
           // integral over phi - no phi dependence so this just gives a factor of 2pi
           SE_Num_sum += 2*pi*(P_SE_Num/2. + S_SE_Num/2.)*sin(thetaI)*dtheta;
           SE_Den_sum += 2*pi*(P_SE_Den/2. + S_SE_Den/2.)*sin(thetaI)*dtheta;
           BBsum += 2*pi*BB_theta*sin(thetaI)*dtheta;
           //printf("  %12.10f  %12.10f  %12.10f %12.10f\n",SE_Num_sum, SE_Den_sum, BBsum,  BB_int);

 
         }
   
         SE = SE_Num_sum/SE_Den_sum;
         PU = SE_Num_sum; 
         printf("  %f  %f  %f    %8.6f    %8.6f    %8.6f %f    %i     %12.10f  %12.10e\n",
                 dalloy,fac1, fac2, d1*fac1, d2*fac2, vf1,  Temp, Nlayer, SE,     PU);
         SEA[TK*numVf*numFac*numFac+F1*numFac*numFac+F2*numFac+NI] = SE;
         SDA[TK*numVf*numFac*numFac+F1*numFac*numFac+F2*numFac+NI] = PU;

       }
     }
   }
 }


  FILE *pf;
  pf = fopen("paretto_Bruggenman_GaussLegendre_Grid.txt","w");
  int id;

  for (TK=0; TK<numVf; TK++) {

  for (F1=0; F1<numFac; F1++) {

  for (F2=0; F2<numFac; F2++) {
  
  for (NI=0; NI<numNlayers; NI++) {

    i = TK*numVf*numFac*numFac+F1*numFac*numFac+F2*numFac+NI;
    // id is 1 if member i is dominated by at least one other member j!=i
    // a member is pareto optimal only if it is NOT dominated 
    id = IsDominated(i, numVf*numFac*numNlayers*numFac, SEA, SDA);
    if (id) PF[i] = 0;
    else {
      PF[i] = 1;
      fprintf(pf,"  %f  %f     %f       %f   %i     %12.10f  %12.10e\n",
                  VFa[TK],SFAC[F1],SFAC[F2],Temp,NLa[NI],SEA[i],SDA[i]);
    }


  }
  }
  }
  }

fclose(fp);
fclose(pf);
return 0;

}


// Functions
int *VEC_INT(int dim){
  int *v,i;
  v = (int *)malloc(dim*sizeof(int));
  if (v==NULL) {
     printf("\n\nVEC_INT: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0;
  return v;
}
double *VEC_DOUBLE(int dim){
  int i;
  double *v;
  v = (double *)malloc(dim*sizeof(double));
  if (v==NULL) {
     printf("\n\nVEC_DOUBLE: Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0;
  return v;
}

double complex *VEC_CDOUBLE(int dim) {
  int i;
  double complex *v;
  v = (double complex *)malloc(dim*sizeof(double complex));
  if (v==NULL) {
     printf("\n\nVEC_CDOUBLE:  Memory allocation error\n\n");
     exit(0);
  }
  for (i=0; i<dim; i++) v[i] = 0.0 + I*0.;
  return v;
}

char *VEC_CHAR(int dim){
  char *v;
  v = (char *)malloc(dim*sizeof(char));
  if (v==NULL) {
     printf("\n\nVEC_CHAR: Memory allocation error\n\n");
     exit(0);
  }
  return v;
}

void TransferMatrix(double thetaI, double k0, double complex *rind, double *d, 
double complex *cosL, double *beta, double *alpha, double complex *m11, double complex *m21) { 
  
  int i, j, k, indx;
  double complex *kz, *phiL, *D, *Dinv, *Pl, *phil;
  double complex *EM, ctheta, tmp, *tmp2, *tmp3, c0, c1, ci, kx;

  kz   = VEC_CDOUBLE(Nlayer);
  phil = VEC_CDOUBLE(Nlayer);
  D    = VEC_CDOUBLE(4*Nlayer);
  Dinv = VEC_CDOUBLE(4*Nlayer);
  Pl   = VEC_CDOUBLE(4*Nlayer);
  EM   = VEC_CDOUBLE(4);
  tmp2 = VEC_CDOUBLE(4); 
  tmp3 = VEC_CDOUBLE(4);

  c0 = 0. + I*0.;
  c1 = 1. + I*0.;
  ci = 0. + I*1.;

  //  x-component of incident wavevector...
  //  should be in dielectric material, so the imaginary 
  //  component should be 0.
  kx = k0*rind[0]*sin(thetaI);

  //  Now get the z-components of the wavevector in each layer
  for (i=0; i<Nlayer; i++) {
     kz[i] = (rind[i]*k0)*(rind[i]*k0) - kx*kx;
     kz[i] = csqrt(kz[i]);
     // Want to make sure the square root returns the positive imaginary branch
     if (cimag(kz[i])<0.)  {
        kz[i] = creal(kz[i]) - cimag(kz[i]);
     }
   }



   //  Calculate the P matrix
   for (i=1; i<Nlayer-1; i++) {
     phil[i]=kz[i]*d[i];

     //  Upper left (diagonal 1)
     Pl[i*4] = cexp(-ci*phil[i]);  
     //  upper right (off diagonal 1)
     Pl[i*4+1] = c0;
     //  lower left (off diagonal 2)
     Pl[i*4+2] = c0;
     //  lower right (diagonal 2)
     Pl[i*4+3] = cexp(ci*phil[i]);

   }

 
   //  Calculate the D and Dinv matrices
   for (i=0; i<Nlayer; i++) {
     ctheta = kz[i]/(rind[i]*k0);
     //  p-polarized incident waves
     if (polflag==1) {  

       //  Upper left (diagonal 1)
       D[i*4] = ctheta;
       // upper right
       D[i*4+1] = ctheta;
       // lower left
       D[i*4+2] = rind[i];
       // lower right
       D[i*4+3] = -rind[i];

     } 
     //  s-polarized incident waves
     if (polflag==2) {

       // upper left
       D[i*4] = 1;
       // upper right
       D[i*4+1] = 1;
       // lower left
       D[i*4+2] = rind[i]*ctheta;
       // lower right
       D[i*4+3] = -1*rind[i]*ctheta;

     }
     //  Now compute inverse
     //  Compute determinant of each D matrix
     tmp = D[i*4]*D[i*4+3]-D[i*4+1]*D[i*4+2];
     tmp = 1./tmp;

     //printf("  tmp is %12.10f  %12.10f\n",creal(tmp),cimag(tmp));    
     Dinv[i*4]=tmp*D[i*4+3];
     Dinv[i*4+1]=-1*tmp*D[i*4+1];
     Dinv[i*4+2]=-1*tmp*D[i*4+2];
     Dinv[i*4+3]=tmp*D[i*4];
 
   }


   // Initial EM matrix
   EM[0] = c1;
   EM[1] = c0;
   EM[2] = c0;
   EM[3] = c1;
   for (i=Nlayer-2; i>0; i--) {
     CMatMult2x2(i, Pl  , i,  Dinv, 0, tmp2);
     CMatMult2x2(i, D   , 0, tmp2,  0, tmp3);
     CMatMult2x2(0, tmp3, 0, EM  ,  0, tmp2); 

     for (j=0; j<2; j++) {
       for (k=0; k<2; k++) {
          EM[2*j+k] = tmp2[2*j+k];
       }
     }
   }
   CMatMult2x2(0, EM  , Nlayer-1, D   , 0, tmp2);
   CMatMult2x2(0, Dinv, 0, tmp2, 0, EM); 


   //  Finally, collect all the quantities we wish 
   //  to have available after this function is called
   *m11 = EM[0*2+0];  //  
   *m21 = EM[1*2+0];
   *beta = creal(kx);
   *alpha = cimag(kx);
   *cosL = ctheta;

   free(kz);  
   free(phil);
   free(D);
   free(Dinv);
   free(Pl);
   free(EM);
   free(tmp2);
   free(tmp3);

}

 

void CMatMult2x2(int Aidx, double complex *A, int Bidx, double complex *B, int Cidx, double complex *C) {
     int i, j, k, m, n;

     double complex sum;

     for (k=0; k<2; k++) {
       for (i=0; i<2; i++) {

          sum = 0. + 0*I;
          for (j=0; j<2; j++) {

            m = 2*i + j;
            n = 2*j + k;           
            sum += A[Aidx*4+m]*B[Bidx*4+n];

          }
          
          C[Cidx*4 + (2*i+k)] = sum;
        }
      }

}

double EvaluateMaterial(double *d, double complex *rind, double lambda_low, double lambda_hi, int numLambda, double lambda_bg) {

  int i, j;
  double k0, thetaI, beta, alpha, dLambda;
  double complex cosL, m11, m21, r, t;
  double R, T, A;

  double *lambda, *emissivity;

  lambda      = VEC_DOUBLE(numLambda);
  emissivity  = VEC_DOUBLE(numLambda);

  dLambda = (lambda_hi-lambda_low)/numLambda;
  thetaI = 0.;

  for (i=0; i<numLambda; i++) {

    lambda[i] = lambda_low+i*dLambda;
    // wavenumber
    k0 = 1000./lambda[i];
    // solve transfer matrix equations
    TransferMatrix(thetaI, k0, rind, d, &cosL, &beta, &alpha, &m11, &m21);
    // Definition of reflection amplitude
    r = m21/m11;
    // Fresnel transmission coefficient (also complex if there are absorbing layers)
    t = 1./m11;
    R = creal(r*conj(r));
    T = creal(rind[Nlayer]*cosL/(rind[0]*cos(thetaI))*t*conj(t));

    
  }
  free(lambda);
  free(emissivity);
}



void Bruggenman(double f, double epsD, double complex epsM, double *eta, double *kappa) {
  // medium 1 is surrounding medium (dielectric)
  // medium 2 is inclusion (W) - f passed to function is volume fraction of inclusion
  double f1, f2;
  double complex b, eps1, eps2, epsBG;
  eps1 = epsD + 0.*I;
  eps2 = epsM;


  f1 = (1 - f);
  f2 = f;
  b = (2*f1 - f2)*eps1 + (2*f2 - f1)*eps2;

  epsBG = (b + csqrt(8.*eps1*eps2 + b*b))/4.;

  // test to see that epsBG satisfy Bruggenman condition
  double complex test;
   *eta   = creal(csqrt(epsBG));
   *kappa = cimag(csqrt(epsBG));

}


void MaxwellGarnett(double f, double epsD, double complex epsM, double *eta, double *kappa) {
   double complex num, denom;

   num   = epsD*(2*f*(epsM - epsD) + epsM + 2*epsD);
   denom = 2*epsD + epsM + f*(epsD-epsM); 

   *eta   = creal(csqrt(num/denom));
   *kappa = cimag(csqrt(num/denom));

}

//  Evaluates real and imaginary part of refractive index from 
//  the Lorent oscillator model given omega_0, gamma_0, and omega
void Lorentz(double we, double de, double w, double *nreal, double *nimag) {

  double complex epsilon;
  double complex n;

  
  epsilon = 1 + pow(we,2)/(pow(we,2) - 2*I*de*w - pow(w,2));

  //printf("  w:  %12.10e  we:  %12.10f  de:  %12.10f  epsr:  %12.10f  epsi:  %12.10f\n",w,we,de,creal(epsilon),cimag(epsilon));
  n = csqrt(epsilon);

  *nreal = creal(n);
  *nimag = cimag(n);

}

int ReadDielectric(char *file, double *lambda, double complex *epsM) {
   int i;
   FILE *fp;
   double lam, epsr, epsi;

   fp = fopen(file,"r");

   i=0;
   while(!feof(fp)) {

     fscanf(fp, "%lf",&lam);
     fscanf(fp, "%lf",&epsr);
     fscanf(fp, "%lf",&epsi);

     lambda[i] = lam;
     epsM[i]   = epsr + I*epsi;

     i++;
   }

   printf("#  There are %i elements in file %s\n",i,file);
   fflush(stdout);
   return i;
   fclose(fp);
}

void ReadBRRind(int numBR, double lambda, double *BRlambda, double complex *BRind, double *n, double *k) {
  int i, fdx, bdx, die;
  double temp, eta, kappa;

  // The wavelength we are interested in is smaller than any in the range of data
  if (lambda<BRlambda[0]) {

    *n = creal(BRind[0]) + (lambda - BRlambda[0])*((creal(BRind[1]) - creal(BRind[0]))/(BRlambda[1] - BRlambda[0]));
    *k = cimag(BRind[0]) + (lambda - BRlambda[0])*((cimag(BRind[1]) - cimag(BRind[0]))/(BRlambda[1] - BRlambda[0]));


  }
  // The wavelength we are interested in is larger than any in the range of data
  else if (lambda>BRlambda[numBR-2]) {

    *n = creal(BRind[numBR-2]) +(lambda - BRlambda[numBR-2])*((creal(BRind[numBR-2]) - creal(BRind[numBR-3]))/(BRlambda[numBR-2] - BRlambda[numBR-3]));
    *k = cimag(BRind[numBR-2]) +(lambda - BRlambda[numBR-2])*((cimag(BRind[numBR-2]) - cimag(BRind[numBR-3]))/(BRlambda[numBR-2] - BRlambda[numBR-3]));


  }
  // We need to scan the data to find the BRlambda for two lambdas that straddle the lambda of interest
  else {

    i=0; 
    die=1;
    do {

      temp = BRlambda[i];
      if (temp>lambda) {
      
        die=0;
        fdx = i;
        bdx = i-1; 

      }
      else i++; 

    }while(die);

    *n = creal(BRind[bdx]) + (lambda - BRlambda[fdx])*((creal(BRind[fdx]) - creal(BRind[bdx]))/(BRlambda[fdx] - BRlambda[bdx]));
    *k = cimag(BRind[bdx]) + (lambda - BRlambda[fdx])*((cimag(BRind[fdx]) - cimag(BRind[bdx]))/(BRlambda[fdx] - BRlambda[bdx]));
  
  }

}


int IsDominated(int idx, int LENGTH, double *O1,double *O2) {
  int i, is, rval;
  double Val1, Val2;

  Val1 = O1[idx];
  Val2 = O2[idx];

  // start by assuming solution is NOT dominated
  rval = 0;
  for (i=0; i<LENGTH; i++)

      if (i!=idx) {

        // Trying to maximize the function, xi dominates xidx if 
        // fj(xi) >= fj(xidx) for all j and fj(xi) < fj(xidx) for at least one j
        if ((O1[i]>=Val1 && O2[i]>=Val2) && (O1[i]>Val1 || O2[i]>Val2)) {

          //printf("  x%i is dominated by x%i\n",idx,i);
          //printf("  f1(%i):  %12.10f  f1(%i): %12.10f  f2(%i): %12.10f  f2(%i):  %12.10f\n",
          //idx,Val1,i,O1[i],idx,Val2,i,O2[i]);
          printf("  terminating early!  i is %i out of %i\n",i,LENGTH);
          i=LENGTH;
          rval = 1;

        }
      }
  return rval;
}

//  Legendre-Gaussian quadrature functions
/******************************************************************************/

void legendre_compute_glr ( int n, double x[], double w[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   19 October 2009
 *
 *                     Author:
 *
 *                         Original C version by Nick Hale.
 *                             This C version by John Burkardt.
 *
 *                               Reference:
 *
 *                                   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                       A fast algorithm for the calculation of the roots of special functions, 
 *                                           SIAM Journal on Scientific Computing,
 *                                               Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                 Parameters:
 *
 *                                                     Input, int N, the order.
 *
 *                                                         Output, double X[N], the abscissas.
 *
 *                                                             Output, double W[N], the weights.
 *                                                             */
{
  int i;
  double p;
  double pp;
  double w_sum;
/*
 *   Get the value and derivative of the N-th Legendre polynomial at 0.
 *   */
  legendre_compute_glr0 ( n, &p, &pp );
/*
 *   Either zero is a root, or we have to call a function to find the first root.
 *   */  
  if ( n % 2 == 1 )
  {
    x[(n-1)/2] = p;
    w[(n-1)/2] = pp;
  }
  else
  {
    legendre_compute_glr2 ( p, n, &x[n/2], &w[n/2] );
  }
/*
 *   Get the complete set of roots and derivatives.
 *   */
  legendre_compute_glr1 ( n, x, w );
/*
 *   Compute the weights.
 *   */
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( 1.0 - x[i] ) / ( 1.0 + x[i] ) / w[i] / w[i];
  }
  w_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr0 ( int n, double *p, double *pp )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   19 October 2009
 *
 *                     Author:
 *
 *                         Original C version by Nick Hale.
 *                             This C version by John Burkardt.
 *
 *                               Reference:
 *
 *                                   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                       A fast algorithm for the calculation of the roots of special functions, 
 *                                           SIAM Journal on Scientific Computing,
 *                                               Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                 Parameters:
 *
 *                                                     Input, int N, the order of the Legendre polynomial.
 *
 *                                                         Output, double *P, *PP, the value of the N-th Legendre polynomial
 *                                                             and its derivative at 0.
 *                                                             */
{
  double dk;
  int k;
  double pm1;
  double pm2;
  double ppm1;
  double ppm2;

  pm2 = 0.0;
  pm1 = 1.0;
  ppm2 = 0.0;
  ppm1 = 0.0;

  for ( k = 0; k < n; k++ )
  {
    dk = ( double ) k;
    *p = - dk * pm2 / ( dk + 1.0 );
    *pp = ( ( 2.0 * dk + 1.0 ) * pm1 - dk * ppm2 ) / ( dk + 1.0 );
    pm2 = pm1;
    pm1 = *p;
    ppm2 = ppm1;
    ppm1 = *pp;
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr1 ( int n, double *x, double *ders )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
 *
 *         Discussion:
 *
 *             This routine requires that a starting estimate be provided for one
 *                 root and its derivative.  This information will be stored in entry
 *                     (N+1)/2 if N is odd, or N/2 if N is even, of ROOTS and DERS.
 *
 *                       Licensing:
 *
 *                           This code is distributed under the GNU LGPL license. 
 *
 *                             Modified:
 *
 *                                 19 October 2009
 *
 *                                   Author:
 *
 *                                       Original C version by Nick Hale.
 *                                           This C version by John Burkardt.
 *
 *                                             Reference:
 *
 *                                                 Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                                     A fast algorithm for the calculation of the roots of special functions, 
 *                                                         SIAM Journal on Scientific Computing,
 *                                                             Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                               Parameters:
 *
 *                                                                   Input, int N, the order of the Legendre polynomial.
 *
 *                                                                       Input/output, double X[N].  On input, a starting value
 *                                                                           has been set in one entry.  On output, the roots of the Legendre 
 *                                                                               polynomial.
 *
 *                                                                                   Input/output, double DERS[N].  On input, a starting value
 *                                                                                       has been set in one entry.  On output, the derivatives of the Legendre 
 *                                                                                           polynomial at the zeros.
 *
 *                                                                                             Local Parameters:
 *
 *                                                                                                 Local, int M, the number of terms in the Taylor expansion.
 *                                                                                                 */
{
  double dk;
  double dn;
  double h;
  int j;
  int k;
  int l;
  int m = 30;
  int n2;
  const double pi = 3.141592653589793;
  int s;
  double *u;
  double *up;
  double xp;

  if ( n % 2 == 1 )
  {
    n2 = ( n - 1 ) / 2;
    s = 1;
  }
  else
  {
    n2 = n / 2;
    s = 0;
  }

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;

  for ( j = n2; j < n - 1; j++ )
  {
    xp = x[j];

    h = rk2_leg ( pi/2.0, -pi/2.0, xp, n ) - xp;

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = ders[j];

    up[0] = 0.0;
    up[1] = u[2];

    for ( k = 0; k <= m - 2; k++ )
    {
      dk = ( double ) k;

      u[k+3] = 
      ( 
        2.0 * xp * ( dk + 1.0 ) * u[k+2]
        + ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1] / ( dk + 1.0 )
      ) / ( 1.0 - xp ) / ( 1.0 + xp ) / ( dk + 2.0 );

      up[k+2] = ( dk + 2.0 ) * u[k+3];
    }

    for ( l = 0; l < 5; l++ )
    { 
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 );
    }

    x[j+1] = xp + h;
    ders[j+1] = ts_mult ( up, h, m-1 );
  }

  free ( u );
  free ( up );

  for ( k = 0; k < n2 + s; k++ )
  {
    x[k] = - x[n-k-1];
    ders[k] = ders[n-k-1];
  }
  return;
}
/******************************************************************************/

void legendre_compute_glr2 ( double pn0, int n, double *x1,  double *d1 )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_COMPUTE_GLR2 finds the first real root.
 *
 *         Discussion:
 *
 *             This routine is only called if N is even.
 *
 *                 Thanks to Morten Welinder, for pointing out a typographical error
 *                     in indexing, 17 May 2013.
 *
 *                       Licensing:
 *
 *                           This code is distributed under the GNU LGPL license. 
 *
 *                             Modified:
 *
 *                                 17 May 2013
 *
 *                                   Author:
 *
 *                                       Original C version by Nick Hale.
 *                                           This C version by John Burkardt.
 *
 *                                             Reference:
 *
 *                                                 Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                                     A fast algorithm for the calculation of the roots of special functions, 
 *                                                         SIAM Journal on Scientific Computing,
 *                                                             Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                               Parameters:
 *
 *                                                                   Input, double PN0, the value of the N-th Legendre polynomial at 0.
 *
 *                                                                       Input, int N, the order of the Legendre polynomial.
 *
 *                                                                           Output, double *X1, the first real root.
 *
 *                                                                               Output, double *D1, the derivative at X1.
 *
 *                                                                                 Local Parameters:
 *
 *                                                                                     Local, int M, the number of terms in the Taylor expansion.
 *                                                                                     */
{
  double dk;
  double dn;
  int k;
  int l;
  int m = 30;
  const double pi = 3.141592653589793;
  double t;
  double *u;
  double *up;

  t = 0.0;
  *x1 = rk2_leg ( t, -pi/2.0, 0.0, n );

  u = ( double * ) malloc ( ( m + 2 ) * sizeof ( double ) );
  up = ( double * ) malloc ( ( m + 1 ) * sizeof ( double ) );

  dn = ( double ) n;
/*
 *   U[0] and UP[0] are never used.
 *     U[M+1] is set, but not used, and UP[M] is set and not used.
 *       What gives?
 *       */
  u[0] = 0.0;
  u[1] = pn0;

  up[0] = 0.0;
 
  for ( k = 0; k <= m - 2; k = k + 2 )
  {
    dk = ( double ) k;

    u[k+2] = 0.0;
    u[k+3] = ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1]
      / ( dk + 1.0 ) / ( dk + 2.0 );
 
    up[k+1] = 0.0;
    up[k+2] = ( dk + 2.0 ) * u[k+3];
  }
  
  for ( l = 0; l < 5; l++ )
  {
    *x1 = *x1 - ts_mult ( u, *x1, m ) / ts_mult ( up, *x1, m-1 );
  }
  *d1 = ts_mult ( up, *x1, m-1 );

  free ( u );
  free ( up) ;

  return;
}
/******************************************************************************/

void legendre_handle ( int n, double a, double b )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   22 October 2009
 *
 *                     Author:
 *
 *                         John Burkardt
 *
 *                           Parameters:
 *
 *                               Input, int N, the order of the rule.
 *
 *                                   Input, double A, B, the left and right endpoints.
 *                                   */ 
{
  int i;
  char output_r[255];
  char output_w[255];
  char output_x[255];
  double *r;
  double t;
  double *w;
  double *x;

  r = ( double * ) malloc ( 2 * sizeof ( double ) );
  w = ( double * ) malloc ( n * sizeof ( double ) );
  x = ( double * ) malloc ( n * sizeof ( double ) );

  r[0] = a;
  r[1] = b;
/*
 *   Compute the rule.
 *   */
  t = wtime ( );
  legendre_compute_glr ( n, x, w );
  t = wtime ( ) - t;

  printf ( "\n" );
  printf ( "  Elapsed time during computation was %g seconds.\n", t );
/*
 *   Rescale the rule to [A,B].
 *   */
  rescale ( a, b, n, x, w );
/*
 *   Write the rule to 3 files.
 *   */
  sprintf ( output_w, "leg_o%d_w.txt", n );
  sprintf ( output_x, "leg_o%d_x.txt", n );
  sprintf ( output_r, "leg_o%d_r.txt", n );

  printf ( "\n" );
  printf ( "  Weight file will be   \"%s\".\n", output_w );
  printf ( "  Abscissa file will be \"%s\".\n", output_x );
  printf ( "  Region file will be   \"%s\".\n", output_r );
            
  r8mat_write ( output_w, 1, n, w );
  r8mat_write ( output_x, 1, n, x );
  r8mat_write ( output_r, 1, 2, r );

  free ( r );
  free ( w );
  free ( x );

  return;
}
/******************************************************************************/

void r8mat_write ( char *output_filename, int m, int n, double table[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       R8MAT_WRITE writes an R8MAT file.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   01 June 2009
 *
 *                     Author:
 *
 *                         John Burkardt
 *
 *                           Parameters:
 *
 *                               Input, char *OUTPUT_FILENAME, the output filename.
 *
 *                                   Input, int M, the spatial dimension.
 *
 *                                       Input, int N, the number of points.
 *
 *                                           Input, double TABLE[M*N], the table data.
 *                                           */
{
  int i;
  int j;
  FILE *output;
/*
 *   Open the file.
 *   */
  output = fopen ( output_filename, "wt" );

  if ( !output )
  {
    printf ( "\n" );
    printf ( "R8MAT_WRITE - Fatal error!\n" );
    printf ( "  Could not open the output file.\n" );
    return;
  }
/*
 *   Write the data.
 *   */
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      fprintf ( output, "  %24.16g", table[i+j*m] );
    }
    fprintf ( output, "\n" );
  }
/*
 *   Close the file.
 *   */
  fclose ( output );

  return;
}
/******************************************************************************/

void rescale ( double a, double b, int n, double x[], double w[] )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   22 October 2009
 *
 *                     Author:
 *
 *                         Original MATLAB version by Nick Hale.
 *                             C version by John Burkardt.
 *
 *                               Reference:
 *
 *                                   Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
 *                                       A fast algorithm for the calculation of the roots of special functions, 
 *                                           SIAM Journal on Scientific Computing,
 *                                               Volume 29, Number 4, pages 1420-1438, 2007.
 *
 *                                                 Parameters:
 *
 *                                                     Input, double A, B, the endpoints of the new interval.
 *
 *                                                         Input, int N, the order.
 *
 *                                                             Input/output, double X[N], on input, the abscissas for [-1,+1].
 *                                                                 On output, the abscissas for [A,B].
 *
 *                                                                     Input/output, double W[N], on input, the weights for [-1,+1].
 *                                                                         On output, the weights for [A,B].
 *                                                                         */
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}
/******************************************************************************/

double rk2_leg ( double t1, double t2, double x, int n )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       RK2_LEG advances the value of X(T) using a Runge-Kutta method.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   22 October 2009
 *
 *                     Author:
 *
 *                         Original C version by Nick Hale.
 *                             This C version by John Burkardt.
 *
 *                               Parameters:
 *
 *                                   Input, double T1, T2, the range of the integration interval.
 *
 *                                       Input, double X, the value of X at T1.
 *
 *                                           Input, int N, the number of steps to take.
 *
 *                                               Output, double RK2_LEG, the value of X at T2.
 *                                               */
{
  double f;
  double h;
  int j;
  double k1;
  double k2;
  int m = 10;
  double snn1;
  double t;

  h = ( t2 - t1 ) / ( double ) m;
  snn1 = sqrt ( ( double ) ( n * ( n + 1 ) ) );

  t = t1;

  for ( j = 0; j < m; j++ )
  {
    f = ( 1.0 - x ) * ( 1.0 + x );
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + k1;

    t = t + h;

    f = ( 1.0 - x ) * ( 1.0 + x );
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );   
    x = x + 0.5 * ( k2 - k1 );
  }
  return x;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       TIMESTAMP prints the current YMDHMS date as a time stamp.
 *
 *         Example:
 *
 *             31 May 2001 09:45:54 AM
 *
 *               Licensing:
 *
 *                   This code is distributed under the GNU LGPL license. 
 *
 *                     Modified:
 *
 *                         24 September 2003
 *
 *                           Author:
 *
 *                               John Burkardt
 *
 *                                 Parameters:
 *
 *                                     None
 *                                     */
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
/******************************************************************************/

double ts_mult ( double *u, double h, int n )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       TS_MULT evaluates a polynomial.
 *
 *         Discussion:
 *
 *             TS_MULT = U[1] + U[2] * H + ... + U[N] * H^(N-1).
 *
 *               Licensing:
 *
 *                   This code is distributed under the GNU LGPL license. 
 *
 *                     Modified:
 *
 *                         17 May 2013
 *
 *                           Author:
 *
 *                               Original C version by Nick Hale.
 *                                   This C version by John Burkardt.
 *
 *                                     Parameters:
 *
 *                                         Input, double U[N+1], the polynomial coefficients.
 *                                             U[0] is ignored.
 *
 *                                                 Input, double H, the polynomial argument.
 *
 *                                                     Input, int N, the number of terms to compute.
 *
 *                                                         Output, double TS_MULT, the value of the polynomial.
 *                                                         */
{
  double hk;
  int k;
  double ts;
  
  ts = 0.0;
  hk = 1.0;
  for ( k = 1; k<= n; k++ )
  {
    ts = ts + u[k] * hk;
    hk = hk * h;
  }
  return ts;
}
/******************************************************************************/

double wtime ( void )

/******************************************************************************/
/*
 *   Purpose:
 *
 *       WTIME estimates the elapsed wall clock time.
 *
 *         Licensing:
 *
 *             This code is distributed under the GNU LGPL license. 
 *
 *               Modified:
 *
 *                   21 October 2009
 *
 *                     Author:
 *
 *                         John Burkardt
 *
 *                           Parameters:
 *
 *                               Output, double WTIME, the current elapsed wall clock time.
 *                               */
{
  double now;

  now = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC; 

  return now;
}

double g_pq(double p, double q, double x) {

  int d;
  d = (int)(fabs(p-q));
  double g;
  g = 0.;
  if (p == q && p == 0) {

    g = 1 - x;

  }
  else if ( p == q && p > 0 ) {

    g = (1 - x)*cos(p*pi*x)/2. - sin(p*pi*x)/(2*p*pi);

  }
  else if ( (d % 2)==0) {

    g = (q*sin(q*pi*x) - p*sin(p*pi*x))/((p*p-q*q)*pi);
  }
  else g = 0.;

  return g;
}


