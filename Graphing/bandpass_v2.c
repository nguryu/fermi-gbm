#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>


double tukey(double alpha, int n, int N)
{
  int nmin = (int)floor(alpha*((double)N-1.0)/2.0);
  int nmax = (int)floor((double)(N-1)*(1. - alpha/2.));

  if (n>=0 && n<nmin)
  {
    return 0.5*( 1.0 + cos(M_PI*( (2.*(double)n)/(alpha*((double)N-1.0)) - 1.0 )) );
  }

  else if (n>=nmin && n<nmax)
  {
    return 1.0;
  }

  else if (n>=nmax && n<N)
  {
    return 0.5*( 1.0 + cos(M_PI*( (2.*(double)n)/(alpha*((double)N-1.0)) - 2.0/alpha + 1.0 )) );
  }

  else return 0.0;
}

int main(int argc, char* argv[])
{

  printf("\nProgram Name: %s\n",argv[0]);

  FILE *galFile = fopen("Galaxy_XAE.dat","r");
  FILE *resFile = fopen("Galaxy_XAE_R1.dat","r");
  FILE *conFile = fopen("Confusion_XAE_1.dat","r");

  int NGAL = 2516582;//1258291;
  int NCON = 2510293;//1255148;

  double TOBS = 62914560;//31457280;
  double minf, maxf;

  minf = atof(argv[1]);
  maxf = atof(argv[2]);

    printf("min: %f",minf);
    printf("\nmax: %f\n",maxf);


  int n = (int)floor(log2((double)NGAL));
  n++;

  int NFFT = (int)pow(2,n);

  printf("NFFT = %i = 2^%i\n",NFFT,n);

  double **resmatrix = malloc(NGAL*sizeof(double *));
  for(n=0; n<NGAL; n++) resmatrix[n] = malloc(3*sizeof(double));

  double **galmatrix = malloc(NGAL*sizeof(double *));
  for(n=0; n<NGAL; n++) galmatrix[n] = malloc(3*sizeof(double));

  double **conmatrix = malloc(NCON*sizeof(double *));
  for(n=0; n<NCON; n++) conmatrix[n] = malloc(2*sizeof(double));

  double junk;

  for(n=0; n<NGAL; n++)
  {
    fscanf(galFile,"%lg %lg %lg %lg %lg %lg %lg",&galmatrix[n][0],&junk,&junk,&galmatrix[n][1],&galmatrix[n][2],&junk,&junk);
  }


  for(n=0; n<NGAL; n++)
  {
    fscanf(resFile,"%lg %lg %lg %lg %lg %lg %lg",&resmatrix[n][0],&junk,&junk,&resmatrix[n][1],&resmatrix[n][2],&junk,&junk);
  }


  for(n=0; n<NCON; n++)
  {
    fscanf(conFile,"%lg %lg %lg %lg %lg",&conmatrix[n][0],&junk,&junk,&conmatrix[n][1],&junk);
  }


  //which bin does the confusion noise file start?
  int qmin = (int)floor(conmatrix[0][0]*TOBS);
  int qmax = (int)floor(conmatrix[NCON-1][0]*TOBS);


  FILE *whitefile = fopen("Galaxy_XAE_white.dat","w");
  FILE *whitefile2  = fopen("Galaxy_XAE_R1_white.dat","w");
  double f,re,im,re2,im2;

//  int nmin = (int)floor(.0004*TOBS);
//  int nmax = (int)floor(.004*TOBS);

  int nmin = (int)floor(minf*TOBS);
  int nmax = (int)floor(maxf*TOBS);

  int Nband = nmax-nmin;

  for(n=0; n<NGAL; n++)
  {

    f = (double)n/TOBS;

    if(n<qmin)
    {
      re = 0.0;
      im = 0.0;
      re2= 0.0;
      im2 = 0.0;
    }
    else if(n<qmax)
    {
      re  = galmatrix[n][1]/sqrt(conmatrix[n-qmin][1]);
      im  = galmatrix[n][2]/sqrt(conmatrix[n-qmin][1]);
      re2 = resmatrix[n][1]/sqrt(conmatrix[n-qmin][1]);
      im2 = resmatrix[n][2]/sqrt(conmatrix[n-qmin][1]);
    }
    else
    {
      re = 0.0;
      im = 0.0;
      re2=0.0;
      im2=0.0;
    }

    //bandpass
      double window = tukey(0.1,n-nmin,Nband);
      re *= window;
      im *= window;
      re2*= window;
      im2*= window;

    if(re!=re ||re2!=re2 ||im!=im ||im2!=im2)
    {
      printf("NAN!\n");
      printf("   galmatrix[n][1]=%lg\n",galmatrix[n][1]);
      printf("   galmatrix[n][2]=%lg\n",galmatrix[n][2]);
      printf("   resmatrix[n][1]=%lg\n",resmatrix[n][1]);
      printf("   resmatrix[n][2]=%lg\n",resmatrix[n][2]);
      printf("   conmatrix[%i][1]=%lg\n",n-qmin,conmatrix[n-qmin][1]);
    }

    fprintf(whitefile,"%.12g %.12g %.12g\n",f,re,im);
    fprintf(whitefile2,"%.12g %.12g %.12g\n",f,re2,im2);

  }
  for(n=NGAL; n<NFFT; n++)
  {

    f = (double)n/TOBS;
    re = 0.0;
    im = 0.0;
    fprintf(whitefile,"%.12g %.12g %.12g\n",f,re,im);
    fprintf(whitefile2,"%.12g %.12g %.12g\n",f,re2,im2);

  }
  fclose(whitefile);
  fclose(whitefile2);


	return 0;
}
