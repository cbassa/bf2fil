#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "bf2fil.h"
#include "skz.h"

int read_hdf5_header(char *filename,header *h, int verbose);
int write_fil_header(char *filename,header *h);

int main(int argc,char *argv[])
{
  header h;
  int nx,ny,my,m,iblock;
  int64_t i,j,k,l,nread;
  float nd=12.0,sigma=4.0;
  double sk_lim[2];
  float zsig[]={3.0,5.0};
  int *mask,nmask,*nzc;
  clock_t startclock;
  float *z,*zsm,*zss,*zcm,*zcs,*zoffset,*zscale,ztmp;
  unsigned char *cz;
  FILE *infile,*outfile;
  int status;

  // Read HDF5 header
  read_hdf5_header(argv[1],&h,0);

  // Half channel offset
  h.fch1+=0.5*h.foff;

  // Write FIL header
  write_fil_header("test.fil",&h);

  // Setup sizes
  nx=h.nchan;
  m=1024;
  ny=120*m;
  my=(int) ceil(ny/(float) m);

  // Mask
  mask=(int *) malloc(sizeof(int)*nx*my);

  // Spectrum mean and standard deviation
  zsm=(float *) malloc(sizeof(int)*my);
  zss=(float *) malloc(sizeof(int)*my);

  // Channel mean and standard deviation
  zcm=(float *) malloc(sizeof(int)*nx);
  zcs=(float *) malloc(sizeof(int)*nx);
  nzc=(int *) malloc(sizeof(int)*nx);
  zoffset=(float *) malloc(sizeof(int)*nx);
  zscale=(float *) malloc(sizeof(int)*nx);

  // Allocate
  z=(float *) malloc(sizeof(float)*nx*ny);
  cz=(unsigned char *) malloc(sizeof(unsigned char)*nx*ny);

  // Compute spectral kurtosis limits
  status=sk_threshold3(m,(double) sigma,(double) nd,sk_lim);
  printf("Block size: %d MB, averaged spectra: %d, d: %g, sigma: %.1f\nSK limits: [%f,%f]\n",nx*ny*sizeof(float)/(1<<20),m,nd,sigma,sk_lim[0],sk_lim[1]);

  // Open RAW file
  infile=fopen(argv[2],"r");

  // Open FIL file
  outfile=fopen("test.fil","a");

  // Loop over file
  for (iblock=0;;iblock++) {
    // Read block
    startclock=clock();
    nread=fread(z,sizeof(float),nx*ny,infile);
    printf("Read block %d, %d MB in %.2f s\n",iblock,nread*sizeof(float)/(1<<20),(float) (clock()-startclock)/CLOCKS_PER_SEC);

    // Compute SK mask
    startclock=clock();
    nmask=compute_mask(z,nx,ny,my,m,(float) nd,sigma,sk_lim[0],sk_lim[1],mask,zcm,zcs,zsm,zss,nzc);
    printf("Computed mask (%d values masked) in %.2f s\n",nmask,(float) (clock()-startclock)/CLOCKS_PER_SEC);

    // Compute channel offsets and scales
    for (i=0;i<nx;i++) {
      zoffset[i]=zcm[i]-zsig[0]*zcs[i];
      zscale[i]=(zsig[0]+zsig[1])*zcs[i];
    }

    // Redigitize
    startclock=clock();
    for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
	k=i+nx*j;
        z[k]-=zoffset[i];
	z[k]*=256.0/zscale[i];
      }
    }
    printf("Rescaled in %.2f s\n",(float) (clock()-startclock)/CLOCKS_PER_SEC);
    startclock=clock();
    for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
	k=i+nx*j;
	l=(nx-i-1)+nx*j;
	ztmp=z[k];
	cz[l]=(unsigned char) ztmp;
	if (ztmp<0.0) cz[l]=0;
	if (ztmp>255.0) cz[l]=255;
      }
    }
    printf("Converted in %.2f s\n",(float) (clock()-startclock)/CLOCKS_PER_SEC);

    // Write
    startclock=clock();
    fwrite(cz,sizeof(unsigned char),nread,outfile);
    printf("Wrote block %d, %d MB in %.2f s\n",iblock,nread*sizeof(float)/(1<<20),(float) (clock()-startclock)/CLOCKS_PER_SEC);

    // Exit on last block
    if (nread<nx*ny)
      break;
  }


  // Close files
  fclose(infile);
  fclose(outfile);

  // Free
  free(z);
  free(cz);
  free(mask);
  free(zsm);
  free(zss);
  free(zcm);
  free(zcs);
  free(nzc);
  free(zoffset);
  free(zscale);

  return 0;
}
