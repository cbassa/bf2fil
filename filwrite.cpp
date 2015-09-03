#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bf2fil.h"
#include "filterbank.h"

int write_fil_header(char *filename,header *h)
{
  FILE *file;
  int nsamples;

  // Copy parameters
  strcpy(source_name,h->source_name);
  machine_id=h->machine_id;
  telescope_id=h->telescope_id;
  data_type=1;
  nchans=h->nchan;
  nsamples=h->nsamp;
  nbits=8;
  obits=8;
  nifs=1;
  scan_number=0;
  barycentric=0,pulsarcentric=0;
  headerless=0;
  tstart=h->tstart;
  mjdobs=0;
  tsamp=h->tsamp;
  fch1=h->fch1;
  foff=h->foff;
  refdm=0.0,az_start=0.0,za_start=0.0;
  src_raj=h->src_raj;
  src_dej=h->src_dej;
  sumifs=1;

  file=fopen("test.fil","w");
  filterbank_header(file);
  fclose(file);

  return 0;
}
