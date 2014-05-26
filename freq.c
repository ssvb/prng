/*
 * Measure whether all values appear with equal frequency
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "prng.h"

typedef  unsigned char      u1;
typedef  unsigned int       u4;
typedef  unsigned long long u8;

#define BUCKETS (1<<8)

/* initialize the data collection array */
static void datainit( u8 *data)
{
  u4 i;
  for (i=0; i<BUCKETS; ++i) data[i] = (u8)0;
}

/* gather statistics on len overlapping subsequences of length 5 each */
#define INCREMENT 0x10000
static void gather( ranctx *r, u8 *data, u8 len)
{
  u8 i;
  u4 j;
  for (i=0; i<len; i+=INCREMENT) {
    for (j=0; j<INCREMENT; j+=4) {
      ++data[ranval(r)&(BUCKETS-1)];
      ++data[ranval(r)&(BUCKETS-1)];
      ++data[ranval(r)&(BUCKETS-1)];
      ++data[ranval(r)&(BUCKETS-1)];
    }
  }
}

static void chi( u8 *data, u8 len)
{
  u4 i;
  double var = 0.0;         /* total variance */
  double temp;              /* used to calculate variance of a bucket */
  double expect = ((double)len)/BUCKETS;
  
  for (i=0; i<BUCKETS; ++i) {
    double temp = (double)data[i] - expect;
    var += temp*temp/expect;
  }

  /* calculate the total variance and chi-square measure */
  printf("expected variance: %11.4f   got: %11.4f   chi-square: %6.4f\n",
         (float)(BUCKETS-1), (float)var, 
	 (float)((var-(BUCKETS-1))/sqrt((float)(BUCKETS-1))));
}


int main( int argc, char **argv)
{
  u8 len;
  u8 data[BUCKETS];
  int loglen = 0;
  ranctx r;
  time_t a,z;

  assert(sizeof(u1) == 1);
  assert(sizeof(u4) == 4);
  assert(sizeof(u8) == 8);

  time(&a);
  if (argc == 2) {
    sscanf(argv[1], "%d", &loglen);
    printf("log_2 sequence length: %d\n", loglen);
    len = (((u8)1)<<loglen);
  } else {
    fprintf(stderr, "usage: \"cou 24\" means run for 2^^24 values\n");
    return 1;
  }

  datainit(data);
  raninit(&r, 0);
  gather(&r, data, len);
  chi(data, len);

  time(&z);
  printf("number of seconds: %6d\n", (int)(z-a));
}
