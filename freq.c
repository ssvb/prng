/*
 * Measure whether all values appear with equal frequency
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

typedef  unsigned char      u1;
typedef  unsigned long      u4;
typedef  unsigned long long u8;

#define BUCKETS (1<<8)

typedef struct ranctx { u4 a; u4 b; u4 c; u4 d;} ranctx;

#define rot(x,k) ((x<<k)|(x>>(32-k)))

static u4 iii = 0;

static u4 ranval( ranctx *x ) {
  u4 e = x->a - rot(x->b, 27);
  x->a = x->b ^ rot(x->c, 17);
  x->b = x->c + x->d;
  x->c = x->d + e;
  x->d = e + x->a;
  return x->d; 
}

static void raninit( ranctx *x, u4 seed ) {
  u4 i;
  x->a = 0xf1ea5eed; 
  x->b = x->c = x->d = seed;
  for (i=0; i<20; ++i) {
    (void)ranval(x);
  }
}

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
  u4 loglen = 0;
  ranctx r;
  time_t a,z;
  
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
  printf("number of seconds: %6d\n", (size_t)(z-a));
}
