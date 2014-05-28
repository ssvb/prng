/*
 * By Bob Jenkins, public domain
 *
 * With a 4-term state, results are w, x+stuff, y+stuff, z+stuff, w+stuff.
 * Make sure we've mixed the state well enough that 1-bit differences at 
 * w are pretty different by the time we report w+stuff.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include "prng.h"

int nrounds = 4;
int step_i  = 1;
int step_j  = 1;
int step_k  = 1;

#define BUCKETS (nbits * 4)
#define LOGLEN  16

/* count how many bits are set in a 64-bit integer, returns 0..64 */
static inline uint64_t count(uint64_t x)
{
  uint32_t c1 = x, c2 = x >> 32;
  c1 = (c1 & 0x55555555) + ((c1>>1 ) & 0x55555555);
  c1 = (c1 & 0x33333333) + ((c1>>2 ) & 0x33333333);
  c1 = (c1 & 0x0f0f0f0f) + ((c1>>4 ) & 0x0f0f0f0f);
  c1 = (c1 & 0x00ff00ff) + ((c1>>8 ) & 0x00ff00ff);
  c1 = (c1 & 0x0000ffff) + ((c1>>16) & 0x0000ffff);

  c2 = (c2 & 0x55555555) + ((c2>>1 ) & 0x55555555);
  c2 = (c2 & 0x33333333) + ((c2>>2 ) & 0x33333333);
  c2 = (c2 & 0x0f0f0f0f) + ((c2>>4 ) & 0x0f0f0f0f);
  c2 = (c2 & 0x00ff00ff) + ((c2>>8 ) & 0x00ff00ff);
  c2 = (c2 & 0x0000ffff) + ((c2>>16) & 0x0000ffff);

  return c1 + c2;
}

/* initialize the data collection array */
static void datainit(uint64_t *data, uint64_t *data2)
{
  int i;
  for (i=0; i<BUCKETS; ++i) {
    data[i] = 0;   /* look for poor XOR mixing */
    data2[i] = 0;  /* look for poor additive mixing */
  }
}

/* gather statistics on len overlapping subsequences of length 5 each */
static void gather(ranctx *x, uint64_t *data, uint64_t *data2, uint64_t length)
{
  uint64_t i, j, h;
  uint64_t k;
  ranctx y;

  for (i=0; i<BUCKETS; ++i) {
    for (k=0; k<length; ++k) {
      y = *x;
      if (i < nbits)
	y.a ^= ((uint64_t)1 << i);
      else if (i < nbits * 2)
	y.b ^= ((uint64_t)1 << (i - nbits));
      else if (i < nbits * 3)
	y.c ^= ((uint64_t)1 << (i - nbits * 2));
      else if (i < nbits * 4)
	y.d ^= ((uint64_t)1 << (i - nbits * 3));
      for (j = 0; j < nrounds; ++j) {
	h = (ranval(x) ^ ranval(&y)) & x->bitmask; /* look for poor mixing */
      }
      data[i] += count(h);
      h ^= (h<<1);     /* graycode to look for poor additive mixing */
      h &= x->bitmask;
      data2[i] += count(h);
    }
  }
}

static double calc_worst_helper(uint64_t *data, uint64_t *data2, uint64_t length)
{
  int i;
  double worst = data[0];
  for (i=1; i<BUCKETS; ++i) {
    if (worst > data[i]) {
      worst = data[i];
    }
    if (worst > data2[i]) {
      worst = data2[i];
    }
  }
  worst /= length;
  return worst;
}

double calc_worst(double cutoff)
{
  int i;
  uint64_t data[BUCKETS];
  uint64_t data2[BUCKETS];
  ranctx r;
  double worst;

  (void)raninit(&r, 0);
  datainit(data, data2);
  gather(&r, data, data2, (1<<6));
  for (i=6; i<LOGLEN; ++i) {
    double adjusted_cutoff = cutoff - ((i+1)==LOGLEN ? 0 : 0.1);
    gather(&r, data, data2, (1<<i));
    worst = calc_worst_helper(data, data2, (1<<(i+1)));
    if (worst < adjusted_cutoff)
        return 0;
  }
  return worst;
}

typedef struct prng_result { int i, j, k; double worst; } prng_result;

int get_results(prng_result *results, int maxresults, double cutoff)
{
  double worst;
  int i, j, k;
  int nresults = 0;

  for (i = 0; i < nbits; i += step_i) {
    for (j = 0; j < nbits; j += step_j) {
      jjj = j;
      iii = i;
      for (k = 0; k < nbits; k += step_k) {
        kkk = k;
        worst = calc_worst(cutoff);
        if (worst >= cutoff) {
          results[nresults].i = i;
          results[nresults].j = j;
          results[nresults].k = k;
          results[nresults].worst = worst;
          nresults++;
        }
        if (nresults == maxresults)
            return nresults;
      }
    }
  }
  return nresults;
}

#define MAX_RESULTS    25
#define WANTED_RESULTS 10

void driver(double initial_cutoff)
{
  /* try to find the best cutoff point */
  double maxcutoff = nbits / 2;
  double mincutoff = 0;
  double cutoff = (maxcutoff + mincutoff) / 2;
  int loopcount = 0;

  if (initial_cutoff > mincutoff &&
      initial_cutoff < maxcutoff) {
    cutoff = initial_cutoff;
  }

  int nresults, i;
  prng_result results[MAX_RESULTS];
  while (1) {
    loopcount++;
    printf("trying cutoff=%7.4f ...", cutoff);
    fflush(stdout);
    nresults = get_results(results, MAX_RESULTS, cutoff);
    printf(" got %d results\n", nresults);
    if (nresults >= MAX_RESULTS) {
        /* cutoff is too low */
        mincutoff = cutoff;
        if (loopcount >= 10)
          break;
    } else if (nresults < WANTED_RESULTS) {
        /* cutoff is too high */
        maxcutoff = cutoff;
        if (loopcount >= 11)
          break;
    } else {
      break;
    }
    cutoff = (maxcutoff + mincutoff) / 2;
  }
  printf("-\n");
  for (i = 0; i < nresults; i++) {
    double forward_avalanche, reverse_avalanche;
    iii = results[i].i;
    jjj = results[i].j;
    kkk = results[i].k;
    if (!reverse) {
      forward_avalanche = results[i].worst;
      reverse = !reverse;
      reverse_avalanche = calc_worst(0);
      reverse = !reverse;
    } else {
      reverse_avalanche = results[i].worst;
      reverse = !reverse;
      forward_avalanche = calc_worst(0);
      reverse = !reverse;
    }
    /* This output can be used by the cycle-statistics program */
    printf("  {.shifts={%2d, %2d, %2d}, .nbits=%2d, .fwd_avalanche=%7.4f, .rev_avalanche=%7.4f},\n",
           results[i].i, results[i].j, results[i].k,
           nbits, forward_avalanche, reverse_avalanche);
  }
}

int main(int argc, char *argv[])
{
  int i, j, k;
  time_t a,z;
  double initial_cutoff;

  if (argc < 4) {
    printf("Usage: rngav <nbits> <cutoff> [nrounds] [reverse] [step_i] [step_j] [step_k]\n");
    printf("Where: nbits     - size of each generated value in bits (up to 64)\n");
    printf("       cutoff    - minimum avalanche threshold\n");
    printf("       nrounds   - rounds before calculating avalanche (default 4)\n");
    printf("       reverse   - test the reverse generator if nonzero\n");
    exit(1);
  }

  nbits = atoi(argv[1]);
  initial_cutoff = atof(argv[2]);
  if (argc >= 4)
    nrounds = atoi(argv[3]);
  if (argc >= 5)
    reverse = atoi(argv[4]);
  if (argc >= 6)
    step_i = atoi(argv[5]);
  if (argc >= 7)
    step_j = atoi(argv[6]);
  if (argc >= 8)
    step_k = atoi(argv[7]);

  if (step_i < 1) step_i = nbits;
  if (step_j < 1) step_j = nbits;
  if (step_k < 1) step_k = nbits;

  nbits = (nbits < 4) ? 4 : nbits;
  nbits = (nbits > 64) ? 64 : nbits;

  printf("nbits   = %d\n", nbits);
  printf("cutoff  = %.2f\n", initial_cutoff);
  printf("nrounds = %d\n", nrounds);
  printf("reverse = %d\n", reverse);
  printf("step_i  = %d\n", step_i);
  printf("step_j  = %d\n", step_j);
  printf("step_k  = %d\n", step_k);
  printf("-\n");

  time(&a);

  driver(initial_cutoff);

  time(&z);

  printf("number of seconds: %6d\n", (int)(z-a));
}
