/*
 * By Siarhei Siamashka, public domain
 *
 * This program tries to analyze the minimum cycle length
 * when playing with the "toy" PRNG settings. Such settings
 * only use something like 8 bits per value/term (32 bits
 * per state).
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <assert.h>
#include "prng.h"

typedef unsigned long long u8;

static uint32_t count_bits(uint32_t x)
{
  uint32_t c = x;
  c = (c & 0x55555555) + ((c>>1 ) & 0x55555555);
  c = (c & 0x33333333) + ((c>>2 ) & 0x33333333);
  c = (c & 0x0f0f0f0f) + ((c>>4 ) & 0x0f0f0f0f);
  c = (c & 0x00ff00ff) + ((c>>8 ) & 0x00ff00ff);
  c = (c & 0x0000ffff) + ((c>>16) & 0x0000ffff);
  return c;
}

static inline uint64_t get_stateidx(ranctx *r)
{
    return ((r->a & r->bitmask) << (r->nbits * 0)) |
           ((r->b & r->bitmask) << (r->nbits * 1)) |
           ((r->c & r->bitmask) << (r->nbits * 2)) |
           ((r->d & r->bitmask) << (r->nbits * 3));
}

void driver(int i, int j, int k)
{
  uint32_t *buffer;
  uint64_t bufsize = ((uint64_t)1 << (nbits * 4)) / 8;
  uint64_t maxseed = (uint64_t)1 << nbits;
  uint64_t seed;
  uint64_t idx;
  uint64_t bit_counter = 0;
  u8 mincyclesize = (u8)-1;
  buffer = malloc(bufsize);
  memset(buffer, 0, bufsize);

  iii = i;
  jjj = j;
  kkk = k;

  for (seed = 0; seed < maxseed; seed++) {
    ranctx r;
    u8 cyclesize = 0;
    raninit(&r, seed);
    uint64_t stateidx = get_stateidx(&r);
    /* check if we have already seen this sequence */
    if (buffer[stateidx / 32] & (1 << (stateidx % 32)))
       continue;
    while (1) {
      stateidx = get_stateidx(&r);
      if (buffer[stateidx / 32] & (1 << (stateidx % 32)))
        break;
      buffer[stateidx / 32] |= 1 << (stateidx % 32);
      cyclesize++;
      ranval(&r);
    }
    if (cyclesize < mincyclesize)
      mincyclesize = cyclesize;
  }

  /* Calculate the percentage of reachable states (using 'nbit' sized seed) */
  for (idx = 0; idx < bufsize / 4; idx++) {
    bit_counter += count_bits(buffer[idx]);
  }

  printf("%2d-bit - {%2d, %2d, %2d} - smallest reachable cycle: %10llu, reachable states: %.4f%%\n",
         nbits, i, j, k, mincyclesize, (double)bit_counter / (bufsize * 8) * 100.);
  free(buffer);
}

typedef struct prng_info {
  struct {int i, j, k;} shifts;
  int nbits;
  double fwd_avalanche, rev_avalanche;
} prng_info;

prng_info tbl[] = {
  /* 6-bit variants */
  {.shifts={ 0,  2,  0}, .nbits= 6, .fwd_avalanche= 2.6770, .rev_avalanche= 2.8257},
  {.shifts={ 0,  3,  0}, .nbits= 6, .fwd_avalanche= 2.8265, .rev_avalanche= 3.0240},
  {.shifts={ 0,  4,  0}, .nbits= 6, .fwd_avalanche= 2.3596, .rev_avalanche= 2.4459},
  {.shifts={ 1,  1,  0}, .nbits= 6, .fwd_avalanche= 2.3739, .rev_avalanche= 2.8436},
  {.shifts={ 1,  2,  0}, .nbits= 6, .fwd_avalanche= 2.9189, .rev_avalanche= 2.9612},
  {.shifts={ 1,  3,  0}, .nbits= 6, .fwd_avalanche= 2.6622, .rev_avalanche= 2.9243},
  {.shifts={ 1,  4,  0}, .nbits= 6, .fwd_avalanche= 2.4111, .rev_avalanche= 2.9471},
  {.shifts={ 1,  5,  0}, .nbits= 6, .fwd_avalanche= 2.4036, .rev_avalanche= 1.9119},
  {.shifts={ 2,  1,  0}, .nbits= 6, .fwd_avalanche= 2.7524, .rev_avalanche= 2.9200},
  {.shifts={ 2,  2,  0}, .nbits= 6, .fwd_avalanche= 2.5965, .rev_avalanche= 2.9471},
  {.shifts={ 2,  3,  0}, .nbits= 6, .fwd_avalanche= 2.8779, .rev_avalanche= 2.9749},
  {.shifts={ 2,  4,  0}, .nbits= 6, .fwd_avalanche= 2.8296, .rev_avalanche= 2.9046},
  {.shifts={ 2,  5,  0}, .nbits= 6, .fwd_avalanche= 2.7168, .rev_avalanche= 2.9540},
  {.shifts={ 3,  1,  0}, .nbits= 6, .fwd_avalanche= 2.6709, .rev_avalanche= 2.9663},
  {.shifts={ 3,  2,  0}, .nbits= 6, .fwd_avalanche= 2.4210, .rev_avalanche= 2.9752},
  {.shifts={ 3,  3,  0}, .nbits= 6, .fwd_avalanche= 2.6925, .rev_avalanche= 2.4798},
  {.shifts={ 3,  4,  0}, .nbits= 6, .fwd_avalanche= 2.7191, .rev_avalanche= 2.9283},
  {.shifts={ 3,  5,  0}, .nbits= 6, .fwd_avalanche= 2.4842, .rev_avalanche= 2.9823},
  {.shifts={ 4,  1,  0}, .nbits= 6, .fwd_avalanche= 2.9082, .rev_avalanche= 2.9604},
  {.shifts={ 4,  2,  0}, .nbits= 6, .fwd_avalanche= 2.9215, .rev_avalanche= 2.8949},
  {.shifts={ 4,  3,  0}, .nbits= 6, .fwd_avalanche= 2.5426, .rev_avalanche= 2.9590},
  {.shifts={ 4,  5,  0}, .nbits= 6, .fwd_avalanche= 2.5329, .rev_avalanche= 2.9767},
  {.shifts={ 5,  2,  0}, .nbits= 6, .fwd_avalanche= 2.7266, .rev_avalanche= 2.7948},
  {.shifts={ 5,  3,  0}, .nbits= 6, .fwd_avalanche= 2.9112, .rev_avalanche= 2.9629},
};

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(a) (sizeof((a)) / sizeof((a)[0]))
#endif

int main(int argc, char *argv[])
{
  time_t a,z;
  int i;

  time(&a);
  for (i = 0; i < ARRAY_SIZE(tbl); i++) {
    nbits = tbl[i].nbits;
    driver(tbl[i].shifts.i, tbl[i].shifts.j, tbl[i].shifts.k);
  }
  time(&z);
  printf("number of seconds: %6d\n", (int)(z-a));
}
