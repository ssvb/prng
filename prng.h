/*
 * Original code by Bob Jenkins, public domain
 * Modified by Siarhei Siamashka, public domain
 */

#include <stdint.h>

extern int iii, jjj, kkk;
extern int nbits;
extern int reverse;

typedef struct {
  uint64_t a, b, c, d; /* prng state */
  int iii, jjj, kkk;   /* configurable rotate constants */
  int nbits;           /* number of bits */
  int reverse;         /* reverse the direction of prng numbers generation */
  uint64_t bitmask;    /* (1 << nbits) - 1 */
} ranctx;

void raninit(ranctx *x, uint64_t seed);
uint64_t ranval(ranctx *x);
