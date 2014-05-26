/*
 * Original code by Bob Jenkins, public domain
 * Modified by Siarhei Siamashka, public domain
 */

#include <assert.h>
#include "prng.h"

int iii = 23, jjj = 16, kkk = 11;
int nbits = 32;
int reverse = 0;

#define rot(x,k,nbits) ((x<<(k)) | (x>>(nbits-(k))))

uint64_t ranval(ranctx *x)
{
  uint64_t e;
  if (!x->reverse) {
    x->b &= x->bitmask;
    x->c &= x->bitmask;
    x->d &= x->bitmask;
    e = x->a - rot(x->b, x->iii, x->nbits);
    x->a = x->b ^ rot(x->c, x->jjj, x->nbits);
    x->b = x->c + rot(x->d, x->kkk, x->nbits);
    x->c = x->d + e;
    x->d = (e + x->a) & x->bitmask;
  } else {
    e = x->d - x->a;
    x->d = (x->c - e) & x->bitmask;
    x->c = (x->b - rot(x->d, x->kkk, x->nbits)) & x->bitmask;
    x->b = (x->a ^ rot(x->c, x->jjj, x->nbits)) & x->bitmask;
    x->a = e + rot(x->b, x->iii, x->nbits);
  }
  return x->d;
}

void raninit(ranctx *x, uint64_t seed)
{
  int i;
  x->nbits = nbits;
  x->bitmask = ((uint64_t)-1) >> (64 - nbits);
  x->reverse = reverse;
  x->iii = iii;
  x->jjj = jjj;
  x->kkk = kkk;
  x->a = x->b = x->c = 0xf1ea5eed & x->bitmask;
  x->d = (seed - x->a) & x->bitmask;

#ifndef NDEBUG
  x->reverse ^= 1;
  for (i = 0; i< 20; ++i)
    ranval(x);
  x->reverse ^= 1;
  for (i = 0; i< 20; ++i)
    ranval(x);
  assert((x->a & x->bitmask) == (0xf1ea5eed & x->bitmask));
  assert((x->b & x->bitmask) == (0xf1ea5eed & x->bitmask));
  assert((x->c & x->bitmask) == (0xf1ea5eed & x->bitmask));
  assert((x->d & x->bitmask) == ((seed - 0xf1ea5eed) & x->bitmask));
#endif

  for (i = 0; i< 20; ++i)
    ranval(x);
}
