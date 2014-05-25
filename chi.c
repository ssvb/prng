#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
/*
----------------------------------------------------------------------
By Bob Jenkins, amateur generator of random number generators, 1994
If you find a bug, you'll have to fix it yourself.  Public Domain.

Run chi-square tests on a random number generator rng.
  uni -- are values uniform; do they occur equally often
  gap -- are the gaps between values of the expected length
  run -- how long are strictly increasing subsequences?
  njk -- uni, but by conglomerating NJKBIT of multiple values.
  njg -- gap, but by conglomerating NJKBIT of multiple values.

Instructions: If "get" is within "expect +- 3*sqrt(expect)", then the
  tests pass.  Reduce ALPHA, OMEGA until the tests fail.  Increase
  MYRUNS until the tests fail.  Where does the bias come from?
----------------------------------------------------------------------
*/

typedef  unsigned      char u1;   /* u1 is unsigned, 1 byte  */
typedef  unsigned long int  u4;   /* u4 is unsigned, 4 bytes */
typedef  unsigned long long u8;   /* u8 is unsigned, 8 bytes */
static FILE *xx;

#define ALPHA 2                   /* arrays of size 2^alpha */
#define OMEGA 32                    /* omega bits per value */
#define OMICRON 4        /* number of bits in values tested */
#define RSIZE  ((u4)8)    /* maximum run length to look for */
#define NJKBIT 0x1 /* which bit to conglomerate for njk,njg */
#define SHIFT 4                    /* amount of barrelshift */
#define SIZE  ((u4)1<<ALPHA)
#define WSIZE ((u4)1<<OMEGA)
#define USIZE ((u4)1<<OMICRON)
#define GSIZE (SIZE<<3)
#define WMASK /* (WSIZE-1) */ 0xffffffff
#define UMASK (USIZE-1)
#define MYRUNS (((u8)1)<<28)
#define DUMP  "chi.txt"

static u4 cyca, cycb, cycz;
static u4 mem[SIZE];              /* m: pool of secret memory */
static u4 ss[SIZE];
static u4 tt[SIZE];
static u4 rsl[SIZE];              /* r: results given to the user */
static u4 cycr[SIZE];              /* r: results given to the user */
static u4 cycm[SIZE];
static u4 iii,jjj,kkk,lll;

static u8     unidata[USIZE];
static u8     njkdata[USIZE];
static double uniprob[USIZE];

static u8     gapdata[GSIZE];
static u4     gaplast[USIZE];
static u4     gapoff;
static u8     gap2data[GSIZE];
static u4     gap2last[USIZE];
static u4     gap2off;
static double gapprob[GSIZE];

static u8     rundata[RSIZE];
static u8     run2data[RSIZE];
static double runprob[RSIZE];
static u4     rcount, rlast;
static u4     rcount2, rlast2;

/* 
---------------------------------------------------------------------
PRNG() is a reliable outside source of random numbers
---------------------------------------------------------------------
*/

#define NNL (8)
#define NNN ((u4)1<<NNL)
static u4 mem2[NNN];              /* m: pool of secret memory */
static u4 rsl2[NNN];              /* r: results given to the user */
static u4 prngi=(NNN-1);                /* counter */

#define rotate(a)  (((a)<<19)^((a)>>(32-19)))
#define f(mm,x)    (*(u4 *)(((u1 *)(mm))+((x)&(4*NNN-4))))

static void prng_loop(mm,rr,aa,bb)
u4 *mm;   /* secret memory             */
u4 *rr;   /* results given to the user */
u4 *aa;  /* accumulator */
u4 *bb;  /* the old m[SIZE-1] */
{
  register u4 a,b,i,x,y;
 
  a = *aa;  b = *bb;
  for (i=0; i<NNN; ++i)
  {
    x = mm[i];
    a = rotate(a) + mm[((i+(NNN/2))&(NNN-1))];
    mm[i] = y = f(mm,x) + a + b;
    rr[i] = b = f(mm,(y>>NNL)) + x;
  }
  *bb = b;  *aa = a;
}
 

#define  prng(i,m,r,aa,bb) \
  if (1) { \
    if (++prngi>=NNN) {prng_loop((m),(r),(aa),(bb)); prngi=0;} \
    (i) = (r)[prngi]; \
  }


/*
------------------------------------------------------------------
Experimental RNG -- place your RNG here
Each call must fill rr[0..SIZE-1] with ALPHA-bit values.
------------------------------------------------------------------
*/
#define rot(d, k) (((d<<k)&WMASK)|((d>>(OMEGA-k))&WMASK))

static void rng(mm,rr,aa,bb,cc,m2,r2,aa2,bb2,ss,tt)
u4 *mm;   /* secret memory             */
u4 *rr;   /* results given to the user */
u4 *aa;  /* accumulator */
u4 *bb;  /* the old m[SIZE-1] */
u4 *cc;  /* the old m[SIZE-2] */
u4 *m2;   /* secret memory             */
u4 *r2;   /* results given to the user */
u4 *aa2;  /* accumulator */
u4 *bb2;  /* the old m[SIZE-1] */
u4 *ss;   /* extra state */
u4 *tt;   /* more extra state */
{
  register u4 a,b,c,d,x,y,z,i,j;

  a = *aa; b = *bb; c = *cc; d = mm[0];

  for (i=0; i<SIZE; ++i)
  {
    /* xxx */
    x = a;
    a = b;
    b = (c<<19) + (c>>13) + d;
    c = d ^ a;
    d = x + b;

    rr[i] = (c&0x3)|(((c>>19)&0x3)<<2);
  }

  *aa = a; *bb = b; *cc = c; mm[0] = d;
}

/*
------------------------------------------------------------------
cyc() assumes that the rng is a permutation, that is, the initial
state is guaranteed to be repeated.
cyc() relies on m, not r, so rng() needs to fill in m.
------------------------------------------------------------------
*/
static int  cyc( m, r, a, b, c, m2, r2, a2, b2, ss, tt, i8 )
u4 *m;  /* random numbers */
u4 *r;  /* random numbers */
u4  a;
u4  b;
u4  c;
u4 *m2;
u4 *r2;
u4 *a2;
u4 *b2;
u4 *ss;
u4 *tt;
u8  i8;
{
  int i;
#ifdef NEVER
  if ((i8&1) == 1) rng(cycm,cycr,&cyca,&cycb,&cycz,m2,r2,a2,b2,ss,tt);
#endif
  if (a != cyca) return 0;
  if (b != cycb) return 0;
  if (c != cycz) return 0;
  for (i=0; i<SIZE; ++i) if (m[i] != cycm[i]) return 0;
  return 1;
}

/*
------------------------------------------------------------------
This gathers statistics for njk and njg.  
Proposed by Niels Jo/rgen Kruse, hence the name.
Correlated values have OMICRON bits, just like normal values,
  except there is only one value per call to rng().  
If SIZE<OMICRON, then njk and njg results will be garbage.
------------------------------------------------------------------
*/
static void njk( r, rl, d, dl, d2, d2l)
u4 *r;  /* random numbers */
u4  rl; /* length of r */
u8 *d;  /* where to put measurements */
u4  dl; /* length of d */
u8 *d2;  /* where to put measurements */
u4  d2l; /* length of d */
{
  u4  i,j,x=0;
  for (i=0; i<OMICRON; ++i) x = (x<<1)|(!(r[i]&NJKBIT));
  ++d[x];
  j = (gap2off+i-1)-gap2last[x];
  if (j<d2l) ++d2[j];
  else ++d2[d2l-1];
  gap2last[x] = i+gap2off;
  gap2off = gap2off + 1;
}

/* How many times does each value appear */
static void uniform( r, rl, d, dl )
u4 *r;  /* random numbers */
u4  rl; /* length of r */
u8 *d;  /* where to put measurements */
u4  dl; /* length of d */
{
  u4  i;
  for (i=0; i<rl; ++i) ++d[r[i]];
}

static void uprob( prob, len)
double *prob;
u4      len;
{
  u4     i;
  double p=1/(double)len;

  for (i=0; i<len; ++i) prob[i]=p;
}


/* For each value, measure gaps between occurances of that value */
static void gap( r, rl, d, dl )
u4 *r;  /* random numbers */
u4  rl; /* length of r */
u8 *d;  /* where to put measurements */
u4  dl; /* length of d */
{
  register u4 i,j,x;
  
  for (i=0; i<rl; ++i)
  {
    if (r[i] < USIZE)
    {
      x = r[i];
      j = (gapoff+i-1)-gaplast[x];
      if (j<dl) ++d[j];
      else ++d[dl-1];
      gaplast[x] = i+gapoff;
    }
  }
  gapoff = gapoff + SIZE;
}

static void gprob( prob, len)
double *prob;
u4      len;
{
  register  u4     i,j,k;
  register  double P, R;
  double           z[16];

  P = 1.0/((double)USIZE);
  z[0] = (1-P);
  for (i=1; i<16; ++i) z[i] = z[i-1]*z[i-1];
  for (i=0; i<len; ++i)
  {
    for (R=P, k=0, j=i; j; j>>=1, ++k) if (j&1) R *= z[k];
    prob[i] = R;
  }
  prob[len-1] *= (1/P);
}



/* Count runs of strictly decreasing sequences */
static void run( r, rl, d, dl )
u4 *r;  /* random numbers */
u4  rl; /* length of r */
u8 *d;  /* where to put measurements */
u4  dl; /* length of d */
{
  u4 i,c,x;
  
  x = rlast; c = rcount;
  for (i=0; i<rl; ++i)
  {
    if (c > 1000) {c=0; x=r[i];}
    else if (x == r[i]) ;
    else if (x > r[i]) {++c; x=r[i];}
    else {if (c > dl-1) c=dl-1;   ++d[c]; c=1001;}
  }
  rcount = c; rlast = x;
}

/* Count runs of strictly increasing sequences */
static void run2( r, rl, d, dl )
u4 *r;  /* random numbers */
u4  rl; /* length of r */
u8 *d;  /* where to put measurements */
u4  dl; /* length of d */
{
  u4 i,c,x;
  
  x = rlast2; c = rcount2;
  for (i=0; i<rl; ++i)
  {
    if (c > 1000) {c=0; x=r[i];}
    else if (x == r[i]) ;
    else if (x < r[i]) {++c; x=r[i];}
    else {if (c > dl-1) c=dl-1;   ++d[c]; c=1001;}
  }
  rcount2 = c; rlast2 = x;
}


static void rprob( prob, len)
double *prob;
u4      len;
{
   u4  i,j,k;

   prob[0] = 1;
   for (i=1; i<=RSIZE; ++i)
   prob[i] = prob[i-1]*((double)USIZE-i)/((i+1)*((double)USIZE-1));
   for (i=0; i<RSIZE-1; ++i)
   prob[i] -= prob[i+1];
}



/* If data[] is random with distribution prob[], the result should be
   len +- 3sqrt(len) unless len is really small */
void chisquare( name, data, prob, len)
char   *name;
u8     *data;
double *prob;
u4      len;
{
  register u4     i,j,k;
  register double V=0, S=0, Q, R;

  for (i=0; i<len; ++i) S += data[i];
  for (i=0; i<len; ++i)
  {
    R = S*prob[i];
    Q = ((double)data[i] - R);
    V += (Q*Q)/R;
  }

  printf("%s: expect %5ld  get %12.4f  chi %12.4f\n",
         name, len-1, V, (V-(len-1))/sqrt((float)(len-1)));
  for (i=0; i<len; ++i) { 
    if (!(i&7)) printf("\n");
    R = S*prob[i];
    Q = ((double)data[i] - R);
    printf("%9.4f", (float)(Q*Q/R));
  }
  printf("\n");

}

void mydump(u8 i8, u4 *mem, u4 *rsl, u4 a, u4 b, u4 z, 
	    u4 *mem2, u4 *rsl2, u4 a2, u4 b2, u4 *ss, u4 *tt, 
	    u8 *unidata, u8 *gapdata, u8 *rundata, u8 *run2data, 
	    u8 *njkdata, u8 *gap2data, 
	    u4 *cycm, u4 cyca, u4 cycb, u4 cycz, 
	    u4 *gaplast, u4 *gap2last)
{
  u4 i;
  FILE *f;
  f = fopen(DUMP, "w");
  if (!f) {
    fprintf(stderr, "could not write checkpoint %llx\n", i8);
    return;
  }
  fprintf(f, "i8: %lld\n", i8);
  fprintf(f, "mem:\n");
  for (i=0; i<SIZE; ++i)
    fprintf(f, "  %.8x\n", mem[i]);
  fprintf(f, "rsl:\n");
  for (i=0; i<SIZE; ++i)
    fprintf(f, "  %.8x\n", rsl[i]);
  fprintf(f, "a: %ld\n", a);
  fprintf(f, "b: %ld\n", b);
  fprintf(f, "z: %ld\n", z);
  fprintf(f, "mem2:\n");
  for (i=0; i<NNN; ++i)
    fprintf(f, "  %.8lx\n", mem2[i]);
  fprintf(f, "rsl2:\n");
  for (i=0; i<NNN; ++i)
    fprintf(f, "  %.8lx\n", rsl2[i]);
  fprintf(f, "a2: %ld\n", a2);
  fprintf(f, "b2: %ld\n", b2);
  fprintf(f, "ss:\n");
  for (i=0; i<SIZE; ++i)
    fprintf(f, "  %.8lx\n", ss[i]);
  fprintf(f, "tt:\n");
  for (i=0; i<SIZE; ++i)
    fprintf(f, "  %.8lx\n", tt[i]);
  fprintf(f, "unidata:\n");
  for (i=0; i<USIZE; ++i)
    fprintf(f, "  %.16llx\n", unidata[i]);
  fprintf(f, "gapdata:\n");
  for (i=0; i<GSIZE; ++i)
    fprintf(f, "  %.16llx\n", gapdata[i]);
  fprintf(f, "rundata:\n");
  for (i=0; i<RSIZE; ++i)
    fprintf(f, "  %.16llx\n", rundata[i]);
  fprintf(f, "run2data:\n");
  for (i=0; i<RSIZE; ++i)
    fprintf(f, "  %.16llx\n", run2data[i]);
  fprintf(f, "njkdata:\n");
  for (i=0; i<USIZE; ++i)
    fprintf(f, "  %.16llx\n", njkdata[i]);
  fprintf(f, "gap2data:\n");
  for (i=0; i<GSIZE; ++i)
    fprintf(f, "  %.16llx\n", gap2data[i]);
  fprintf(f, "cycm:\n");
  for (i=0; i<SIZE; ++i)
    fprintf(f, "  %.8lx\n", cycm[i]);
  fprintf(f, "cyca: %ld\n", cyca);
  fprintf(f, "cycb: %ld\n", cycb);
  fprintf(f, "cycz: %ld\n", cycz);
  fprintf(f, "gaplast:\n");
  for (i=0; i<USIZE; ++i)
    fprintf(f, "  %.8lx\n", gaplast[i]);
  fprintf(f, "njkdata:\n");
  for (i=0; i<USIZE; ++i)
    fprintf(f, "  %.8lx\n", gap2last[i]);
  fclose(f);
}

void myload(u8 *i8, u4 *mem, u4 *rsl, u4 *a, u4 *b, u4 *z, 
	    u4 *mem2, u4 *rsl2, u4 *a2, u4 *b2, u4 *ss, u4 *tt, 
	    u8 *unidata, u8 *gapdata, u8 *rundata, u8 *run2data, 
	    u8 *njkdata, u8 *gap2data, 
	    u4 *cycm, u4 *cyca, u4 *cycb, u4 *cycz, 
	    u4 *gaplast, u4 *gap2last)
{
  u4 i;
  FILE *f;
  f = fopen(DUMP, "r");
  if (!f) {
    fprintf(stderr, "could not read checkpoint %s\n", DUMP);
    return;
  }
  fscanf(f, "i8: %lld\n", i8);
  fscanf(f, "mem:\n");
  for (i=0; i<SIZE; ++i)
    fscanf(f, "  %lx\n", &mem[i]);
  fscanf(f, "rsl:\n");
  for (i=0; i<SIZE; ++i)
    fscanf(f, "  %lx\n", &rsl[i]);
  fscanf(f, "a: %ld\n", a);
  fscanf(f, "b: %ld\n", b);
  fscanf(f, "z: %ld\n", z);
  fscanf(f, "mem2:\n");
  for (i=0; i<NNN; ++i)
    fscanf(f, "  %lx\n", &mem2[i]);
  fscanf(f, "rsl2:\n");
  for (i=0; i<NNN; ++i)
    fscanf(f, "  %lx\n", &rsl2[i]);
  fscanf(f, "a2: %ld\n", a2);
  fscanf(f, "b2: %ld\n", b2);
  fscanf(f, "ss:\n");
  for (i=0; i<SIZE; ++i)
    fscanf(f, "  %x\n", &ss[i]);
  fscanf(f, "tt:\n");
  for (i=0; i<SIZE; ++i)
    fscanf(f, "  %lx\n", &tt[i]);
  fscanf(f, "unidata:\n");
  for (i=0; i<USIZE; ++i)
    fscanf(f, "  %llx\n", &unidata[i]);
  fscanf(f, "gapdata:\n");
  for (i=0; i<GSIZE; ++i)
    fscanf(f, "  %llx\n", &gapdata[i]);
  fscanf(f, "rundata:\n");
  for (i=0; i<RSIZE; ++i)
    fscanf(f, "  %llx\n", &rundata[i]);
  fscanf(f, "run2data:\n");
  for (i=0; i<RSIZE; ++i)
    fscanf(f, "  %llx\n", &run2data[i]);
  fscanf(f, "njkdata:\n");
  for (i=0; i<USIZE; ++i)
    fscanf(f, "  %llx\n", &njkdata[i]);
  fscanf(f, "gap2data:\n");
  for (i=0; i<GSIZE; ++i)
    fscanf(f, "  %llx\n", &gap2data[i]);
  fscanf(f, "cycm:\n");
  for (i=0; i<SIZE; ++i)
    fscanf(f, "  %lx\n", &cycm[i]);
  fscanf(f, "cyca: %ld\n", cyca);
  fscanf(f, "cycb: %ld\n", cycb);
  fscanf(f, "cycz: %ld\n", cycz);
  fscanf(f, "gaplast:\n");
  for (i=0; i<USIZE; ++i)
    fscanf(f, "  %lx\n", &gaplast[i]);
  fscanf(f, "njkdata:\n");
  for (i=0; i<USIZE; ++i)
    fscanf(f, "  %lx\n", &gap2last[i]);
  fclose(f);
}

/*
------------------------------------------------------------------
Control -- initialize counters, run tests, and report results
------------------------------------------------------------------
*/
static void driver(u4 need_load)
{
  u4     i,expected;
  u8     i8;
  u4     a=2,b=2,z=2,a2=2,b2=2;
  double actual;

  for (i=0; i<NNN; ++i)   mem2[i] = rsl2[i] = i;
  for (i=0; i<SIZE; ++i)  ss[i] = tt[i] = mem2[i] = mem[i];
  for (i=0; i<USIZE; ++i) unidata[i] = gaplast[i] = gap2last[i] = 0;
  for (i=0; i<USIZE; ++i) njkdata[i] = 0;
  for (i=0; i<GSIZE; ++i) gapdata[i] = gap2data[i] = 0;
  for (i=0; i<RSIZE; ++i) rundata[i] = 0;
  for (i=0; i<RSIZE; ++i) run2data[i] = 0;
  gapoff = gap2off = 1;
  rcount = rlast = 0;

  /* generate billions of bits */
  rng( mem, rsl, &a, &b, &z, mem2, rsl2, &a2, &b2, ss, tt);
  for (i=0; i<10; ++i) 
  {
    prng_loop(mem2,rsl2, &a2, &b2);
    rng( mem, rsl, &a, &b, &z, mem2, rsl2, &a2, &b2, ss, tt);
  }
  for (i=0; i<SIZE; ++i) cycm[i] = mem[i];
  cyca = a; cycb = b; cycz = z;
  i8 = 0;
  if (need_load)
  {
      myload(&i8, mem, rsl, &a, &b, &z, mem2, rsl2, &a2, &b2, ss, tt, 
	     unidata, gapdata, rundata, run2data, njkdata, gap2data, 
	     cycm, &cyca, &cycb, &cycz, gaplast, gap2last);
      mydump(i8, mem, rsl, a, b, z, mem2, rsl2, a2, b2, ss, tt, 
	     unidata, gapdata, rundata, run2data, njkdata, gap2data, 
	     cycm, cyca, cycb, cycz, gaplast, gap2last);
  }
  for (; i8< MYRUNS; ++i8) 
  {
    if ((i8&0xffffffff) == 0xffffffff) {
      mydump(i8, mem, rsl, a, b, z, mem2, rsl2, a2, b2, ss, tt, 
	     unidata, gapdata, rundata, run2data, njkdata, gap2data, 
	     cycm, cyca, cycb, cycz, gaplast, gap2last);
    }
    rng( mem, rsl, &a, &b, &z, mem2, rsl2, &a2, &b2, ss, tt);
    uniform( rsl, SIZE, unidata, USIZE);
    gap( rsl, SIZE, gapdata, GSIZE);
    run( rsl, SIZE, rundata, RSIZE);
    run2( rsl, SIZE, run2data, RSIZE);
#ifdef NEVER
    njk( rsl, SIZE, njkdata, USIZE, gap2data, GSIZE); 
    if (cyc( mem, rsl, a, b, z, mem2, rsl2, &a2, &b2, ss, tt, i8)) {
      printf("cycle: ");
      break;
    }
#endif
  }
  printf("length %lld\n", i8);

  chisquare( "frequency", unidata, uniprob, USIZE);
  chisquare( "gap      ", gapdata, gapprob, GSIZE);
  chisquare( "run      ", rundata, runprob, RSIZE);
  chisquare( "run2     ", run2data, runprob, RSIZE);
#ifdef NEVER
  chisquare( "corr freq", njkdata, uniprob, USIZE);
  chisquare( "corr gap ", gap2data, gapprob, GSIZE);
#endif
}

/*
------------------------------------------------------------------
Try different seeds
------------------------------------------------------------------
*/
int main( int argc, char **argv)
{
  u4     i,j,k,expected, need_load = 0;
  double actual;
  need_load = (argc == 2 && strcmp(argv[1], DUMP) == 0);

  uprob( uniprob, USIZE);
  gprob( gapprob, GSIZE);
  rprob( runprob, USIZE);
  for (i=0; i<SIZE; ++i) mem[i] = (i+1)&UMASK;
  for (jjj=0; jjj<1; ++jjj)
    driver(need_load);
  return 1;
}
