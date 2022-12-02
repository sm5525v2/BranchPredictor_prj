/***************************************************************************************
*    Title: BATAGE simulator
*    Author: Pierre Michaud
*    Date: may 2018
*    Code version: 1
*    Availability: https://team.inria.fr/pacap/members/pierre-michaud/
*
***************************************************************************************/


// Author: Pierre Michaud
// Release 1, may 2018

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <algorithm>


#include <stdint.h>
#include <vector>

#define SIMFASTER

#ifdef SIMFASTER
#define ASSERT(cond)
#else
#define ASSERT(cond) if (!(cond)) {fprintf(stderr,"file %s assert line %d\n",__FILE__,__LINE__); abort();}
#endif


// ARM instructions are 4-byte aligned ==> shift PC by 2 bits
#define PC_SHIFT 2

// NUMG: number of tagged banks (also, number of non-null global history lengths)
// MAXHIST: greatest global history length
// MINHIST: smallest non-null global history length
// MAXPATH: path history length
// PATHBITS: number of target bits pushed into the path history per branch
// LOGB: log2 of the number of bimodal prediction bits
// LOGB2: log2 of the number of bimodal hysteresis entries
// LOGG: log2 of the number of entries in one tagged bank
// TAGBITS: tag size
// BHYSTBITS: number of bimodal hysteresis bits 
#define NUMG 26
#define MAXHIST 700
#define MINHIST 4
#define MAXPATH 30
#define PATHBITS 6
#define LOGB 12
#define LOGB2 10
#define LOGG 7
#define TAGBITS 12
#define BHYSTBITS 2

// SKIPMAX: maximum number of banks skipped on allocation
// if you change NUMG, you must re-tune SKIPMAX
#define SKIPMAX 5

// meta predictor, for a small accuracy boost (not in the BATAGE paper)
// inspired from the meta predictor in TAGE, but used differently
#define USE_META

// controlled allocation throttling (CAT)
// CATR = CATR_NUM / CATR_DEN
// CATR, CATMAX and MINAP must be re-tuned if you change NUMG or LOGG 
// for NUMG<~20, a smaller NUMG needs a smaller CATR
// a larger predictor may need a larger CATMAX
#define CATR_NUM 2
#define CATR_DEN 3
#define CATMAX ((CATR_NUM<<15)-1)
#define MINAP 16

// controlled decay, for a tiny accuracy boost (not in the BATAGE paper)
// CDR = CDR_NUM / CDR_DEN
// CDR must be greater than CATR
// CDR, CDMAX and MINDP must be re-tuned if you change NUMG or LOGG 
// in particular, if you decrease CATR, you should probably decrease CDR too
#define USE_CD
#ifdef USE_CD
#define CDR_NUM 6
#define CDR_DEN 5
#define CDMAX ((CDR_NUM<<10)-1)
#define MINDP 7
#endif

// bank interleaving, inspired from Seznec's TAGE-SC-L (CBP 2016)
#define BANK_INTERLEAVING
#ifdef BANK_INTERLEAVING
// taking MIDBANK=(NUMG-1)*0.4 is probably close to optimal
// take GHGBITS=1 for NUMG<~10
#define MIDBANK 10
#define GHGBITS 2
#endif


using namespace std;


uint32_t
rando()
{
  // Marsaglia's xorshift
  static uint32_t x = 2463534242;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  return x;
}


int
gcd(int a, int b)
{
  if (a == b) {
    return a;
  } else if (a > b) {
    return gcd(b,a-b);
  } else{
    return gcd(a,b-a);
  }
}


int
lcm(int a, int b)
{
  int g = gcd(a,b);
  ASSERT(g>0);
  return (a*b)/g;
}


uint32_t
reverse(uint32_t x, int nbits)
{
  // reverse the 16 rightmost bit (Strachey's method)
  ASSERT(nbits <= 16);
  x &= 0xFFFF;
  x |= (x & 0x000000FF) << 16;
  x = (x & 0xF0F0F0F0) | ((x & 0x0F0F0F0F) << 8);
  x = (x & 0xCCCCCCCC) | ((x & 0x33333333) << 4);
  x = (x & 0xAAAAAAAA) | ((x & 0x55555555) << 2);
  return x >> (31-nbits);
}

class dualcounter {
 public:
  static const int nmax = 7;
  int n[2];
dualcounter(int b1, int b2)
{
  ASSERT(BHYSTBITS>=1);
  n[b1^1] = (1<<(BHYSTBITS-1)) - 1;
  n[b1] = n[b1^1] + (1<<BHYSTBITS) - b2;
}


dualcounter()
{
  reset();
}


void
reset()
{
  n[0] = n[1] = 0;
}


int
pred()
{
  return (n[1] > n[0]);
}


bool
mediumconf()
{
  return (n[1] == (2*n[0]+1)) || (n[0] == (2*n[1]+1));
}


bool
lowconf()
{
  return (n[1] < (2*n[0]+1)) && (n[0] < (2*n[1]+1));
}


bool
highconf()
{
  return ! mediumconf() && ! lowconf();
}


bool
veryhighconf()
{
  ASSERT(nmax==7);
  // the formula below works with nmax=7 (see the BATAGE paper)
  return ((n[0]==0) && (n[1]>=4)) || ((n[1]==0) && (n[0]>=4));
}


int
sum()
{
  return n[1] + n[0];
}


int
diff()
{
  return abs(n[1]-n[0]);
}


bool
saturated()
{
  return (n[0] == nmax) || (n[1] == nmax);
}


int
conflevel(int meta)
{
  if (n[1]==n[0]) return 3;
  bool low = lowconf();
  bool promotemed = (meta>=0) && (sum()==1);
  bool med = ! promotemed && mediumconf();
  return (low<<1) | med;
}


void
decay()
{
  if (n[1] > n[0]) n[1]--;
  if (n[0] > n[1]) n[0]--;
}


void
update(bool taken)
{
  ASSERT((n[0]>=0) && (n[0]<=nmax));
  ASSERT((n[1]>=0) && (n[1]<=nmax));
  int d = (taken)? 1:0;
  if (n[d] < nmax) {
    n[d]++;
  } else if (n[d^1] > 0) {
    n[d^1]--;
  }
}


static int
size()
{
  unsigned x = nmax;
  int nbits = 0;
  while (x) {
    nbits++;
    x>>=1;
  }
  return 2*nbits;
}
};


class path_history {
public:
  int ptr; 
  int hlength;
  unsigned * h;
void
init(int hlen)
{
  hlength = hlen;
  h = new unsigned [hlen];
  for (int i=0; i<hlength; i++) {
    h[i] = 0;
  }
  ptr = 0;
}


void 
insert(unsigned val)
{
  ptr--;
  if (ptr == (-1)) {
    ptr = hlength-1;
  }
  h[ptr] = val;
}


unsigned & 
operator [] (int n)
{
  ASSERT((n>=0) && (n<hlength));
  int k = ptr + n;
  if (k >= hlength) {
    k -= hlength;
  }
  ASSERT((k>=0) && (k<hlength));
  return h[k];
}
};


// folded global history, for speeding up hash computations
// I introduced it in the 2001 tech report "A Comprehensive Study of Dynamic Global History Branch Prediction" (https://hal.inria.fr/inria-00072400)
// I advertised it with the PPM-like predictor (CBP 2004)
// see also the paper by Schlais and Lipasti at ICCD 2016

class folded_history {
public:
  uint32_t fold;
  int clength;
  int olength;
  int outpoint;
  uint32_t mask1;
  uint32_t mask2;
void 
init(int original_length, int compressed_length, int injected_bits)
{
  olength = original_length;
  clength = compressed_length;
  outpoint = olength % clength;
  ASSERT(clength < 32);
  injected_bits = min(injected_bits,clength);
  mask1 = (1<<clength)-1;
  mask2 = (1<<injected_bits)-1;
  fold = 0; // must be consistent with path_history::init()
}


uint32_t
rotateleft(uint32_t x, int m)
{
  ASSERT(m < clength);
  ASSERT((x>>clength) == 0);
  uint32_t y = x >> (clength-m);
  x = (x << m) | y;
  return x & mask1;
}


void 
update(path_history & ph)
{
  fold = rotateleft(fold,1);
  unsigned inbits = ph[0] & mask2;
  unsigned outbits = ph[olength] & mask2;
  outbits = rotateleft(outbits,outpoint);
  fold ^= inbits ^ outbits;
}
};


class histories {
 public:
  path_history bh; // global history of branch directions
  path_history ph; // path history (target address bits)
  folded_history * chg;
  folded_history * chgg;
  folded_history * cht;
  folded_history * chtt;
histories()
{
  int * hist = new int [NUMG];
  int prevh = 0;
  for (int i=0; i<NUMG; i++) {
    // geometric history lengths
    int h = MINHIST * pow((double)MAXHIST/MINHIST,(double)i/(NUMG-1));
    h = max(prevh+1,h);
    hist[NUMG-1-i] = h; // hist[0] = longest history
    prevh = h;
  }
  bh.init(MAXHIST+1);
  ph.init(MAXPATH+1);
  chg = new folded_history [NUMG];
  chgg = new folded_history [NUMG];
  cht = new folded_history [NUMG];
  chtt = new folded_history [NUMG];
  for (int i=0; i<NUMG; i++) {
    chg[i].init(hist[i],LOGG,1);
    cht[i].init(hist[i],TAGBITS,1);
    int hashparam = 1;
    if (LOGG == TAGBITS) {
      hashparam = (lcm(LOGG,LOGG-3) > lcm(LOGG,LOGG-2))? 3 : 2;
    }
    if (hist[i] <= MAXPATH) {
      chgg[i].init(hist[i],LOGG-hashparam,PATHBITS);
      chtt[i].init(hist[i],TAGBITS-1,PATHBITS);
    } else {
      chgg[i].init(hist[i],LOGG-hashparam,1);
      chtt[i].init(hist[i],TAGBITS-1,1);
    }
  }
}


void 
update(uint32_t targetpc, bool taken)
{
#ifdef PC_SHIFT
  targetpc ^= targetpc >> PC_SHIFT;
#endif
  bh.insert(taken);
  ph.insert(targetpc);
  for (int i=0; i<NUMG; i++) {
    chg[i].update(bh);
    cht[i].update(bh);
    if (chgg[i].olength <= MAXPATH) {
      chgg[i].update(ph);
      chtt[i].update(ph);
    } else {
      chgg[i].update(bh);
      chtt[i].update(bh);
    }      
  }
}

// the hash functions below are somewhat more complex than what would be implemented in real processors

int
gindex(uint32_t pc, int i)
{
  ASSERT((i>=0) && (i<NUMG));
  uint32_t hash = pc ^ i ^ chg[i].fold ^ (chgg[i].fold << (chg[i].clength-chgg[i].clength));
  return hash & ((1<<LOGG)-1);
}


int
gtag(uint32_t pc, int i)
{
  ASSERT((i>=0) && (i<NUMG));
  uint32_t hash = (pc+i) ^ reverse(cht[i].fold,cht[i].clength) ^ (chtt[i].fold << (cht[i].clength-chtt[i].clength));
#ifdef GHGBITS
  // when bank interleaving is enabled, introducing 1 or 2 bits in the tag
  // for identifying the path length generally reduces tag aliasing when NUMG is large
  hash = (hash << GHGBITS) | ghg(i);
#endif
  return hash & ((1<<TAGBITS)-1);
}


#ifdef BANK_INTERLEAVING
// inspired from Seznec's TAGE-SC-L (CBP 2016), but different:
// interleaving is global, unlike in TAGE-SC-L ==> unique tag size
int
phybank(int i)
{
  ASSERT((i>=0) && (i<NUMG));
  unsigned pos;
  if (i >= (NUMG-MIDBANK)) {
    pos = i;
  } else {
    // on some workloads, the shortest non-null global history length does not generate
    // enough entropy, which may lead to uneven utilization of banks, hence MIDBANK
    pos = (chgg[NUMG-MIDBANK-1].fold + i) % (NUMG-MIDBANK);
  }
  return (chgg[NUMG-1].fold + pos) % NUMG;
}


#ifdef GHGBITS
int
ghg(int i)
{
  return ((NUMG-1-i) << GHGBITS) / NUMG;
}
#endif

#endif // BANK_INTERLEAVING


void
printconfig()
{
  printf("history lengths: ");
  for (int i=NUMG-1; i>=0; i--) {
    printf("%d ",chg[i].olength);
  }
  printf("\n");
}


int
size()
{
  return MAXHIST + MAXPATH * PATHBITS; // number of bits
  // The storage for folded histories is ignored here
  // If folded histories are implemented in hardware, they must be checkpointed (cf. Schlais & Lipasti, ICCD 2016)
  // An implementation without folded histories is possible (cf. Seznec's GEHL, ISCA 2005)
}
};


// tagged entry
class tagged_entry {
 public:
  int tag;
  dualcounter dualc; 
  tagged_entry()
{
  tag = 0;
  dualc.reset();
}
};


class batage {
 public:
  int b [1<<LOGB];   // bimodal predictions
  int b2 [1<<LOGB2]; // bimodal hystereses
  // eg) 2 bit msb: prediction, lsb: hyeteresis 11, 10, (prediction 1) 01, 00 (prediction 0) 

  tagged_entry ** g; // tagged entries
  int bi; // hash for the bimodal prediction
  int bi2; // hash for the bimodal hysteresis
  int * gi; // hashes for the tagged banks
  vector<int> hit; // tell which banks have a hit
  vector<dualcounter> s; // dual-counters for the hitting banks
  int bp; // dual-counter providing the final BATAGE prediction
  int cat; // CAT counter
  int meta; // for a small accuracy boost
#ifdef USE_CD
  int cd; // for a small accuracy boost
#endif
#ifdef BANK_INTERLEAVING
  int * bank;
  bool * check;
#endif
 batage()
{
  g = new tagged_entry * [NUMG];
  for (int i=0; i<NUMG; i++) {
    g[i] = new tagged_entry [1<<LOGG];
  }
  gi = new int [NUMG];
  for (int i=0; i<(1<<LOGB); i++) {
    b[i] = 0; // not-taken prediction
  }
  for (int i=0; i<(1<<LOGB2); i++) {
    b2[i] = (1<<BHYSTBITS)-1; // weak state
  }
  cat = 0;
  meta = -1;
#ifdef USE_META
  meta = 0;
#endif
#ifdef USE_CD
  cd = 0;
#endif
#ifdef BANK_INTERLEAVING
  bank = new int [NUMG];
  check = new bool [NUMG]; 
#endif   
}


#ifdef BANK_INTERLEAVING
void
check_bank_conflicts()
{
  for (int i=0; i<NUMG; i++) {
    check[i] = false;
  }
  for (int i=0; i<NUMG; i++) {
    if (check[bank[i]]) {
      fprintf(stderr,"BANK CONFLICT\n");
      ASSERT(0);
    }
    check[bank[i]] = true;
  }
}
#endif 


tagged_entry &
getg(int i)
{
  ASSERT((i>=0) && (i<NUMG));
#ifdef BANK_INTERLEAVING
  ASSERT((bank[i]>=0) && (bank[i]<NUMG));
  return g[bank[i]][gi[i]];
#else
  return g[i][gi[i]];
#endif
}



bool 
predict(uint32_t pc, histories & p)
{
#ifdef PC_SHIFT
  pc ^= pc >> PC_SHIFT;
#endif
  
  hit.clear();
  s.clear();
  for (int i=0; i<NUMG; i++) {
#ifdef BANK_INTERLEAVING
    bank[i] = p.phybank(i); //for bank interleaving, inspired from Seznec's TAGE-SC-L (CBP 2016)
#endif
    gi[i] = p.gindex(pc,i); //gi: hashed for tagged banks, gindex return index using address, history, i
    if (getg(i).tag == p.gtag(pc,i)) { //getg get tagged_entry
      hit.push_back(i);
      s.push_back(getg(i).dualc); 
      //getg returns g[i][gi[i]]; i = length of history, g[i] == length i history table?
      //push dualCounter of g[i][gi[i]] tagged_entry to s
    }
  }
  //now s has all dual counter in tagged entries
  
#ifdef BANK_INTERLEAVING
#ifndef SIMFASTER
  check_bank_conflicts();
#endif
#endif

  //int bi : hash for the bimodal prediction
  //int bi2 : hash for the bimodal hysteresis
  //then why push it to s making dualcounter?
  bi = pc & ((1<<LOGB)-1);
  bi2 = bi & ((1<<LOGB2)-1);
  s.push_back(dualcounter(b[bi],b2[bi2]));

  bp = 0;
  for (int i=1; i<(int)s.size(); i++) {
    if (s[i].conflevel(meta) < s[bp].conflevel(meta)) {
      bp = i;
      //find lowest conflevel index
    }
  }
  
  //return prediction ( = n[1] > n[0]) of lowest confident level dual counter
  return s[bp].pred();
}


void
update_bimodal(bool taken)
{

  // see Loh et al., "Exploiting bias in the hysteresis bit of 2-bit saturating counters in branch predictors", Journal of ILP, 2003   
  if (b[bi] == taken) {
    if (b2[bi2]>0) b2[bi2]--;
  } else {
    if (b2[bi2] < ((1<<BHYSTBITS)-1)) {
      b2[bi2]++;
    } else {
      b[bi] = taken;
    }
  }  
}


void
update_entry(int i, bool taken)
{
  ASSERT(i<s.size());
  if (i<(int)hit.size()) {
    getg(hit[i]).dualc.update(taken);
  } else {
    update_bimodal(taken);
  }
}


void
update(uint32_t pc, bool taken, histories & p, bool noalloc = false)
{
#ifdef PC_SHIFT
  pc ^= pc >> PC_SHIFT;
#endif

#ifdef USE_META
  if ((s.size()>1) && (s[0].sum()==1) && s[1].highconf() && (s[0].pred() != s[1].pred())) {
    if (s[0].pred() == taken) {
      if (meta<15) meta++;
    } else {
      if (meta>(-16)) meta--;
    }
  }
#endif
  
  //bp: previous prediction tagged length

  // update from 0 to bp-1
  for (int i=0; i<bp; i++) {
    if ((meta>=0) || s[i].lowconf() || (s[i].pred() != s[bp].pred()) || ((rando() % 8)==0)) {
      getg(hit[i]).dualc.update(taken);
    }
  }
  // update at bp
  if ((bp<(int)hit.size()) && s[bp].highconf() && s[bp+1].highconf() && (s[bp+1].pred()==taken) && ((s[bp].pred()==taken) || (cat>=(CATMAX/2)))) {
    if (! s[bp].saturated() || ((meta<0) && ((rando() % 8)==0))) { 
      getg(hit[bp]).dualc.decay();
    }
  } else {
    update_entry(bp,taken);
  }
  // update at bp+1
  if (! s[bp].highconf() && (bp<(int)hit.size())) {
    update_entry(bp+1,taken);
  }
 
  // ALLOCATE

  bool allocate = ! noalloc && (s[bp].pred() != taken);

  if (allocate && ((int)(rando() % MINAP) >= ((cat * MINAP) / (CATMAX+1)))) {
    int i = (hit.size()>0)? hit[0] : NUMG;
    i -= rando() % (1+s[0].diff()*SKIPMAX/dualcounter::nmax);
    int mhc = 0;
    while (--i >= 0) {
      if (getg(i).dualc.highconf()) { //if index entry is high confidence
    #ifdef USE_CD
      if ((int)(rando() % MINDP) >= ((cd*MINDP)/(CDMAX+1))) getg(i).dualc.decay();
    #else
      if ((rando() % 4)==0) getg(i).dualc.decay(); // decay dual counter
    #endif
      if (! getg(i).dualc.veryhighconf()) mhc++; // high but not very high confidence
      } else {  //not high confidence
        //reset tagged_entry
        getg(i).tag = p.gtag(pc,i); //get gtag from history and set to tagged_entry
        getg(i).dualc.reset(); //reset
        getg(i).dualc.update(taken);
        cat += CATR_NUM - mhc * CATR_DEN;
        cat = min(CATMAX,max(0,cat));
      #ifdef USE_CD
        cd += CDR_NUM - mhc * CDR_DEN;
        cd = min(CDMAX,max(0,cd));
      #endif
        break;
      }
    }
  }

}


int
size()
{
  int totsize = (1 << LOGB) + (BHYSTBITS << LOGB2);
  totsize += NUMG * ((dualcounter::size() + TAGBITS) << LOGG);
  return totsize; // number of bits
  // the storage for counters 'cat', 'meta' and 'cd' is neglected here
}
};


class batage_update : public branch_update {
public:
    //add variable if need
};

class batage_predictor : public branch_predictor{

 private:

  batage pred;
  histories hist;

  batage_update u; //predict.cc use branch_update to get result
  long long int PC; //PC is need to predict and update batage, so store it in member variable and pass to pred.predict method

 public:
  batage_predictor(void)
  {
    hist.printconfig();
    printf("total bits = %d\n",pred.size()+hist.size());
  }

  branch_update* predict (branch_info & b) {
      PC = b.address;
      u.direction_prediction(pred.predict(PC, hist));
      return &u;
  }

  void update (branch_update *u, bool taken, unsigned int target) {
      pred.update(PC,taken,hist, false);
      hist.update(target,taken);
  }
};