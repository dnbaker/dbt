#ifndef NEW_SCAN_H__
#define NEW_SCAN_H__
#include "klib/khash.h"
#include <string>
#include <mutex>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <sys/stat.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <vector>
#include <map>
#ifdef GZSTREAM
#include <gzstream.h>
#endif
extern "C" {
#include "utils.h"
#include "xerrors.h"
}
static std::mutex map_mutex;

// =============== algorithm limits =================== 
// maximum number of distinct words
#define MAX_DISTINCT_WORDS (INT64_MAX -1)
typedef uint64_t word_int_t;
// maximum number of occurrences of a single word
#define MAX_WORD_OCC (UINT64_MAX)
typedef uint64_t occ_int_t;

// values of the wordFreq map: word, its number of occurrences, and its rank
struct word_stats {
  char *str;
  occ_int_t occ;
  word_int_t rank=0;
};
KHASH_MAP_INIT_INT64(stats, word_stats)
using std::string;
using std::cout;
using std::cerr;

// -------------------------------------------------------------
// struct containing command line parameters and other globals
struct Args {
   string inputFileName = "";
   string parse0ext = EXTPARS0;    // extension tmp parse file 
   string parseExt =  EXTPARSE;    // extension final parse file  
   string occExt =    EXTOCC;      // extension occurrences file  
   string dictExt =   EXTDICT;     // extension dictionary file  
   string lastExt =   EXTLST;      // extension file containing last chars   
   string saExt =     EXTSAI;      // extension file containing sa info   
   int w = 10;            // sliding window size and its default 
   int p = 100;           // modulus for establishing stopping w-tuples 
   bool SAinfo = false;   // compute SA information 
   int th=0;              // number of helper threads
   int verbose=0;         // verbosity level 
};

struct mt_data {
  khash_t(stats) *wordFreq; // shared dictionary
  Args *arg;       // command line input
  long start, end; // input
  long skipped, parsed, words;  // output
  FILE *parse, *last, *sa;
};

template<typename T, typename T2>
const T modexp(T base, T exp, T2 mod) {
	T2 result = 1;
	while (exp) {
		if (exp & 1) result = (result * base) % mod;
		exp >>= 1;
		base *= base;
	}
	return result;
}
static void die(const string s) {
  perror(s.c_str());
  exit(1);
}
static void save_update_word(string& w, unsigned int minsize, khash_t(stats) *freq, FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos);
// -----------------------------------------------------------------
// class to maintain a window in a string and its KR fingerprint
struct KR_window {
  int wsize;
  int asize;
  uint8_t *window;
  static constexpr uint64_t prime = 1999999973;
  uint64_t hash;
  uint64_t tot_char;
  uint64_t asize_pot;   // asize^(wsize-1) mod prime 
  
  KR_window(int w): wsize(w), asize(256), window(new int[w]), asize_pot(modexp(asize, wsize - 1, prime)) {
    // alloc and clear window
    reset();     
  }
  int       *begin()       {return window;}
  const int *begin() const {return window;}
  int       *end()       {return window + wsize;}
  const int *end() const {return window + wsize;}
  
  // init window, hash, and tot_char 
  void reset() {
    std::fill(begin(), end(), 0);
    // init hash value and related values
    hash=tot_char=0;    
  }
  
  uint64_t addchar(int c) {
    int k = tot_char++ % wsize;
    // complex expression to avoid negative numbers 
    hash += (prime - (window[k]*asize_pot) % prime); // remove window[k] contribution  
    window[k] = c;
    hash = (asize*hash + c) % prime;      //  add char i 
    // cerr << get_window() << " ~~ " << window << " --> " << hash << endl;
    return hash; 
  }
  // debug only 
  string get_window() const {
    string w;
    int k = (tot_char-1) % wsize;
    for(int i=k+1;i<k+1+wsize;i++)
      w.append(1,window[i%wsize]);
    return w;
  }
  
  ~KR_window() {
    delete[] window;
  } 

};


void *mt_parse(void *dx)
{
  // extract input data
  mt_data *d = (mt_data *) dx;
  Args *arg = d->arg;
  khash_t(stats) *wordFreq = d->wordFreq;
  uint64_t pos;

  if(arg->verbose>1)
    fprintf(stderr, "Scanning from %ld, size %ld\n",d->start,d->end-d->start);

  // open input file
  std::ifstream f(arg->inputFileName);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open file " + arg->inputFileName);
  }

  // prepare for parsing
  f.seekg(d->start); // move to the begining of assigned region
  KR_window krw(arg->w);
  int c; string word = "";
  d->skipped = d->parsed = d->words = 0;
  if(d->start==0) {
    word += static_cast<char>(Dollar); // no need to reach the next kr-window
  }
  else {   // reach the next breaking window
    while( (c = f.get()) != EOF ) {
      if(c<=Dollar) die("Invalid char found in input file. Exiting...");
      if(d->start + ++d->skipped == d->end + arg->w) {goto end;}
      word += static_cast<char>(c);
      uint64_t hash = krw.addchar(c);
      if(hash%arg->p==0 && d->skipped >= arg->w) break;
    }
    if(c==EOF) {f.close(); goto end;} // reached EOF without finding a breaking point nothing to do
    d->parsed = arg->w;   // the kr-window is part of the next word
    d->skipped -= arg->w; // ... so w less chars have been skipped
    word.erase(0,word.size() - arg->w);// keep only the last w chars
  }
  // cout << "Skipped: " << d->skipped << endl;

  // there is some parsing to do
  {
    pos = d->start;             // ending position+1 in text of previous word
    if(pos>0) pos+= d->skipped+ arg->w;  // or 0 for the first word
    assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file
    while( (c = f.get()) != EOF ) {
      if(c<=Dollar) die("Invalid char found in input file. Exiting...");
      word += static_cast<char>(c);
      const uint64_t hash = krw.addchar(c);
      d->parsed++;
      if(hash%arg->p==0 && d->parsed>arg->w) {
        // end of word, save it and write its full hash to the output file
        // pos is the ending position+1 of previous word and is updated in the next call
        save_update_word(word,arg->w,wordFreq,d->parse,d->last,d->sa,pos);
        d->words++;
        if(d->start+d->skipped+d->parsed>=d->end+arg->w) goto end;
      }
    }
  }
  // end of file reached
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg->w,Dollar);
  save_update_word(word,arg->w,wordFreq,d->parse,d->last,d->sa,pos);
  // close input file and return
  end:
  return NULL;
}


// prefix free parse of file fnam. w is the window size, p is the modulus
// use a KR-hash as the word ID that is written to the parse file
#if !NOTHREADS
uint64_t mt_process_file(Args& arg, khash_t(stats) *wf)
{
  // get input file size
  std::ifstream f(arg.inputFileName, std::ifstream::ate);
  if(!f.is_open()) {
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " +arg.inputFileName);
  }
  long size = f.tellg();
  f.close();

  // prepare and execute threads
  assert(arg.th>0);
  pthread_t t[arg.th];
  mt_data td[arg.th];
  for(int i=0;i<arg.th;i++) {
    td[i].wordFreq = wf;
    td[i].arg = &arg;
    td[i].start = i*(size/arg.th); // range start
    td[i].end = (i+1==arg.th) ? size : (i+1)*(size/arg.th); // range end
    assert(td[i].end<=size);
    // open the 1st pass parsing file
    td[i].parse = open_aux_file_num(arg.inputFileName.c_str(),arg.parse0ext.c_str(),i,"wb");
    // open output file containing the char at position -(w+1) of each word
    td[i].last = open_aux_file_num(arg.inputFileName.c_str(),arg.lastExt.c_str(),i,"wb");
    // if requested open file containing the ending position+1 of each word
    td[i].sa = arg.SAinfo ?open_aux_file_num(arg.inputFileName.c_str(),arg.saExt.c_str(),i,"wb") : NULL;
    xpthread_create(&t[i],NULL,&mt_parse,&td[i]);
  }

#if 0
  // open the 1st pass parsing file
  FILE *parse = open_aux_file(arg.inputFileName.c_str(),arg.parse0ext.c_str(),"wb");
  // open output file containing the char at position -(w+1) of each word
  FILE *last = open_aux_file(arg.inputFileName.c_str(),arg.lastExt.c_str(),"wb");
  // if requested open file containing the ending position+1 of each word
  FILE *sa = NULL;
  if(arg.SAinfo)
    sa = open_aux_file(arg.inputFileName.c_str(),arg.saExt.c_str(),"wb");
#endif

  // wait for the threads to finish (in order) and close output files
  long tot_char=0;
  for(int i=0;i<arg.th;i++) {
    xpthread_join(t[i],NULL);
    if(arg.verbose) {
      cout << "s:" << td[i].start << "  e:" << td[i].end << "  pa:";
      cout << td[i].parsed << "  sk:" << td[i].skipped << "  wo:" << td[i].words << '\n';
    }
    fclose(td[i].parse);
    fclose(td[i].last);
    if(td[i].sa) fclose(td[i].sa);
    if(td[i].words>0) {
      // extra check
      assert(td[i].parsed>arg.w);
      tot_char += td[i].parsed - (i!=0? arg.w: 0); //parsed - overlapping
    }
    else assert(i>0); // the first thread must produce some words
  }
  assert(tot_char==size);
  // close output files
  #if 0
  if(sa) if(fclose(sa)!=0) die("Error closing SA file");
  if(fclose(last)!=0) die("Error closing last file");
  if(fclose(parse)!=0) die("Error closing parse file");
  #endif
  return size;
}
#endif /* #if !NOTHREADS */

#endif /* #ifndef NEW_SCAN_H__ */
