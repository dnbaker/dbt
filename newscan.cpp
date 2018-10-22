/* ******************************************************************************
 * newscan.cpp
 * 
 * parsing algorithm for the bwt construction for repetitive sequences based 
 * on prefix free parsing. See:
 *   Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini
 *   Prefix-Free Parsing for Building Big BWTs
 *   [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/)
 * 
 * Usage:
 *   newscan.x wsize modulus file
 * 
 * Accepts any kind of file that does not contain the chars 0x0, 0x1, 0x2 
 * which are used internally. If input file is gzipped use cnewscan.x which 
 * automatcally extracts the content
 * 
 * The parameters wsize and modulus are used to define the prefix free parsing 
 * using KR-fingerprints (see paper)
 * 
 * The algorithm computes the prefix free parsing of 
 *     T = (0x2)file_content(0x2)^wsize
 * in a dictionary of words D and a parsing P of T in terms of the  
 * dictionary words. Note that the words in the parsing overlap by wsize.
 * Let d denote the number of words in D and p the number of phrases in 
 * the parsing P
 * 
 * newscan outputs the following files:
 * 
 * file.dict
 * containing the dictionary words in lexicographic order with a 0x1 at the end of
 * each word and a 0x0 at the end of the file. Size: |D| + d + 1 where
 * |D| is the sum of the word lengths
 * 
 * file.occ
 * the number of occurrences of each word in lexicographic order.
 * We assume the number of occurrences of each word is at most 2^32-1
 * so the size is 4d bytes
 * 
 * file.parse
 * containing the parse P with each word identified with its 1-based lexicographic 
 * rank (ie its position in D). We assume the number of distinct words
 * is at most 2^32-1, so the size is 4p bytes
 * 
 * file.last 
 * contaning the charater in positon w+1 from the end for each dictionary word
 * Size: d
 * 
 * file.sai (if option -s is given on the command line) 
 * containing the ending position +1 of each dictionary word in the original
 * text written using IBYTES bytes for each entry 
 * Size: d*IBYTES
 * 
 * The output of newscan must be processed by bwtparse, which invoked as
 * 
 *    bwtparse file
 * 
 * computes the BWT of file.parse and produces file.ilist of size 4p+4 bytes
 * contaning, for each dictionary word in lexicographic order, the list 
 * of BWT positions where that word appears (ie i\in ilist(w) <=> BWT[i]=w).
 * There is also an entry for the EOF word which is not in the dictionary 
 * but is assumed to be the smallest word.  
 * 
 * In addition, bwtparse permutes file.last according to
 * the BWT permutation and generates file.bwlast such that file.bwlast[i] 
 * is the char from P[SA[i]-2] (if SA[i]==0 , BWT[i]=0 and file.bwlast[i]=0, 
 * if SA[i]==1, BWT[i]=P[0] and file.bwlast[i] is taken from P[n-1], the last 
 * word in the parsing).  
 * 
 * If the option -s is given to bwtparse, it permutes file.sai according
 * to the BWT permutation and generate file.bwsai using again IBYTES
 * per entry.  file.bwsai[i] is the ending position+1 ofBWT[i] in the 
 * original text 
 * 
 * The output of bwtparse (the files .ilist .bwlast) together with the
 * dictionary itself (file .dict) and the number of occurrences
 * of each word (file .occ) are used to compute the final BWT by the 
 * pfbwt algorithm.
 * 
 */
#include "newscan.hpp"
#include <cstring>

using namespace std;
using namespace __gnu_cxx;




// -----------------------------------------------------------



// compute 64-bit KR hash of a string 
// to avoid overflows in 64 bit aritmethic the prime is taken < 2**55
uint64_t kr_hash(string s) {
    uint64_t hash = 0;
    //const uint64_t prime = 3355443229;     // next prime(2**31+2**30+2**27)
    static constexpr uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    for(size_t k=0;k<s.size();k++) {
      int c = (unsigned char) s[k];
      assert(c>=0 && c< 256);
      hash = (256*hash + c) % prime;    //  add char k
    } 
    return hash; 
}
uint64_t kr_hash(const char *s) {
    static constexpr uint64_t prime = 27162335252586509; // next prime (2**54 + 2**53 + 2**47 + 2**13)
    uint64_t hash = (unsigned char) *s++;
    while(*s) {
        int c = (unsigned char) *s++;
        hash <<= 8;
        hash |= c;
        hash %= prime;
    }
    return hash; 
}



// save current word in the freq map and update it leaving only the 
// last minsize chars which is the overlap with next word  
static void save_update_word(string& w, unsigned int minsize, khash_t(stats) *freq, FILE *tmp_parse_file, FILE *last, FILE *sa, uint64_t &pos)
{
  assert(pos==0 || w.size() > minsize);
  if(w.size() <= minsize) return;
  // get the hash value and write it to the temporary parse file
  uint64_t hash = kr_hash(w);
  if(fwrite(&hash,sizeof(hash),1,tmp_parse_file)!=1) die("parse write error");

  // update frequency table for current hash
  khiter_t ki = kh_get(stats, freq, hash);
  if(ki == kh_end(freq)) {
      int khr;
#ifndef NOTHREADS
      std::lock_guard<std::mutex> lock(map_mutex);
#endif
      ki = kh_put(stats, freq, hash, &khr);
      kh_val(freq, ki).occ = 1;
      kh_val(freq, ki).rank = 0;
      kh_val(freq, ki).str = static_cast<char *>(std::malloc(w.size() + 1));
      std::memcpy(kh_val(freq, ki).str, w.data(), w.size() + 1);
  }
  else {
#ifndef NOTHREADS
      std::lock_guard<std::mutex> lock(map_mutex);
#endif
      if(__builtin_expect(++kh_val(freq, ki).occ <= 0, 0)) throw std::runtime_error(std::string("Emergency exit: Maximum number of occurrences ") + std::to_string(MAX_WORD_OCC) + " exceeded.");
      if(__builtin_expect(std::strcmp(kh_val(freq, ki).str, w.data()), 0)) {
        throw std::runtime_error(std::string("Emergency exit! Hash collision for strings: ") + kh_val(freq, ki).str);
      }
  }
  // output char w+1 from the end
  if(fputc(w[w.size()- minsize-1],last)==EOF) die("Error writing to .last file");
  // compute ending position +1 of current word and write it to sa file 
  // pos is the ending position+1 of the previous word and is updated here 
  if(pos==0) pos = w.size()-1; // -1 is for the initial $ of the first word
  else pos += w.size() -minsize; 
  if(sa) if(fwrite(&pos,IBYTES,1,sa)!=1) die("Error writing to sa info file");
  // keep only the overlapping part of the window
  w.erase(0,w.size() - minsize);
}



// prefix free parse of file fnam. w is the window size, p is the modulus 
// use a KR-hash as the word ID that is immediately written to the parse file
uint64_t process_file(Args& arg, khash_t(stats) *wordFreq)
{
  //open a, possibly compressed, input file
  string fnam = arg.inputFileName;
  #ifdef GZSTREAM 
  igzstream f(fnam.c_str());
  #else
  ifstream f(fnam);
  #endif    
  if(!f.rdbuf()->is_open()) {// is_open does not work on igzstreams 
    perror(__func__);
    throw new std::runtime_error("Cannot open input file " + fnam);
  }

  // open the 1st pass parsing file 
  FILE *g = open_aux_file(arg.inputFileName.c_str(),arg.parse0ext.c_str(),"wb");
  // open output file containing the char at position -(w+1) of each word
  FILE *last_file = open_aux_file(arg.inputFileName.c_str(),arg.lastExt.c_str(),"wb");  
  // if requested open file containing the ending position+1 of each word
  FILE *sa_file = NULL;
  if(arg.SAinfo) 
    sa_file = open_aux_file(arg.inputFileName.c_str(),arg.saExt.c_str(),"wb");
  
  // main loop on the chars of the input file
  int c;
  uint64_t pos = 0; // ending position +1 of previous word in the original text, used for computing sa_info 
  assert(IBYTES<=sizeof(pos)); // IBYTES bytes of pos are written to the sa info file 
  // init first word in the parsing with a Dollar char 
  string word("");
  word.append(1,Dollar);
  // init empty KR window: constructor only needs window size
  KR_window krw(arg.w);
  while( (c = f.get()) != EOF ) {
    if(c<=Dollar) {cerr << "Invalid char found in input file: no additional chars will be read\n"; break;}
    word += static_cast<char>(c);
    uint64_t hash = krw.addchar(c);
    if(hash%arg.p==0) {
      // end of word, save it and write its full hash to the output file
      // cerr << "~"<< c << "~ " << hash << " ~~ <" << word << "> ~~ <" << krw.get_window() << ">" <<  endl;
      save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
    }    
  }
  // virtually add w null chars at the end of the file and add the last word in the dict
  word.append(arg.w,Dollar);
  save_update_word(word,arg.w,wordFreq,g,last_file,sa_file,pos);
  // close input and output files 
  if(sa_file) if(fclose(sa_file)!=0) die("Error closing SA file");
  if(fclose(last_file)!=0) die("Error closing last file");  
  if(fclose(g)!=0) die("Error closing parse file");
  if(pos!=krw.tot_char+arg.w) cerr << "Pos: " << pos << " tot " << krw.tot_char << endl;
  f.close();
  return krw.tot_char;
}

// function used to compare two string pointers
bool pstringCompare(const string *a, const string *b)
{
  return *a <= *b;
}

// given the sorted dictionary and the frequency map write the dictionary and occ files
// also compute the 1-based rank for each hash
void writeDictOcc(Args &arg, khash_t(stats) *wfreq, vector<const char *> &sortedDict)
{
  assert(sortedDict.size() == kh_size(wfreq));
  // open dictionary and occ files 
  string fdictname = arg.inputFileName + "." + arg.dictExt;
  FILE *fdict = fopen(fdictname.c_str(),"wb");
  if(fdict==NULL) die("Cannot open " + fdictname);
  string foccname = arg.inputFileName + "." + arg.occExt;
  FILE *focc = fopen(foccname.c_str(),"wb");
  if(focc==NULL) die("Cannot open " + foccname);
  
  word_int_t wrank = 1; // current word rank (1 based)
  for(auto x: sortedDict) {
    size_t sl = std::strlen(x);
    size_t s = fwrite(x, 1, sl, fdict);
    if(s!=sl) die("Error writing to " + fdictname);
    if(fputc(EndOfWord,fdict)==EOF) die("Error writing EndOfWord to " + fdictname);
    uint64_t hash = kr_hash(x);
    khiter_t ki = kh_get(stats, wfreq, hash);
    if(ki == kh_end(wfreq)) throw 1; // TODO: better error message
    auto &wf = kh_val(wfreq, ki);
    assert(wf.occ>0);
    s = fwrite(&wf.occ,sizeof(wf.occ),1, focc);
    if(s!=1) die("Error writing to " + foccname);
    assert(wf.rank==0);
    wf.rank = wrank++;
  }
  if(fputc(EndOfDict,fdict)==EOF) die("Error writing EndOfDict to " + fdictname);
  if(fclose(focc)!=0) die("Error closing " + foccname);
  if(fclose(fdict)!=0) die("Error closing" + fdictname);
}

void remapParse(Args &arg, khash_t(stats) *wfreq)
{
  // open parse files. the old parse can be stored in a single file or in multiple files
  mFile *moldp = mopen_aux_file(arg.inputFileName.c_str(), arg.parse0ext.c_str(), arg.th);
  FILE *newp = open_aux_file(arg.inputFileName.c_str(), arg.parseExt.c_str(), "wb");

  // recompute occ as an extra check 
  vector<occ_int_t> occ(kh_size(wfreq)+1,0); // ranks are zero based 
  uint64_t hash;
  while(true) {
    size_t s = mfread(&hash,sizeof(hash),1,moldp);
    if(s==0) break;
    if(s!=1) die("Unexpected parse EOF");
    word_int_t rank;
    {
        auto it = kh_get(stats, wfreq, hash);
        rank = kh_val(wfreq, it).rank;
    }
    ++occ[rank];
    s = fwrite(&rank,sizeof(rank),1,newp);
    if(s!=1) die("Error writing to new parse file");
  }
  if(fclose(newp)!=0) die("Error closing new parse file");
  if(mfclose(moldp)!=0) die("Error closing old parse segment");
  // check old and recomputed occ coincide 
#if !NDEBUG
  for(khiter_t ki = 0; ki < kh_size(wfreq); ++ki)
    if(kh_exist(wfreq, ki))
        assert(kh_val(wfreq, ki).rank == occ[kh_val(wfreq, ki).rank]);
#endif
}
 



void print_help(char** argv, Args &args) {
  cout << "Usage: " << argv[ 0 ] << " <input filename> [options]" << endl;
  cout << "  Options: " << endl
        << "\t-w W\tsliding window size, def. " << args.w << endl
        << "\t-p M\tmodulo for defining phrases, def. " << args.p << endl
        #ifndef NOTHREADS
        << "\t-t M\tnumber of helper threads, def. none " << endl
        #endif        
        << "\t-h  \tshow help and exit" << endl
        << "\t-s  \tcompute suffix array info" << endl;
  #ifdef GZSTREAM
  cout << "If the input file is gzipped it is automatically extracted\n";
  #endif
  exit(1);
}

void parseArgs( int argc, char** argv, Args& arg ) {
   int c;
   extern char *optarg;
   extern int optind;

  puts("==== Command line:");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  puts("");

   string sarg;
   while ((c = getopt( argc, argv, "p:w:sht:v") ) != -1) {
      switch(c) {
        case 's':
        arg.SAinfo = true; break;
        case 'w':
        sarg.assign( optarg );
        arg.w = stoi( sarg ); break;
        case 'p':
        sarg.assign( optarg );
        arg.p = stoi( sarg ); break;
        case 't':
        sarg.assign( optarg );
        arg.th = stoi( sarg ); break;
        case 'v':
           arg.verbose++; break;
        case 'h':
           print_help(argv, arg); exit(1);
        case '?':
        cout << "Unknown option. Use -h for help." << endl;
        exit(1);
      }
   }
   // the only input parameter is the file name 
   if (argc == optind+1) {
     arg.inputFileName.assign( argv[optind] );
   }
   else {
      cout << "Invalid number of arguments" << endl;
      print_help(argv,arg);
   }
   // check algorithm parameters 
   if(arg.w <4) {
     cout << "Windows size must be at least 4\n";
     exit(1);
   }
   if(arg.p<10) {
     cout << "Modulus must be at leas 10\n";
     exit(1);
   }
   #ifdef NOTHREADS
   if(arg.th!=0) {
     cout << "The NT version cannot use threads\n";
     exit(1);
   }
   #else
   if(arg.th<0) {
     cout << "Number of threads cannot be negative\n";
     exit(1);
   }
   #endif   
}



int main(int argc, char** argv)
{
  // translate command line parameters
  Args arg;
  parseArgs(argc, argv, arg);
  cout << "Windows size: " << arg.w << endl;
  cout << "Stop word modulus: " << arg.p << endl;  

  // measure elapsed wall clock time
  time_t start_main = time(NULL);
  time_t start_wc = start_main;  
  // init sorted map counting the number of occurrences of each word
  khash_t(stats) *wordFreq = kh_init(stats);
  uint64_t totChar;

  // ------------ parsing input file 
  try {
      if(arg.th==0)
        totChar = process_file(arg,wordFreq);
      else {
        #ifdef NOTHREADS
        cerr << "Sorry, this is the no-threads executable and you requested " << arg.th << " threads\n";
        exit(EXIT_FAILURE);
        #else
        totChar = mt_process_file(arg,wordFreq);
        #endif
      }
  }
  catch(const std::bad_alloc&) {
      cout << "Out of memory (parsing phase)... emergency exit\n";
      die("bad alloc exception");
  }
  // first report 
  uint64_t totDWord = kh_size(wordFreq);
  cout << "Total input symbols: " << totChar << endl;
  cout << "Found " << totDWord << " distinct words" <<endl;
  cout << "Parsing took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  // check # distinct words
  if(totDWord>MAX_DISTINCT_WORDS) {
    cerr << "Emergency exit! The number of distinc words (" << totDWord << ")\n";
    cerr << "is larger than the current limit (" << MAX_DISTINCT_WORDS << ")\n";
    exit(1);
  }

  // -------------- second pass  
  start_wc = time(NULL);
  // create array of dictionary words
  vector<const char *> dictArray;
  dictArray.reserve(totDWord);
  // fill array
  uint64_t sumLen = 0;
  uint64_t totWord = 0;
  for(khiter_t ki = 0; ki < kh_end(wordFreq); ++ki) {
    sumLen += std::strlen(kh_val(wordFreq, ki).str);
    totWord += kh_val(wordFreq, ki).occ;
    dictArray.push_back(kh_val(wordFreq, ki).str);
  }
  assert(dictArray.size()==totDWord);
  cout << "Sum of lengths of dictionary words: " << sumLen << endl; 
  cout << "Total number of words: " << totWord << endl; 
  // sort dictionary
  sort(dictArray.begin(), dictArray.end(),[](const char *x, const char *y) {return std::strcmp(x, y) < 0;});
  // write plain dictionary and occ file, also compute rank for each hash 
  cout << "Writing plain dictionary and occ file\n";
  writeDictOcc(arg, wordFreq, dictArray);
  dictArray.clear(); // reclaim memory
  cout << "Dictionary construction took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
    
  // remap parse file
  start_wc = time(NULL);
  cout << "Generating remapped parse file\n";
  remapParse(arg, wordFreq);
  cout << "Remapping parse file took: " << difftime(time(NULL),start_wc) << " wall clock seconds\n";  
  cout << "==== Elapsed time: " << difftime(time(NULL),start_main) << " wall clock seconds\n";        
  for(khiter_t ki = 0; ki < kh_end(wordFreq); ++ki)
    if(kh_exist(wordFreq, ki))
      std::free(kh_val(wordFreq, ki).str);
  return 0;
}

