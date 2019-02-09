//#include <google/sparse_hash_map>
//#include <google/dense_hash_map>
//#include "MurmurHash3.cpp"
#include <iostream>
#include <climits>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <time.h>
#include "metrohash64.cpp"
#include <stdint.h>
#include <unordered_map>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <string.h>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <stack>
#include <limits.h>
#include <map>
#include <bitset>
#include <ctime>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include <random>
#include <cinttypes>
//#include "dna_test.h"
#include "ntHashIterator.hpp"
#include "FastxParser.hpp"
#include <thread>
#include <mutex>
#include <chrono>
#define _POSIX_C_SOURCE 200809L


#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"

using spp::sparse_hash_map;

using namespace std;
//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)
std::map<char, char> mapp = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}};

static const int MultiplyDeBruijnBitPosition[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};
unsigned trailing_zeros(unsigned n) {
    return n ? __builtin_ctz(n) : -1;
}

static const char basemap[255] =
    {
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };

unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}

void printHelp()
{

    cout << "KmerEst [options] -f <fasta/fastq> -k <k-mer length>  -s <sample size> -o <output file>"    << endl
    << "  -h               help"                                   << endl
    << "  -f <file>       Input sequence file "                << endl
    << "  -k <k-mer size >        kmer size (default 31) "        << endl
    << "  -s <sample size>        sample size (default 25m)"        << endl
     << "  -c coverage>       coverage (default 64)"        << endl
    << "  -o         	  Prefix of the Output file " << endl;

    exit(0);
}

int main(int argc, char** argv)
{

    if(argc == 1){
      cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -o <out.txt>" << endl;
      exit(0);
    }
    int n=31;
    int s = 25000000;
    int cov = 64;
    string f = "", outf = "";
    for (int c = 1; c < argc; c++)
        {

            if (!strcmp(argv[c], "-h"))       { printHelp(); }
            else if (!strcmp(argv[c], "-k"))     { n = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-f"))    { f = argv[c+1]; c++; }
            else if (!strcmp(argv[c], "-s"))    { s = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-c"))    { cov = atoi(argv[c+1]); c++; }
            else if (!strcmp(argv[c], "-o")) { outf = argv[c+1]; c++; }
        }

       if (f.empty()  || outf.empty())
        {
          printHelp();
        }

    int k = s;
    mutex locks[64],count_lock,th_lock;
    typedef sparse_hash_map<uint64_t, uint32_t> SMap;
    vector<SMap> MAP(64);


    // nt is the no of reader threads on FQFeeder queue, np is the no of threads writing to the queue
    size_t nt = 4;
    size_t np = 1;

    int th = 0;

    // Declared total_t, no_kmers_t for counting total no of reads, kmers handled by different threads
    uint64_t total_t[nt] = {0}, no_kmers_t[nt] = {0};
    int count = 0;
    uint64_t hash[nt]={0}, fhVal=0, rhVal=0;
    std::vector<std::string> fp;
    fp.push_back(f);

    //Creating instance of FastxParser
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(fp, nt, np);

    parser.start();
    cout << "read the Sequences .. " << endl;

    std::vector<std::thread> readers;
	using namespace std::chrono;
    high_resolution_clock::time_point start = chrono::high_resolution_clock::now();

    for (size_t i = 0; i < nt; ++i) {
    readers.emplace_back([&, i]() {
      auto rg = parser.getReadGroup();
      while (true) {
        if (parser.refill(rg)) {
          for (auto& seq : rg) {
              total_t[i]++;
              ntHashIterator itr(seq.seq, 1, n);
              while (itr != itr.end()) {
                hash[i] = (*itr)[0];
                no_kmers_t[i]++;
                uint8_t tz = trailing_zeros(hash[i]);
                if(tz >= th){
                    // A thread is acquiring lock on the hash table which has k-mers with no of trailing zeros tz
                    locks[tz].lock();
                    if(MAP[tz].find(hash[i]) != MAP[tz].end()) MAP[tz][hash[i]] += 1;
                    else{
                      MAP[tz].insert(make_pair(hash[i], 1));
                      count_lock.lock();
                      ++count;
                      if(count == k){
                        int cnt = MAP[th].size();
                        count = count - cnt;
                        SMap().swap(MAP[th]); // The internal logic of this swap is that it deletes hashmap with number of trailing zeros th
                        th_lock.lock();
                        ++th;
                        th_lock.unlock();

                        cout  << "count: " << count << endl;
                      }
                      count_lock.unlock();
                    }
                    locks[tz].unlock();
                }
                ++itr;
              }

          }
        } else {
          break;
        }
      }
    }
    );
  }

    for (auto& t : readers) {
      t.join();
    }
    parser.stop();
    // Combining values of no_kmers and total no of sequences from all the threads
    uint64_t total = 0,no_kmers=0;
    for(int j=0;j<nt;j++){
       total+=total_t[j];
       no_kmers+=no_kmers_t[j];
    }
    high_resolution_clock::time_point stop = chrono::high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(stop - start);
    cout<< "execution time: "<<time_span.count()<<endl;


    cout << "th: " << th << endl;
    cout << "No. of sequences: " << total << endl;
    FILE *fo = fopen(outf.c_str(), "w");
    uint32_t csize = 0;
    for(int i=th; i<64; i++) csize += MAP[i].size();
    unsigned long F0 = csize * pow(2, (th));
    cout << "F0: " << F0 << endl;
    fprintf(fo, "F1\t%lu\n", no_kmers);
    fprintf(fo, "F0\t%lu\n", F0);
    cout << endl;
    cout << "total: " << total << endl;
    cout << "no_kmer: " << no_kmers << endl;
   unsigned long *freq = new unsigned long[cov];
   for(int i=1; i<=cov; i++) freq[i] = 0;
    unsigned long tot = 0;
    int xx = 0;
    for(int i=th; i<64; i++){
      for(auto& p: MAP[i]){
        if(p.second <= cov) freq[p.second]++;
      }
    }

    cout << "th: " << th << endl;
    for(int i=1; i<=cov; i++){
      unsigned long fff = (freq[i]*pow(2, th));
      fprintf(fo, "f%d\t%lu\n", i, fff);
    }
    fclose(fo);
    return 0;

}