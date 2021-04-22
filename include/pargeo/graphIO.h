#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <string>
#include <cstring>
#include "../parlay/primitives.h"
#include "../parlay/parallel.h"
#include "../parlay/io.h"

namespace pargeo {
namespace graphIO {
  using namespace std;
  using parlay::sequence;
  using parlay::tabulate;
  using parlay::make_slice;

  ///////////////////////////////////////////////////////////////////
  // Write sequence of edges (u,v) to an edge list (SNAP format)
  ///////////////////////////////////////////////////////////////////

  inline int xToStringLen(unsigned long a) { return 21;}
  inline void xToString(char* s, unsigned long a) { sprintf(s,"%lu",a);}

  string emptyHeaderIO = "";

  template <class Seq>
    parlay::sequence<char> edgeSeqToString(Seq const &edges) {
    parlay::sequence<size_t> A(edges.size()*2);
    parlay::parallel_for(0, edges.size(), [&](size_t i){
	A[i*2] = edges[i].u;
	A[i*2+1] = edges[i].v;
      });

    size_t n = A.size();
    auto L = parlay::tabulate(n, [&] (size_t i) -> long {
	size_t x = A[i];
	return xToStringLen(x)+1;});
    size_t m;
    std::tie(L,m) = parlay::scan(std::move(L));

    parlay::sequence<char> B(m+1, (char) 0);
    char* Bs = B.begin();

    parlay::parallel_for(0, n-1, [&] (long i) {
	xToString(Bs + L[i], A[i]);
	if (i % 2 == 1)
	  Bs[L[i+1] - 1] = '\n';
	else
	  Bs[L[i+1] - 1] = ' ';
      });
    xToString(Bs + L[n-1], A[n-1]);
    Bs[m] = Bs[m-1] = '\n';

    parlay::sequence<char> C = parlay::filter(B, [&] (char c) {return c != 0;});
    C[C.size()-1] = 0;
    return C;
  }

  template <class T>
    void writeEdgeSeqToStream(ofstream& os, parlay::sequence<T> const &A) {
    size_t bsize = 10000000;
    size_t offset = 0;
    size_t n = A.size();
    while (offset < n) {
      // Generates a string for a sequence of size at most bsize
      // and then wrties it to the output stream
      parlay::sequence<char> S = edgeSeqToString(A.cut(offset, min(offset + bsize, n)));
      os.write(S.begin(), S.size()-1);
      offset += bsize;
    }
  }

  template <class T>
    int writeEdgeSeqToFile(string header,
			   parlay::sequence<T> const &A,
			   char const *fileName) {
    auto a = A[0];

    ofstream file (fileName, ios::out | ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }

    //file << header << endl; // not writing a header

    writeEdgeSeqToStream(file, A);
    file.close();
    return 0;
  }

  template <class T>
    int writeEdgeSeqToFile(parlay::sequence<T> const &A, char const *fileName) {
    return writeEdgeSeqToFile(emptyHeaderIO, A, fileName);
  }

} // End namespace graphIO
} // End namespace pargeo
