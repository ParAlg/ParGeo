// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#pragma once
#include "../parlay/parallel.h"
#include "../parlay/primitives.h"
#include "geometry.h"
#include "IO.h"

using namespace benchIO;

template <class coord>
inline int xToStringLen(point2d<coord> a) {
  return xToStringLen(a.x) + xToStringLen(a.y) + 1;
}

template <class coord>
inline void xToString(char* s, point2d<coord> a) {
  int l = xToStringLen(a.x);
  xToString(s, a.x);
  s[l] = ' ';
  xToString(s+l+1, a.y);
}

template <class coord>
inline int xToStringLen(point3d<coord> a) {
  return xToStringLen(a.x) + xToStringLen(a.y) + xToStringLen(a.z) + 2;
}

template <class coord>
inline void xToString(char* s, point3d<coord> a) {
  int lx = xToStringLen(a.x);
  int ly = xToStringLen(a.y);
  xToString(s, a.x);
  s[lx] = ' ';
  xToString(s+lx+1, a.y);
  s[lx+ly+1] = ' ';
  xToString(s+lx+ly+2, a.z);
}

namespace benchIO {
  using namespace std;

  string HeaderTriangles = "pbbs_triangles";

  string pbbsHeader(int dim) {
    if (dim < 2 || dim > 9) {
      cout << "Error, unsupported dimension " << dim << ", abort." << endl;
      abort();
    }
    return "pbbs_sequencePoint" + to_string(dim) + "d";
  }

  bool isGenericHeader(std::string line) {
    for (auto c: line) {
      if (!is_number(c) && !is_delim(c)) return true;
    }
    return false;
  }

  int countEntry(std::string line) {
    while (is_delim(line.back()) ||
	   is_space(line.back()) ||
	   is_newline(line.back())) {
      line.pop_back();
    }

    int count = 0;
    for (auto c: line) {
      if (is_delim(c)) count ++;
    }
    return count + 1;
  }

  // returns dim
  int readHeader(char* fileName) {
    ifstream file (fileName);
    if (!file.is_open())
      throw std::runtime_error("Unable to open file");

    std::string line1; std::getline(file, line1);
    if (isGenericHeader(line1)) {
      std::string line2; std::getline(file, line2);
      return countEntry(line2);
    } else {
      return countEntry(line1);
    }
  }

  // todo deprecate
  int readDimensionFromFile(char* fileName) {
    cout << "warning: using deprecated function readDimensionFromFile" << endl;
    return readHeader(fileName);
  }

  template <class pointT>
  int writePointsToFile(parlay::sequence<pointT> const &P, char const *fname) {
    string Header = pbbsHeader(pointT::dim);
    int r = writeSeqToFile(Header, P, fname);
    return r;
  }

  template <class pointT, class Seq>
  parlay::sequence<pointT> parsePoints(Seq W) {
    using coord = double;
    int d = pointT::dim;
    size_t n = W.size()/d;
    auto a = parlay::tabulate(d * n, [&] (size_t i) -> coord {
				       return atof(W[i]);});
    auto points = parlay::tabulate(n, [&] (size_t i) -> pointT {
					return pointT(a.cut(d*i,d*(i + 1)));});
    return points;
  }

  template <class pointT>
  parlay::sequence<pointT> readPointsFromFile(char const *fname) {
    parlay::sequence<char> S = readStringFromFile(fname);
    parlay::sequence<char*> W = stringToWords(S);
    int d = pointT::dim;
    if (W.size() == 0)
      throw std::runtime_error("readPointsFromFile empty file");

    if (isGenericHeader(W[0]))
      return parsePoints<pointT>(W.cut(1,W.size()));
    else
      return parsePoints<pointT>(W.cut(0,W.size()));
  }

  template <class pointT>
  triangles<pointT> readTrianglesFromFile(char const *fname, int offset) {
    int d = pointT::dim;
    parlay::sequence<char> S = readStringFromFile(fname);
    parlay::sequence<char*> W = stringToWords(S);
    if (W[0] != HeaderTriangles) {
      cout << "readTrianglesFromFile wrong file type" << endl;
      abort();
    }

    int headerSize = 3;
    size_t n = atol(W[1]);
    size_t m = atol(W[2]);
    if (W.size() != headerSize + 3 * m + d * n) {
      cout << "readTrianglesFromFile inconsistent length" << endl;
      abort();
    }

    auto pts_slice = W.cut(headerSize, headerSize + d * n);
    auto tri_slice = W.cut(headerSize + d * n, W.size());
    parlay::sequence<pointT> Pts = parsePoints<pointT>(pts_slice);
    auto Tri = parlay::tabulate(m, [&] (size_t i ) -> tri {
				     return {(int) atol(tri_slice[3*i])-offset,
					     (int) atol(tri_slice[3*i+1])-offset,
					     (int) atol(tri_slice[3*i+2])-offset};});
    return triangles<pointT>(Pts,Tri);
  }

  template <class pointT>
  int writeTrianglesToFile(triangles<pointT> Tr, char* fileName) {
    ofstream file (fileName, ios::binary);
    if (!file.is_open()) {
      std::cout << "Unable to open file: " << fileName << std::endl;
      return 1;
    }
    file << HeaderTriangles << endl;
    file << Tr.numPoints() << endl; 
    file << Tr.numTriangles() << endl; 
    writeSeqToStream(file, Tr.P);
    //writeSeqToStream(file, Tr.T);
    auto A = parlay::tabulate(3*Tr.numTriangles(), [&] (size_t i) -> int {
      						     return (Tr.T[i/3])[i%3];});
    writeSeqToStream(file, A);
    file.close();
    return 0;
  }

};

