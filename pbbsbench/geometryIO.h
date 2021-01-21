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

//using namespace geometry;
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

  // string HeaderPoint2d = "pbbs_sequencePoint2d";
  // string HeaderPoint3d = "pbbs_sequencePoint3d";
  string HeaderTriangles = "pbbs_triangles";

  string pointHeader(int dim) {
    if (dim < 2 || dim > 9) {
      cout << "Error, unsupported dimension " << dim << ", abort." << endl;
      abort();
    }
    return "pbbs_sequencePoint" + to_string(dim) + "d";
  }

  inline int readDimensionFromFile(char* fname) {
    parlay::sequence<char> S = readStringFromFile(fname);
    parlay::sequence<char*> W = stringToWords(S);
    return (int)(W[0][18])-48;
  }

  template <class pointT>
  int writePointsToFile(parlay::sequence<pointT> const &P, char const *fname) {
    //string Header = (pointT::dim == 2) ? HeaderPoint2d : HeaderPoint3d;
    string Header = pointHeader(pointT::dim);
    int r = writeSeqToFile(Header, P, fname);
    return r;
  }

  template <class pointT, class Seq>
  parlay::sequence<pointT> parsePoints(Seq W) {
    //using coord = typename Point::coord;
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
    if (W.size() == 0 || W[0] != pointHeader(d)) {
      cout << "readPointsFromFile wrong file type" << endl;
      abort();
    }
    return parsePoints<pointT>(W.cut(1,W.size()));
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

