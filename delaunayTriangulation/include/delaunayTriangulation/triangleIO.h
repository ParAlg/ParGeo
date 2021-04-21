#pragma once

#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "pargeo/IO.h"

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

string HeaderTriangles = "pbbs_triangles";

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
