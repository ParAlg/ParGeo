// Copyright (c) 2020 Yiqiu Wang and the Pargeo Team
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

#include <vector>
#include "pbbs/gettime.h"
#include "pbbs/utils.h"
#include "pbbs/randPerm.h"
#include "prefix.h"
#include "miniDisc.h"
#include "geometry.h"
#include "check.h"
#include "welzl.h"

using namespace std;

static const floatT numericKnob = 1e-6;

template<int dim>
intT ortScanSerial(point<dim> c, floatT rSqr, point<dim>* P, intT n, vector<point<dim>>& support, floatT* dist) {
  typedef point<dim> pointT;

  intT dd = intT(pow(2.0, dim));
  intT idx[dd];
  for(intT i=0; i<dd; ++i) idx[i] = -1;

  intT confCount = 0;
  for (intT i=0; i<n; ++i) {
    floatT dSqr = P[i].distSqr(c);
    if (dSqr > rSqr*(1+numericKnob)) {
      confCount ++;
      intT o = c.quadrant(P[i]);
      if (dSqr > dist[o]*(1+numericKnob)) {
        dist[o] = dSqr;
        idx[o] = i;}
    }
  }
  if (confCount > 0)
    cout << "conflicts = " << confCount << "/" << n << " (" << 100*floatT(confCount)/n << " %)" << endl;

  for(intT i=0; i<dd; ++i) {
    if (idx[i] != -1) {
      support.push_back(P[idx[i]]);
    }
  }
  return confCount;
}

template<int dim>
intT ortScan(point<dim> c, floatT rSqr, point<dim>* A, intT n, vector<point<dim>>& support, floatT* distGlobal) {
  typedef point<dim> pointT;

  if (n<2000) return ortScanSerial(c, rSqr, A, n, support, distGlobal);

  intT dd = intT(pow(2.0, dim));
  intT P = getWorkers()*8;//todo tune
  intT blockSize = (n+P-1)/P;
  intT idx[dd*P];
  floatT dist[dd*P];
  intT confs[P];
  for (intT i=0; i<dd*P; ++i) idx[i] = -1;
  for (intT i=0; i<P; ++i) {
    for (intT j=0; j<dd; ++j) {
      dist[i*dd+j] = distGlobal[j];}
  }

  parallel_for(0, P,
	       [&](intT p) {
		 intT s = p*blockSize;
		 intT e = min((intT)(p+1)*blockSize,n);
		 intT* locIdx = idx + p*dd;
		 floatT* locDist = dist + p*dd;
		 confs[p] = 0;
		 for (intT i=s; i<e; ++i) {
		   floatT dSqr = A[i].distSqr(c);
		   if (dSqr > rSqr*(1+numericKnob)) {
		     confs[p]++;
		     intT o = c.quadrant(A[i]);
		     if (dSqr > locDist[o]*(1+numericKnob)) {
		       locDist[o] = dSqr;
		       locIdx[o] = i;}
		   }
		 }
	       }, 1);

  intT idxGlobal[dd];
  for (intT o=0; o<dd; ++o) idxGlobal[o] = -1;

  intT totalConfs = 0;
  for(intT p=0; p<P; ++p) {
    totalConfs += confs[p];
    for(intT o=0; o<dd; ++o) {
      intT* locIdx = idx + p*dd;
      floatT* locDist = dist + p*dd;
      if (locIdx[o] != -1) {
        if(locDist[o] > distGlobal[o]) {
          idxGlobal[o] = locIdx[o];
          distGlobal[o] = locDist[o];}
      }
    }
  }

  if (totalConfs > 0)
    cout << "conflicts = " << totalConfs << "/" << n << " (" << 100*floatT(totalConfs)/n << " %)" << endl;

  for(intT o=0; o<dd; ++o) {
    if (idxGlobal[o] != -1) {
      support.push_back(A[idxGlobal[o]]);}
  }
  return totalConfs;
}

template<int dim>
ball<dim> miniDiscSamplingOrtSerial(point<dim>* P, intT n) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  intT sample = dim*3;
  ballT B;
  if (sample > n) {
    vector<pointT> support;
    return miniDiscPlainSerial(P, sample, support, B);
  } else {
    vector<pointT> support;
    B = miniDiscPlainSerial(P, sample, support, B);
  }

  intT dd = intT(pow(2.0, dim));
  floatT dist[dd];
  for(intT i=0; i<dd; ++i) dist[i] = -1;

  floatT sampleTime = 0;
  floatT remainTime = 0;
  timing t; t.start();

  intT scanned = 0;
  intT step = n/100;
  while (scanned < n) {
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScanSerial<dim>(B.center(), B.radius()*B.radius(), P+scanned, step, support, dist);
    scanned += step;

    if (!found) {
      break;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlainSerial(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  sampleTime = t.next();

  cout << "---" << endl;

  while (1) {
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScanSerial<dim>(B.center(), B.radius()*B.radius(), P, n, support, dist);

    if (!found) {
      remainTime = t.stop();
      cout << "sample-time = " << 100*sampleTime/(sampleTime+remainTime) << " %" << endl;
      cout << "remain-time = " << 100*remainTime/(sampleTime+remainTime) << " %" << endl;

      return B;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlainSerial(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  remainTime = t.stop();
  cout << "sample-time = " << 100*sampleTime/(sampleTime+remainTime) << " %" << endl;
  cout << "remain-time = " << 100*remainTime/(sampleTime+remainTime) << " %" << endl;

  return B;
}

template<int dim>
ball<dim> miniDiscSamplingOrt(point<dim>* P, intT n) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;
  if (n < 2000)  return miniDiscSamplingOrtSerial(P, n);

  intT sample = dim*3;
  ballT B;
  if (sample > n) {
    vector<pointT> support;
    return miniDiscPlain(P, sample, support, B);
  } else {
    vector<pointT> support;
    B = miniDiscPlain(P, sample, support, B);
  }

  intT dd = intT(pow(2.0, dim));
  floatT dist[dd];
  for(intT i=0; i<dd; ++i) dist[i] = -1;

  floatT sampleTime = 0;
  floatT remainTime = 0;
  timing t; t.start();

  intT scanned = 0;
  intT step = n/100;
  while (scanned < n) {
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScan<dim>(B.center(), B.radius()*B.radius(), P+scanned, step, support, dist);
    scanned += step;

    if (!found) {
      break;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlain(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  sampleTime = t.next();

  cout << "---" << endl;

  while (1) {
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScan<dim>(B.center(), B.radius()*B.radius(), P, n, support, dist);

    if (!found) {
      remainTime = t.stop();
      cout << "sample-time = " << 100*sampleTime/(sampleTime+remainTime) << " %" << endl;
      cout << "remain-time = " << 100*remainTime/(sampleTime+remainTime) << " %" << endl;

      return B;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlain(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  remainTime = t.stop();
  cout << "sample-time = " << 100*sampleTime/(sampleTime+remainTime) << " %" << endl;
  cout << "remain-time = " << 100*remainTime/(sampleTime+remainTime) << " %" << endl;

  return B;
}

template<int dim>
ball<dim> miniDiscSamplingOrt2(point<dim>* P, intT n) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;
  if (n < 2000)  return miniDiscSamplingOrtSerial(P, n);

  intT sample = dim*3;
  ballT B;
  if (sample > n) {
    vector<pointT> support;
    return miniDiscPlain(P, sample, support, B);
  } else {
    vector<pointT> support;
    B = miniDiscPlain(P, sample, support, B);
  }

  intT dd = intT(pow(2.0, dim));
  floatT dist[dd];
  for(intT i=0; i<dd; ++i) dist[i] = -1;

  floatT sampleTime = 0;
  floatT remainTime = 0;
  timing t; t.start();

  intT scanned = 0;
  intT step = 10;//find 10 conflicts each sample
  while (scanned < n/10) {//scan 10% of points
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    intT found = 0;
    intT scanStep = 100;
    while(found < step && scanned < n) {
      intT nn = n-scanned;
      found += ortScan<dim>(B.center(), B.radius()*B.radius(), P+scanned, min(nn, scanStep), support, dist);
      scanned += scanStep;
      scanStep *= 2;
    }
    cout << "found = " << found << "/" << scanned << endl;

    if (found == 0 || scanned >= n) {
      break;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlain(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  sampleTime = t.next();

  cout << "---" << endl;

  while (1) {
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScan<dim>(B.center(), B.radius()*B.radius(), P, n, support, dist);

    if (!found) {
      remainTime = t.stop();
      cout << "sample-time = " << 100*sampleTime/(sampleTime+remainTime) << " %" << endl;
      cout << "remain-time = " << 100*remainTime/(sampleTime+remainTime) << " %" << endl;

      return B;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlain(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  remainTime = t.stop();
  cout << "sample-time = " << 100*sampleTime/(sampleTime+remainTime) << " %" << endl;
  cout << "remain-time = " << 100*remainTime/(sampleTime+remainTime) << " %" << endl;

  return B;
}

template<int dim>
ball<dim> miniDiscOrt(point<dim>* P, intT n) {
  typedef ball<dim> ballT;
  typedef point<dim> pointT;

  if (n<2000) return miniDiscOrtSerial(P, n);

  intT sample = dim*3;
  ballT B;
  if (sample > n) {
    vector<pointT> support;
    return miniDiscPlain(P, sample, support, B);
  } else {
    vector<pointT> support;
    B = miniDiscPlain(P, sample, support, B);
  }

  intT dd = intT(pow(2.0, dim));
  floatT dist[dd];
  for(intT i=0; i<dd; ++i) dist[i] = -1;

  while (1) {
    vector<pointT> support;
    for(intT i=0; i<B.size(); ++i) {
      support.push_back(B.support()[i]);}

    bool found = ortScan<dim>(B.center(), B.radius()*B.radius(), P, n, support, dist);

    if (!found) {
      return B;
    } else {
      auto supportNew = vector<pointT>();
      B = miniDiscPlain(&support[0], support.size(), supportNew, ballT());
      for(intT i=0; i<dd; ++i) dist[i] = -1;
    }
  }
  return B;
}

template<int dim>
void miniDisc(point<dim>* P, intT n) {
  typedef point<dim> pointT;
  typedef circle discT;
  typedef ball<dim> ballT;

  static const bool preprocess = false;
  static const bool serial = false;

  cout << "smallest enclosing disc, " << n << ", dim " << dim << " points" << endl;

  timing t0;t0.start();
  if(preprocess) {
    randPerm(P, n);
    cout << "preprocess-time = " << t0.next() << endl;
  }

  ballT D;
  cout << "method = orthant-scan-sampling" << endl;
  if (serial) {
    D = miniDiscSamplingOrtSerial(P, n);
  } else {
    D = miniDiscSamplingOrt(P, n);
    //D = miniDiscSamplingOrt2(P, n);
  }

  cout << "seb-time = " << t0.stop() << endl;

  cout << D.radius() << ", center = " << D.center() << endl;

  cout << endl;

  check<dim,ballT>(&D, P, n);
}

template void miniDisc<2>(point<2>*, intT);
template void miniDisc<3>(point<3>*, intT);
template void miniDisc<4>(point<4>*, intT);
template void miniDisc<5>(point<5>*, intT);
template void miniDisc<6>(point<6>*, intT);
template void miniDisc<7>(point<7>*, intT);
template void miniDisc<8>(point<8>*, intT);
template void miniDisc<9>(point<9>*, intT);
