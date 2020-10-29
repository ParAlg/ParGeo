#ifndef MINI_DISC_2D_PREFIX_H
#define MINI_DISC_2D_PREFIX_H

#include "pbbs/utils.h"
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/gettime.h"
#include "geometry.h"
#include "miniDisc2d.h"

circle miniDisc2DPrefix(point<2>* A, intT n, bool verbose=true, intT* flag=NULL) {
  typedef circle circleT;

  static const intT PREFIX = 2000000;
  static const intT THRESH = 5000000;
  if (n < THRESH) return miniDisc2DSerial(A, n);

  auto circle = circleT(A[0], A[1]);

  intT prefix;
  intT i = 0;
  intT conflict;

  bool freeFlag = false;
  if (!flag) {
    flag = newA(intT, n+1);
    freeFlag = true;}

  timing t0;
  intT serCount = 0;
  intT serSize = 0;
  floatT serTime = 0;
  intT parCount = 0;
  intT parSize = 0;
  floatT parTime = 0;
  floatT parTime1 = 0;
  floatT parTime2 = 0;
  floatT parTime3 = 0;

  while (i < n) {
    if (i==0) prefix = PREFIX;//prefix < n is known
    else prefix = min(i*3, n);//prefix to be taken as 3*i

    //if (verbose) cout << "prefix = " << prefix << endl;

    if (prefix-i < THRESH) {
      if(verbose) t0.start();
      serCount += 1;
      serSize += prefix - i;
      while (i < prefix) {
        if (!circle.contain(A[i])) {
          circle = miniDisc2DSerial(A, i, A[i]);
          swap(A[2], A[i]);
        }
	i ++;
      }
      if(verbose) serTime += t0.stop();
    } else {
      if(verbose) t0.start();
      parCount += 1;
      parSize += prefix - i;

      par_for (intT j=i; j<prefix; ++j) {
	if (!circle.contain(A[j])) flag[j-i] = 1;//conflict
	else flag[j-i] = 0;
      }
      if(verbose) parTime += t0.next();

      intT numBad = sequence::prefixSum(flag, 0, prefix-i);
      flag[prefix-i] = numBad;
      if(verbose) parTime1 += t0.stop();

      if (numBad > 0) {
	if(verbose) t0.start();
	par_for(intT j=0; j<prefix-i; ++j) {
	  if (flag[j]==0 && flag[j]!=flag[j+1]) conflict = j;}
	i += conflict;
	if(verbose) parTime2 += t0.next();
	//if (verbose) cout << "ci = " << i << endl;

        circle = miniDisc2DSerial(A, i, A[i]);
        if(verbose) parTime3 += t0.stop();
        swap(A[2], A[i]);

      } else {
	i = prefix;
      }

    }
  }

  if(verbose) {
    cout << "serial prefix stats" << endl;
    cout << " count = " << serCount << endl;
    cout << " avg size = " << serSize/(floatT)serCount << endl;
    cout << " total time = " << serTime << endl;
    cout << "parallel prefix stats" << endl;
    cout << " count = " << parCount << endl;
    cout << " avg size = " << parSize/(floatT)parCount << endl;
    cout << " total time 1 = " << parTime << endl;
    cout << " total time 2 = " << parTime1 << endl;
    cout << " total time 3 = " << parTime2 << endl;
    cout << " total time 4 = " << parTime3 << endl;
  }

  if(freeFlag) free(flag);
  return circle;
}

#endif
