#pragma once

#include "pbbs/utils.h"
#include "pbbs/parallel.h"
#include "pbbs/sequence.h"
#include "pbbs/gettime.h"

template<class T, class F, class G>
void serial_prefix(T* A, intT n, F& process, G& cleanUp) {
  for(intT i=0; i<n; ++i) {
    if (process(A[i])) {
      cleanUp(A, i);}
  }
}

template<class T, class F, class G>
void parallel_prefix(T* A, intT n, F& process, G& cleanUp, bool verbose=false, intT* flag=NULL, intT PREFIX=2000, intT THRESH=5000) {
  if (n < THRESH) return serial_prefix(A, n, process, cleanUp);

  intT prefix;
  intT i = 0;
  intT conflict;

  bool freeFlag = false;
  if (!flag) {
    flag = newA(intT, n+1);
    freeFlag = true;}

  /* timing t0; */
  /* intT serCount = 0; */
  /* intT serSize = 0; */
  /* floatT serTime = 0; */
  /* intT parCount = 0; */
  /* intT parSize = 0; */
  /* floatT parTime = 0; */
  /* floatT parTime1 = 0; */
  /* floatT parTime2 = 0; */
  /* floatT parTime3 = 0; */

  while (i < n) {
    if (i==0) prefix = PREFIX;//prefix < n is known
    else prefix = min(i*3, n);//prefix to be taken as 3*i

    //if (verbose) cout << "prefix = " << prefix << endl;

    if (prefix-i < THRESH) {
      //if(verbose) t0.start();
      /* serCount += 1; */
      /* serSize += prefix - i; */
      while (i < prefix) {
	if (process(A[i])) cleanUp(A, i);
	i ++;
      }
      //if(verbose) serTime += t0.stop();
    } else {
      //if(verbose) t0.start();
      parCount += 1;
      parSize += prefix - i;

      auto body = [&](size_t j) {
		    if (process(A[j])) flag[j-i] = 1;//conflict
		    else flag[j-i] = 0;
		  };
      parallel_for(i, prefix, body);
      //if(verbose) parTime += t0.next();

      intT numBad = sequence::prefixSum(flag, 0, prefix-i);
      flag[prefix-i] = numBad;
      //if(verbose) parTime1 += t0.stop();

      if (numBad > 0) {
	//if(verbose) t0.start();
	auto body = [&](size_t j) {
		      if (flag[j]==0 && flag[j]!=flag[j+1]) conflict = j;
		    };
	parallel_for(0, prefix-i, body);
	i += conflict;
	//if(verbose) parTime2 += t0.next();
	cleanUp(A, i);
	//if(verbose) parTime3 += t0.stop();
      } else {
	i = prefix;
      }

    }
  }

/* #ifndef SILENT */
/*   if(verbose) { */
/*     cout << "serial prefix stats" << endl; */
/*     cout << " count = " << serCount << endl; */
/*     cout << " avg size = " << serSize/(floatT)serCount << endl; */
/*     cout << " total time = " << serTime << endl; */
/*     cout << "parallel prefix stats" << endl; */
/*     cout << " count = " << parCount << endl; */
/*     cout << " avg size = " << parSize/(floatT)parCount << endl; */
/*     cout << " total time 1 = " << parTime << endl; */
/*     cout << " total time 2 = " << parTime1 << endl; */
/*     cout << " total time 3 = " << parTime2 << endl; */
/*     cout << " total time 4 = " << parTime3 << endl; */
/*   } */
/* #endif */

  if(freeFlag) free(flag);
}

#endif
