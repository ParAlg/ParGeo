#pragma once

//#include <vector>
#include "pargeo/point.h"
#include "parlay/sequence.h"

/* #include "pbbs/gettime.h" */
/* #include "pbbs/utils.h" */
/* #include "pbbs/randPerm.h" */
// #include "prefix.h"
// #include "miniDisc.h"
// #include "geometry.h"

// using namespace std;

// static size_t counter = 0;

template<int dim>
pargeo::seb::ball<dim> support2Ball(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
				    parlay::sequence<pargeo::point<dim>>& support) {
  using ballT = pargeo::seb::ball<dim>;
  // for (auto s: support) std::cout << s << ", ";
  // std::cout << "\n";

  ballT B;
  if (B.isEmpty()) {
    if (support.size() == 0) {
      //std::cout << "support2Ball-1\n";
      //B = ballT(P, 2);
      B = ballT(P.cut(0, 2));
    } else if (support.size() == 1) {
      //std::cout << "support2Ball-2\n";
      support.push_back(P[0]);
      //B = ballT(&support[0], support.size());
      B = ballT(parlay::make_slice(support));
      support.pop_back();
    } else { // >=2
      //std::cout << "support2Ball-3\n";
      B = ballT(parlay::make_slice(support));
    }
  }
  //std::cout << "ball size = " << B.size() << "\n";
  return B;
}

template<int dim>
pargeo::seb::ball<dim> miniDiscPlainSerial(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
					   parlay::sequence<pargeo::point<dim>>& support,
					   pargeo::seb::ball<dim> B) {
  typedef pargeo::seb::ball<dim> ballT;
  typedef pargeo::point<dim> pointT;

  // if (counter++ > 10000) abort();
  // std::cout << "\ncall " << counter << ": " << P.size() << ", " << B.size() << "\n";

  B = support2Ball(P, support);

  if (B.size() == dim+1) return B;

  for (size_t i=0; i<P.size(); ++i) {
    //std::cout << "\nprocessing i " << i << "\n\n";

    if (!B.contain(P[i])) { // process a conflict, ow keep going
      //std::cout << " -> not contain\n";

      if (support.size() == B.size()) B.grow(P[i]);
      else B = ballT();

      support.push_back(P[i]);

      // compute a ball on the prefix
      // B = miniDiscPlainSerial(P, i, support, B);
      B = miniDiscPlainSerial(P.cut(0, i), support, B);

      support.pop_back();
    }

  }

  return B;
}

template<int dim>
pargeo::seb::ball<dim> miniDiscPlain(parlay::slice<pargeo::point<dim>*, pargeo::point<dim>*> P,
				     parlay::sequence<pargeo::point<dim>>& support,
				     pargeo::seb::ball<dim> B, size_t* flag=NULL) { // todo flag can be removed
  return miniDiscPlainSerial<dim>(P, support, B);
  // todo below

  /* typedef ball<dim> ballT; */
  /* typedef point<dim> pointT; */

  /* if (n < 2000) return miniDiscPlainSerial(P, n, support, B); */

  /* B = support2Ball(P, support); */
  /* if (B.size() == dim+1) return B; */

  /* bool freeFlag = false; */
  /* if (!flag) { */
  /*   freeFlag = true; */
  /*   flag = newA(intT, n+1);} */

  /* auto process = [&](pointT p) { */
  /*                  if (!B.contain(p)) return true; */
  /*                  else return false; */
  /*                }; */
  /* auto cleanUp = [&](pointT* A, intT i) { */
  /*                  if (support.size() == B.size()) B.grow(A[i]); */
  /*                  else B = ballT(); */
  /*                  support.push_back(A[i]); */
  /*                  B = miniDiscPlain(A, i, support, B, flag); */
  /*                  support.pop_back(); */
  /*                }; */
  /* parallel_prefix(P, n, process, cleanUp, freeFlag, flag, 200000, 500000); */

  /* if(freeFlag) free(flag); */
  /* return B; */
}
