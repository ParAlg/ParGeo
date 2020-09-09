#ifndef PARALLEL_H
#define PARALLEL_H

#include <limits>
#include <iostream>
using namespace std;

typedef int intT;
typedef unsigned int uintT;
typedef double floatT;
static intT intMax() {return numeric_limits<intT>::max();}
static uintT uintMax() {return numeric_limits<uintT>::max();}
static floatT floatMax() {return numeric_limits<floatT>::max();}

#if defined(OPENCILK)

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#define parallel_main main
#define par_for cilk_for
#define par_for_1 _Pragma("cilk_grainsize = 1") cilk_for
#define par_for_256 _Pragma("cilk_grainsize = 256") cilk_for

extern "C" int __cilkrts_internal_worker_id(void);

static int getWorkers() {
  return __cilkrts_get_nworkers();}
static int getWorkerId() {return __cilkrts_internal_worker_id();}
static void setWorkers(int n) { }
static void printScheduler() {
  cout << "scheduler = OpenCilk" << endl;
  cout << "num-threads = " << getWorkers() << endl;
}

#elif defined(CILK)

#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#define parallel_main main
#define par_for cilk_for
#define par_for_1 _Pragma("cilk_grainsize = 1") cilk_for
#define par_for_256 _Pragma("cilk_grainsize = 256") cilk_for

static int getWorkers() {
  return __cilkrts_get_nworkers();}
static int getWorkerId() {return __cilkrts_get_worker_number();}
static void setWorkers(int n) { }
static void printScheduler() {
  cout << "scheduler = CilkPlus" << endl;
  cout << "num-threads = " << getWorkers() << endl;
}

#else

#define cilk_spawn
#define cilk_sync
#define parallel_main main
#define par_for for
#define par_for_1 for
#define par_for_256 for

static void printScheduler() {
  cout << "scheduler = sequential" << endl;}
static int getWorkers() {return 1;}
static int getWorkerId() {return 0;}
static void setWorkers(int n) { }

#endif

#endif
