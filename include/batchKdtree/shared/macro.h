// This code is part of the project "Parallel Batch-Dynamic Kd-Trees"
// Copyright (c) 2021-2022 Rahul Yesantharao, Yiqiu Wang, Laxman Dhulipala, Julian Shun
//
// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

namespace pargeo::batchKdTree {

//#define PRINT_CONFIG
//#define PRINT_LOGTREE_TIMINGS
//#define PRINT_DKNN_TIMINGS
//#define PRINT_INSERT_TIMINGS
//#define PRINT_DELETE_TIMINGS
//#define PRINT_KDTREE_TIMINGS
//#define PRINT_COKDTREE_TIMINGS
//#define USE_MEDIAN_SELECTION
//#define PRINT_PARALLEL_PARTITION_TIMINGS
//#define ERASE_SEARCH_TIMES

//#define ALL_USE_BLOOM
//#define LOGTREE_USE_BLOOM
//#define BLOOM_FILTER_BUILD_COPY

#define PARTITION_OBJECT_MEDIAN 0
#define PARTITION_SPATIAL_MEDIAN 1
#ifndef PARTITION_TYPE
#define PARTITION_TYPE PARTITION_OBJECT_MEDIAN
#endif

// DUAL KNN MODE
#define DKNN_ATOMIC_LEAF 0
#define DKNN_NONATOMIC_LEAF 1
#define DKNN_ARRAY 2

#ifndef DUAL_KNN_MODE  // default if not defined in cmake
#define DUAL_KNN_MODE DKNN_NONATOMIC_LEAF
#endif

// LOGTREE BUFFER
#define BHL_BUFFER 0
#define ARR_BUFFER 1
#define LOGTREE_BUFFER BHL_BUFFER

// KNN OPTIMIZATION
#define FEWER_SQRT 1

// LEAF CLUSTER SIZE
#ifndef CLUSTER_SIZE
#define CLUSTER_SIZE 16
#endif

// SERIAL BASE CASES
#ifndef ERASE_BASE_CASE
#define ERASE_BASE_CASE 1000
#endif

#ifndef RANGEQUERY_BASE_CASE
#define RANGEQUERY_BASE_CASE 1000
#endif

#ifndef BOUNDINGBOX_BASE_CASE
#define BOUNDINGBOX_BASE_CASE 1000
#endif

#ifndef DUALKNN_BASE_CASE
#define DUALKNN_BASE_CASE 1000
#endif

#ifndef CO_TOP_BUILD_BASE_CASE
#define CO_TOP_BUILD_BASE_CASE 1000
#endif

#ifndef CO_BOTTOM_BUILD_BASE_CASE
#define CO_BOTTOM_BUILD_BASE_CASE 1000
#endif

#ifndef BHL_BUILD_BASE_CASE
#define BHL_BUILD_BASE_CASE 1000
#endif

#ifdef PRINT_CONFIG
#include <iostream>
void print_config() {
  std::cout << "DUAL_KNN_MODE = " << DUAL_KNN_MODE << ";\n"
            << "PARTITION_TYPE = " << PARTITION_TYPE << ";\n"
            << "LOGTREE_BUFFER = " << LOGTREE_BUFFER << ";\n"
            << "CLUSTER_SIZE = " << CLUSTER_SIZE << ";\n"
            << "ERASE_BASE_CASE = " << ERASE_BASE_CASE << ";\n"
            << "RANGEQUERY_BASE_CASE = " << RANGEQUERY_BASE_CASE << ";\n"
            << "BOUNDINGBOX_BASE_CASE = " << BOUNDINGBOX_BASE_CASE << ";\n"
            << "DUALKNN_BASE_CASE = " << DUALKNN_BASE_CASE << ";\n"
            << "CO_TOP_BUILD_BASE_CASE = " << CO_TOP_BUILD_BASE_CASE << ";\n"
            << "CO_BOTTOM_BUILD_BASE_CASE = " << CO_BOTTOM_BUILD_BASE_CASE << ";\n"
            << "BHL_BUILD_BASE_CASE = " << BHL_BUILD_BASE_CASE << std::endl;
}
#else
void print_config() {}
#endif

} // End namespace batchKdTree