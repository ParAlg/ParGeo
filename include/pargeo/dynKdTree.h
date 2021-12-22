#pragma once

#include <algorithm>
#include <math.h>
#include <queue>
#include <iostream>
#include <vector>

#include "../parlay/parallel.h"
#include "../parlay/sequence.h"

template <typename T>
using container = parlay::sequence<T>;

namespace pargeo {

namespace dynKdTree {

  static const bool spatialMedian = true;

  template <typename In_Seq, typename Bool_Seq>
  size_t split_two(In_Seq const &In, Bool_Seq const &Fl, parlay::flags fl = parlay::no_flag) {
    using namespace parlay;
    using namespace parlay::internal;
    using T = typename In_Seq::value_type;
    size_t n = In.size();
    size_t l = num_blocks(n, _block_size);
    sequence<size_t> Sums(l);
    sliced_for(n, _block_size,
	       [&](size_t i, size_t s, size_t e) {
		 size_t c = 0;
		 for (size_t j = s; j < e; j++) c += (Fl[j] == false);
		 Sums[i] = c;
	       },
	       fl);
    size_t m = scan_inplace(make_slice(Sums), addm<size_t>());
    sequence<T> Out = sequence<T>::uninitialized(n);
    sliced_for(n, _block_size,
	       [&](size_t i, size_t s, size_t e) {
		 size_t c0 = Sums[i];
		 size_t c1 = s + (m - c0);
		 for (size_t j = s; j < e; j++) {
		   if (Fl[j] == false)
		     assign_uninitialized(Out[c0++], In[j]);
		   else
		     assign_uninitialized(Out[c1++], In[j]);
		 }
	       },
	       fl);
    //return std::make_pair(std::move(Out), m);
    parallel_for(0, n, [&](size_t i) {
      In[i] = Out[i];
    });
    return m;
  }

  template<int dim, typename _floatT = double>
  class coordinate {

  protected:

    _floatT data[dim];

  public:

    typedef _floatT floatT;

    coordinate(const _floatT* _data) {

      for (int i = 0; i < dim; ++ i) data[i] = _data[i];

    }

    template<typename T>
    coordinate(T& _data) {

      for (int i = 0; i < dim; ++ i) data[i] = _data[i];

    }

    coordinate() { }

    _floatT& operator[](int i) {

      return data[i];

    }

    template<typename T>
    _floatT& dist(T other) {

      _floatT total = 0.0;

      for (int i = 0; i < dim; ++ i) {
	_floatT tmp = abs(other[i] - data[i]);
	total += tmp * tmp;
      }

      return sqrt(total);
    }

  };


  template<int dim> class boundingBox {

  private:

    coordinate<dim> minCoords, maxCoords;

  public:

    enum relation { exclude, include, overlap };

    coordinate<dim> getMin() { return minCoords; }

    coordinate<dim> getMax() { return maxCoords; }

    void setMin(int i, typename coordinate<dim>::floatT val) { minCoords[i] = val; }

    void setMax(int i, typename coordinate<dim>::floatT val) { maxCoords[i] = val; }

    boundingBox() {

      for (int i = 0; i < dim; ++ i) {
	minCoords[i] = 0;
	maxCoords[i] = 0;
      }

    };

    template<typename T>
    boundingBox(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      if ((e - s) < 1) return;

      minCoords = coordinate<dim>(_input[s]);
      maxCoords = coordinate<dim>(_input[s]);

      for (int j = s; j < e; ++ j) {
	T p = _input[j];
	for (int i = 0; i < dim; ++ i) {
	  minCoords[i] = std::min(p[i], minCoords[i]);
	  maxCoords[i] = std::max(p[i], maxCoords[i]);
	}
      }

    }

    template<typename T>
    void update(T p) {
      for (int i = 0; i < dim; ++ i) {
	minCoords[i] = std::min(p[i], minCoords[i]);
	maxCoords[i] = std::max(p[i], maxCoords[i]);
      }

    }

    template <typename T>
    bool contains(T x) {

      for (int i = 0; i < dim; ++ i) {
	if (x[i] < minCoords[i]) return false;
	if (x[i] > maxCoords[i]) return false;
      }
      return true;

    }

    ~boundingBox() {

    }

    relation compare(boundingBox other) {

      bool exc = false;
      bool inc = true;

      for (int i = 0; i < dim; ++ i) {
	if (maxCoords[i] < other.minCoords[i] || minCoords[i] > other.maxCoords[i])
	  exc = true;

	if (maxCoords[i] < other.maxCoords[i] || minCoords[i] > other.minCoords[i])
	  inc = false;
      }

      if (exc) return exclude;
      else if (inc) return include;
      else return overlap;

    }

  };


  template <typename T, typename _floatT = double>
  class kBuffer: public std::priority_queue<std::pair<_floatT, T>> {

    using queueT = std::priority_queue<std::pair<_floatT, T>>;

  private:
    int k;

  public:

    kBuffer(int _k): queueT(), k(_k) { };

    void insertK(std::pair<_floatT, T> elem) {

      queueT::push(elem);

      if (queueT::size() > k) queueT::pop();

    };

    std::pair<_floatT, T> getK() {

      return queueT::top();

    }

    bool hasK() { return queueT::size() >= k; }
  };


  // template <typename T, typename _floatT = double>
  // class kBuffer {

  // private:
  //   int k;
  //   int writePt;
  //   std::vector<std::pair<_floatT, T>> A;

  // public:

  //   kBuffer(int _k): k(_k), writePt(0) {

  //     A = std::vector<std::pair<_floatT, T>>(k + 1);

  //   };

  //   void insertK(std::pair<_floatT, T> elem) {

  //     A[writePt ++] = elem;

  //     if (writePt >= k + 1) {

  // 	std::nth_element(A.begin(), A.begin() + k, A.end());

  // 	writePt --;

  //     }

  //   };

  //   std::pair<_floatT, T> getK() {

  //     return A[k - 1];

  //   }

  //   bool hasK() { return writePt >= k; }

  //   int size() { return writePt; }

  //   std::pair<_floatT, T> top() {

  //     return A[writePt - 1];

  //   }

  //   void pop() { writePt --; }

  //   void sort() {

  //     std::sort(A.begin(), A.begin() + k);

  //   }

  // };


  template<int dim, typename T, typename _floatT = double>
  class baseNode {

  protected:

    boundingBox<dim> box;

    static const int threshold = 16; // for splitting

  public:

    using floatT = _floatT;

    boundingBox<dim>& getBox() { return box; }

    virtual void setSiblin(baseNode* _siblin) = 0;

    virtual baseNode<dim, T>* getSiblin() = 0;

    virtual bool isRoot() = 0;

    virtual int size() = 0;

    baseNode() {  };

    virtual ~baseNode() { };

    virtual baseNode* insert(container<T>& _input, int s = -1, int e = -1) = 0;

    virtual int erase(container<T>& _input, int s = -1, int e = -1) = 0;

    virtual bool check() = 0;

    virtual void iterate(std::function<void(T)> func) = 0;

    virtual void kNNHelper(T query, kBuffer<T>& buffer) = 0;

    virtual void kNNRangeHelper(T query,
				floatT radius,
				boundingBox<dim> bb,
				kBuffer<T>& buffer) = 0;

    void kNNRange(T query,
		  floatT radius,
		  kBuffer<T>& buffer) {

      boundingBox<dim> bb;

      for (int i = 0; i < dim; ++ i) {

	bb.setMin(i, query[i] - radius);
	bb.setMax(i, query[i] + radius);

      }

      kNNRangeHelper(query, radius, bb, buffer);

    }

  };


  template<int dim, typename T> class internalNode;


  template<int dim, typename T>
  class dataNode: public baseNode<dim, T> { // leaf node

  private:

    baseNode<dim, T>* siblin;

    container<T> data;

    container<char> flag;

    int n;

  public:

    int size() { return n; }

    void setSiblin(baseNode<dim, T>* _siblin) { siblin = _siblin; }

    baseNode<dim, T>* getSiblin() { return siblin; }

    bool isRoot() { return false; }

    dataNode(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      n = e - s;

      baseNode<dim, T>::box = boundingBox<dim>(_input, s, e);

      data = container<T>(n);
      flag = container<char>(n);

      int j = 0;
      for (int i = s; i < e; ++ i) {
	data[j] = _input[i];
	flag[j] = 1;
	j++;
      }

    }

    baseNode<dim, T>* insert(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      if (e - s + size() >= baseNode<dim, T>::threshold) {

	container<T> tmp(e - s + size());

	int i = 0;
	for (int j = s; j < e; ++ j) tmp[i++] = _input[j];

	iterate([&](T p) { tmp[i++] = p; });

	internalNode<dim, T>* newNode = new internalNode<dim, T>(tmp);

	return newNode;

      } else {

	for (int i = s; i < e; ++ i) {

	  baseNode<dim, T>::box.update(_input[i]);

	  data.push_back(_input[i]);

	  flag.push_back(1);
	}

	n += e - s;

	return nullptr;

      }

    }

    int erase(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      if (e - s <= 0) return 0;

      int erased = 0;

      for (int i = s; i < e; ++ i) {
	for (int j = 0; j < data.size(); ++ j) {

	  int k = 0;
	  for (; k < dim; ++ k) {
	    if (data[j][k] != _input[i][k]) {
	      break;
	    }
	  }

	  if (k == dim) {
	    flag[j] = 0;
	    erased ++;
	  }

	}
      }

      n -= erased;

      return erased;

    }

    void iterate(std::function<void(T)> func) {

      int i = 0;
      for (auto exist: flag) {

	if (exist) func(data[i]);
	i++;

      }

    }

    void kNNHelper(T query,
		   kBuffer<T>& buffer) {

      iterate([&](T x) { buffer.insertK({query.dist(x), x}); });

      if (buffer.hasK()) {

	siblin->kNNRange(query, buffer.getK().first, buffer);

      } else {

	siblin->iterate([&](T x) { buffer.insertK({query.dist(x), x}); });

      }

    }

    void kNNRangeHelper(T query,
			typename baseNode<dim, T>::floatT radius,
			boundingBox<dim> bb,
			kBuffer<T>& buffer) {

      iterate([&](T x) {

	typename baseNode<dim, T>::floatT d = query.dist(x);

	if (d <= radius)
	  buffer.insertK({d, x});

      });

    }

    bool check() {

      for (auto x: data) {

	if (!baseNode<dim, T>::box.contains(x)) return false;

      }

      return true;

    }

  };


  template<int dim, typename T>
  class internalNode: public baseNode<dim, T> { // internal node

  protected:

    baseNode<dim, T>* siblin;

    baseNode<dim, T>* left = nullptr;

    baseNode<dim, T>* right = nullptr;

    int splitDim = -1;

    typename baseNode<dim, T>::floatT split = -1;

    int n;

  public:

    int size() { return n; }

    bool isRoot() { return false; }

    void setSiblin(baseNode<dim, T>* _siblin) { siblin = _siblin; }

    baseNode<dim, T>* getSiblin() { return siblin; }

    internalNode(container<T>& _input, int s = -1, int e = -1, int _splitDim = 0):
      splitDim(_splitDim) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      n = e - s;

      if (n < 2)
	throw std::runtime_error("dynKdTree: error, construction requires input size >= 2.");

      baseNode<dim, T>::box = boundingBox<dim>(_input, s, e);

      int leftSize;

      if (spatialMedian) {

	split = (baseNode<dim, T>::box.getMax()[splitDim] +
		 baseNode<dim, T>::box.getMin()[splitDim]) / 2;

	if (e - s < 2000) {
	  auto middle = std::partition(_input.begin() + s, _input.begin() + e, [&](T& elem) {
	    return elem[splitDim] < split;
	  });

	  leftSize = std::distance(_input.begin(), middle) - s;
	} else {
	  parlay::sequence<bool> flag(e - s);

	  parlay::parallel_for(0, e - s, [&](size_t i) {
	    flag[i] = _input[s + i][splitDim] >= split;
	  });

	  leftSize = split_two(_input.cut(s, e), flag);

	}

	if (leftSize == 0 || leftSize == n) {
	  splitDim = (splitDim + 1) % dim;
	  goto computeObjectMedian;
	}

      } else {
      computeObjectMedian:

	if (n < 2000) {
	  std::nth_element(_input.begin() + s,
			   _input.begin() + s + (e - s) / 2,
			   _input.begin() + e,
			   [&](T& a, T& b){
			     return a[splitDim] < b[splitDim];
			   });
	} else {
	  parlay::sort_inplace(_input.cut(s, e), [&](T a, T b){
	    return a[splitDim] < b[splitDim];
	  });
	}

	split = _input[s + (e - s) / 2][splitDim];

	leftSize = s + (e - s) / 2 - s;

      }

      if (leftSize < baseNode<dim, T>::threshold) {

	left = new dataNode<dim, T>(_input, s, s + leftSize);

      } else {

	left = new internalNode(_input, s, s + leftSize, (splitDim + 1) % dim);

      }

      int rightSize = e - s - leftSize;;

      if (rightSize < baseNode<dim, T>::threshold) {

	right = new dataNode<dim, T>(_input, s + leftSize, e);

      } else {

	right = new internalNode(_input, s + leftSize, e, (splitDim + 1) % dim);

      }

      left->setSiblin(right);

      right->setSiblin(left);

    }

    baseNode<dim, T>* insert(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      n += e - s;

      for (int i = s; i < e; ++ i) {

	baseNode<dim, T>::box.update(_input[i]);

      }

      int middleIdx;

      if (e - s < 2000) {
	auto middle = std::partition(_input.begin() + s, _input.begin() + e, [&](T& elem) {
	  return elem[splitDim] < split;
	});

	middleIdx = std::distance(_input.begin(), middle);
      } else {
	parlay::sequence<bool> flag(e - s);

	parlay::parallel_for(0, e - s, [&](size_t i) {
	  flag[i] = _input[s + i][splitDim] >= split;
	});

	middleIdx = s + split_two(_input.cut(s, e), flag);
      }

      auto insertLeft = [&] () {

	auto newLeft = left->insert(_input, s, middleIdx);

	if (newLeft) {
	  delete left;
	  left = newLeft;

	  left->setSiblin(right);
	  right->setSiblin(left);
	}

      };

      auto insertRight = [&] () {

	auto newRight = right->insert(_input, middleIdx, e);

	if (newRight) {
	  delete right;
	  right = newRight;

	  left->setSiblin(right);
	  right->setSiblin(left);
	}

      };

      if (e - s < 2000) {

	insertLeft();
	insertRight();

      } else parlay::par_do(insertLeft, insertRight);

      return nullptr;

    }

    int erase(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      if (e - s <= 0) return 0;

      int middleIdx;

      if (e - s < 2000) {
	auto middle = std::partition(_input.begin() + s, _input.begin() + e, [&](T& elem) {
	  return elem[splitDim] < split;
	});

	middleIdx = std::distance(_input.begin(), middle);
      } else {
	parlay::sequence<bool> flag(e - s);

	parlay::parallel_for(0, e - s, [&](size_t i) {
	  flag[i] = _input[s + i][splitDim] >= split;
	});

	middleIdx = s + split_two(_input.cut(s, e), flag);
      }

      int leftErased, rightErased;

      auto eraseLeft = [&] () {
	leftErased = left->erase(_input, s, middleIdx);
      };

      auto eraseRight = [&] () {
	rightErased = right->erase(_input, middleIdx, e);
      };

      if (e - s < 2000) {

	eraseLeft();
	eraseRight();

      } else parlay::par_do(eraseLeft, eraseRight);

      n -= leftErased + rightErased;

      return leftErased + rightErased;

    }

    ~internalNode() {

      delete left;
      delete right;

    }

    void iterate(std::function<void(T)> func) {

      left->iterate(func);

      right->iterate(func);

    }

    void kNNHelper(T query,
		   kBuffer<T>& buffer) {

      if (query[splitDim] < split) left->kNNHelper(query, buffer);
      else right->kNNHelper(query, buffer);

      if (buffer.hasK()) {

	siblin->kNNRange(query, buffer.getK().first, buffer);

      } else {

	siblin->iterate([&](T x) { buffer.insertK({query.dist(x), x}); });

      }

    }

    void kNNRangeHelper(T query,
			typename baseNode<dim, T>::floatT radius,
			boundingBox<dim> bb,
			kBuffer<T>& buffer) {

      typename boundingBox<dim>::relation rel = bb.compare(baseNode<dim, T>::box);

      if (rel == boundingBox<dim>::relation::include) {

	iterate([&](T x) {

	  typename baseNode<dim, T>::floatT d = query.dist(x);

	  if (d <= radius)
	    buffer.insertK({d, x});

	});

      } else if (rel == boundingBox<dim>::relation::overlap) {

	left->kNNRangeHelper(query, radius, bb, buffer);
	right->kNNRangeHelper(query, radius, bb, buffer);

      }

    }

    bool check() {

      if (left->getSiblin() != right) return false;
      if (right->getSiblin() != left) return false;

      boundingBox<dim> leftBox = left->getBox();
      boundingBox<dim> rightBox = right->getBox();

      if (!baseNode<dim, T>::box.contains(leftBox.getMin()) ||
	  !baseNode<dim, T>::box.contains(leftBox.getMax())) {
	return false;
      }

      return left->check() && right->check();

    }

  };


  template<int dim, typename T>
  class rootNode: public internalNode<dim, T> { // root node

  public:

    baseNode<dim, T>* getSiblin() {

      throw std::runtime_error("dynKdTree: error, cannot get siblin of root\n");

    }

    bool isRoot() { return true; }

    void iterate(std::function<void(T)> func) {

      internalNode<dim, T>::left->iterate(func);

      internalNode<dim, T>::right->iterate(func);

    }

    void kNNHelper(T query,
		   kBuffer<T>& buffer) {

      if (query[internalNode<dim, T>::splitDim] < internalNode<dim, T>::split)
	internalNode<dim, T>::left->kNNHelper(query, buffer);
      else
	internalNode<dim, T>::right->kNNHelper(query, buffer);

    }

    rootNode(container<T>& _input, int s = -1, int e = -1, int _splitDim = 0):
      internalNode<dim, T>(_input, s, e, _splitDim) { };

    container<T> kNN(T query, int k) {

      kBuffer<T> buffer(k);

      kNNHelper(query, buffer);

      int kOut = buffer.size();

      if (kOut < k) {

	std::cout << "dynKdTree: warning, kNN outputs ";
	std::cout << kOut << ", fewer than k.\n";

      }

      auto nns = container<T>(kOut);

      buffer.sort();

      for (int i = 0; i < kOut; ++ i) {
	nns[kOut - 1 - i] = buffer.top().second;
	buffer.pop();
      }

      return nns;
    }

  };


}; // End namespace dynKdTree

}; // End namespace pargeo
