#pragma once

#include <algorithm>
#include <math.h>
#include <queue>

#include "../parlay/parallel.h"
#include "../parlay/sequence.h"

#ifndef PARLAY_PARALLEL_H_

namespace parlay {

  template <typename Lf, typename Rf>
  inline void par_do(Lf left, Rf right) {
    left();
    right();
  }

}

#endif

#ifndef PARLAY_SEQUENCE_H_

#include <vector>

template <typename T>
using container = std::vector<T>;

#else

template <typename T>
using container = parlay::sequence<T>;

#endif

namespace pargeo {

namespace dynKdTree {


  template<int dim, typename _floatT = double>
  class coordinate {

  protected:

    _floatT data[dim];

  public:

    typedef _floatT floatT;

    coordinate(_floatT* _data) {

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
    _floatT& dist(const T& other) {

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
    void update(T& p) {
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


  template <typename T>
  class kBuffer: public std::priority_queue<std::pair<double, T>> {

    using queueT = std::priority_queue<std::pair<double, T>>;

  private:
    int k;

  public:

    kBuffer(int _k): queueT(), k(_k) { };

    void insertK(std::pair<double, T> elem) {

      queueT::push(elem);

      if (queueT::size() > k) queueT::pop();

    };

    std::pair<double, T> getK() {

      return queueT::top();

    }

    bool hasK() { return queueT::size() >= k; }
  };


  template<int dim, typename T> class baseNode {

  protected:

    boundingBox<dim> box;

    static const int threshold = 16; // for splitting

  public:

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

    void kNNRange(T query,
		  double radius,
		  kBuffer<T>& buffer) {

      iterate([&](T x) {

	double d = query.dist(x);

	if (d <= radius)
	  buffer.insertK({d, x});

      });

    }

  };


  template<int dim, typename T> class internalNode;


  template<int dim, typename T>
  class dataNode: public baseNode<dim, T> { // leaf node

  private:

    baseNode<dim, T>* siblin;

    container<T> data;

    container<char> flag;

  public:

    int size() { return data.size(); }

    void setSiblin(baseNode<dim, T>* _siblin) { siblin = _siblin; }

    baseNode<dim, T>* getSiblin() { return siblin; }

    bool isRoot() { return false; }

    dataNode(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      baseNode<dim, T>::box = boundingBox<dim>(_input, s, e);

      data = container<T>();
      flag = container<char>();

      for (int i = s; i < e; ++ i) {
	data.push_back(_input[i]);
	flag.push_back(1);
      }

    }

    baseNode<dim, T>* insert(container<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      for (int i = s; i < e; ++ i) {

	baseNode<dim, T>::box.update(_input[i]);

      }

      if ((e - s) + data.size() >= baseNode<dim, T>::threshold) {

	internalNode<dim, T>* newNode = new internalNode<dim, T>(_input, s, e);

	return newNode;

      } else {

	for (int i = s; i < e; ++ i) {
	  data.push_back(_input[i]);
	  flag.push_back(1);
	}

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

      return erased;

    }

    void iterate(std::function<void(T)> func) {

      int i = 0;
      for (auto exist: flag) {

	if (exist) func(data[i++]);

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

    double split = -1;

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

      std::nth_element(_input.begin() + s,
		       _input.begin() + s + (e - s) / 2,
		       _input.begin() + e,
		       [&](T& a, T& b){
			 return a[splitDim] < b[splitDim];
		       });

      split = _input[s + (e - s) / 2][splitDim];

      int leftSize = s + (e - s) / 2 - s;

      if (leftSize < baseNode<dim, T>::threshold) {

	left = new dataNode<dim, T>(_input, s, s + (e - s) / 2);

      } else {

	left = new internalNode(_input, s, s + (e - s) / 2, (splitDim + 1) % dim);

      }

      int rightSize = e - s + (e - s) / 2;

      if (rightSize < baseNode<dim, T>::threshold) {

	right = new dataNode<dim, T>(_input, s + (e - s) / 2, e);

      } else {

	right = new internalNode(_input, s + (e - s) / 2, e, (splitDim + 1) % dim);

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

      auto middle = std::partition(_input.begin() + s, _input.begin() + e, [&](T& elem) {
	return elem[splitDim] < split;
      });

      auto insertLeft = [&] () {

	auto newLeft = left->insert(_input, s, std::distance(_input.begin(), middle));

	if (newLeft) {
	  delete left;
	  left = newLeft;

	  left->setSiblin(right);
	  right->setSiblin(left);
	}

      };

      auto insertRight = [&] () {

	auto newRight = right->insert(_input, std::distance(_input.begin(), middle), e);

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

      auto middle = std::partition(_input.begin() + s, _input.begin() + e, [&](T& elem) {
	return elem[splitDim] < split;
      });

      int leftErased, rightErased;

      auto eraseLeft = [&] () {
	leftErased = left->erase(_input, s, std::distance(_input.begin(), middle));
      };

      auto eraseRight = [&] () {
	rightErased = right->erase(_input, std::distance(_input.begin(), middle), e);
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

    template<typename F>
    void iterate(F func) {

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

	std::cout << "dynKdTree: warning, kNN outputs fewer than k.\n";

      }

      auto nns = container<T>(kOut);

      for (int i = 0; i < kOut; ++ i) {
	nns[kOut - 1 - i] = buffer.top().second;
	buffer.pop();
      }

      return nns;
    }

  };


}; // End namespace dynKdTree

}; // End namespace pargeo
