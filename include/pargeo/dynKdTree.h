#pragma once

#include <algorithm>
#include <iostream>
#include <vector>

#include "../parlay/parallel.h"


#ifndef PARLAY_PARALLEL_H_

namespace parlay {

  template <typename Lf, typename Rf>
  inline void par_do(Lf left, Rf right) {
    left();
    right();
  }

}

#endif


namespace dynKdTree {


  template<int dim, typename floatT = double>
  class coordinate {

  protected:

    floatT data[dim];

  public:

    coordinate(floatT* _data) {

      for (int i = 0; i < dim; ++ i) data[i] = _data[i];

    }

    template<typename T>
    coordinate(T& _data) {

      for (int i = 0; i < dim; ++ i) data[i] = _data[i];

    }

    coordinate() { }

    floatT& operator[](int i) {

      return data[i];

    }

  };


  template<int dim> class boundingBox {

  public:

    coordinate<dim> topLeft, lowerRight;

    boundingBox() { };

    template<typename T>
    boundingBox(std::vector<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      if ((e - s) < 2) return;

      topLeft = coordinate<dim>(_input[s]);
      lowerRight = coordinate<dim>(_input[s]);

      for (int j = s; j < e; ++ j) {
	T p = _input[j];
	for (int i = 0; i < dim; ++ i) {
	  topLeft[i] = std::min(p[i], topLeft[i]);
	  lowerRight[i] = std::max(p[i], lowerRight[i]);
	}
      }

    }

    template<typename T>
    void update(T& p) {

      for (int i = 0; i < dim; ++ i) {
	topLeft[i] = std::min(p[i], topLeft[i]);
	lowerRight[i] = std::max(p[i], lowerRight[i]);
      }

    }

    ~boundingBox() {

    }

  };


  template<int dim, typename T> class baseNode {

  protected:

    static const int threshold = 16; // for splitting

    boundingBox<dim> box;

  public:

    virtual bool internal() { return true; }

    baseNode() { }

    virtual ~baseNode() { };

    virtual baseNode* insert(std::vector<T>& _input, int s = -1, int e = -1) {

      return nullptr;

    }

    virtual int erase(std::vector<T>& _input, int s = -1, int e = -1) { return 0; };

  };


  template<int dim, typename T> class splitNode;


  template<int dim, typename T>
  class dataNode: public baseNode<dim, T> { // leaf node

  private:

    std::vector<T> data;

    std::vector<char> flag;

  public:

    int size() { return data.size(); }

    bool internal() { return false; }

    dataNode(std::vector<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      baseNode<dim, T>::box = boundingBox<dim>(_input, s, e);

      data = std::vector<T>();
      flag = std::vector<char>();

      for (int i = s; i < e; ++ i) {
	data.push_back(_input[i]);
	flag.push_back(1);
      }

    }

    baseNode<dim, T>* insert(std::vector<T>& _input, int s = -1, int e = -1) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      for (int i = s; i < e; ++ i) {

	baseNode<dim, T>::box.update(_input[i]);

      }

      if ((e - s) + data.size() >= baseNode<dim, T>::threshold) {

	splitNode<dim, T>* newNode = new splitNode<dim, T>(_input, s, e);

	return newNode;

      } else {

	for (int i = s; i < e; ++ i) {
	  data.push_back(_input[i]);
	  flag.push_back(1);
	}

	return nullptr;

      }

    }

    int erase(std::vector<T>& _input, int s = -1, int e = -1) {

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

  };


  template<int dim, typename T>
  class splitNode: public baseNode<dim, T> { // internal node

  private:

    baseNode<dim, T>* left = nullptr;

    baseNode<dim, T>* right = nullptr;

    int splitDim = -1;

    double split = -1;

    int n;

  public:

    int size() { return n; }

    splitNode(std::vector<T>& _input, int s = -1, int e = -1, int _splitDim = 0):
      splitDim(_splitDim) {

      if (s < 0 || e < 0) {
	s = 0;
	e = _input.size();
      }

      n = e - s;

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

	left = new splitNode(_input, s, s + (e - s) / 2, (splitDim + 1) % dim);

      }

      int rightSize = e - s + (e - s) / 2;

      if (rightSize < baseNode<dim, T>::threshold) {

	right = new dataNode<dim, T>(_input, s + (e - s) / 2, e);

      } else {

	right = new splitNode(_input, s + (e - s) / 2, e, (splitDim + 1) % dim);

      }

    }

    baseNode<dim, T>* insert(std::vector<T>& _input, int s = -1, int e = -1) {

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
	}

      };

      auto insertRight = [&] () {

	auto newRight = right->insert(_input, std::distance(_input.begin(), middle), e);

	if (newRight) {
	  delete right;
	  right = newRight;
	}

      };

      if (e - s < 2000) {

	insertLeft();
	insertRight();

      } else parlay::par_do(insertLeft, insertRight);

      return nullptr;

    }

    int erase(std::vector<T>& _input, int s = -1, int e = -1) {

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

    ~splitNode() {

      delete left;
      delete right;

    }

  };


}; // End namespace dynKdTree
