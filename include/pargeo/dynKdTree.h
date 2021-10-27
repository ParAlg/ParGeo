#pragma once

#include <vector>


// Simple dynamic kd-tree

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

    template<typename T>
    boundingBox(std::vector<T> _input) {

      if (_input.size() < 2) return;

      topLeft = coordinate<dim>(_input[0]);
      lowerRight = coordinate<dim>(_input[0]);

      for (auto p: _input) {
	for (int i = 0; i < dim; ++ i) {
	  topLeft[i] = std::min(p[i], topLeft[i]);
	  lowerRight[i] = std::max(p[i], lowerRight[i]);
	}
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

    baseNode() { };

    virtual ~baseNode() { };

  };


  template<int dim, typename T>
  class splitNode: public baseNode<dim, T> { // internal node

    baseNode<dim, T>* left = nullptr;

    baseNode<dim, T>* right = nullptr;

    int splitDim = -1;

  public:

    splitNode(vector<T>& _input, _splitDim = 0):
      splitDim(_splitDim) {

      // nth_element

    }

  };


  template<int dim, typename T>
  class dataNode: public baseNode<dim, T> { // leaf node

    std::vector<T> data;

  public:

    bool internal() { return false; }

    dataNode(const std::vector<T>& _data) {

    }

  };


}; // End namespace dynKdTree
