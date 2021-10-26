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

  };


  template<int dim> class boundingBox {

  public:

    template<typename T>
    boundingBox(std::vector<T> _input) {

    }

    ~boundingBox() {

    }

  };


  template<int dim, typename T> class baseNode {

  protected:

    boundingBox<dim> box;

  public:

    virtual bool internal() { return true; }

    baseNode() { };

    virtual ~baseNode() { };

  };


  template<int dim, typename T>
  class splitNode: public baseNode<dim, T> {

  public:

    splitNode() {

    }

  };


  template<int dim, typename T>
  class dataNode: public baseNode<dim, T> {

    std::vector<T> data;

  public:

    bool internal() { return false; }

    dataNode(const std::vector<T>& _data) {

    }

  };


}; // End namespace dynKdTree
