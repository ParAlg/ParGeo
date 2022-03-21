#ifndef TEST_SHARED_BASICSTRUCTURE2D_H
#define TEST_SHARED_BASICSTRUCTURE2D_H

#include <gtest/gtest.h>
//#include "common/geometryIO.h"
#include "pargeo/pointIO.h"
#include "pargeo/point.h"

template <class T>
static auto KEEP_EVEN(const parlay::sequence<T>& seq) {
  // construct the points to delete
  auto num_even = (seq.size() + 1) / 2;
  parlay::sequence<T> to_remove(num_even);
  for (size_t i = 0; i < num_even; i++) {
    to_remove[i] = seq[2 * i];
  }
  return to_remove;
}

template <class T>
static auto KEEP_ODD(const parlay::sequence<T>& seq) {
  // construct the points to delete
  auto num_odd = (seq.size()) / 2;
  parlay::sequence<T> to_remove(num_odd);
  for (size_t i = 0; i < num_odd; i++) {
    to_remove[i] = seq[2 * i + 1];
  }
  return to_remove;
}
template <class Tree>
class BasicStructure2D : public ::testing::Test {
  static const int dim = 2;

 public:
  static constexpr double POINT_ARR_2[2][2] = {{0, 0}, {1, 1}};
  static constexpr double POINT_ARR_8[8][2] = {
      {0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {7, 7}};

  static Tree CONSTRUCT_2D_TREE(const double point_array[][2], int num_points) {
    std::vector<pargeo::point<dim>> point_vec(num_points);
    for (int i = 0; i < num_points; i++)
      point_vec[i] = pargeo::point<dim>(point_array[i]);
    parlay::sequence<pargeo::point<dim>> points(point_vec.begin(), point_vec.end());

    return Tree(points);
  }

  static Tree CONSTRUCT_2D_SIZE_2() {
    constexpr int tree_height = 2;
    constexpr int num_points = 1 << (tree_height - 1);
    return CONSTRUCT_2D_TREE(POINT_ARR_2, num_points);
  }

  static Tree CONSTRUCT_2D_SIZE_8() {
    constexpr int tree_height = 4;
    constexpr int num_points = 1 << (tree_height - 1);
    return CONSTRUCT_2D_TREE(POINT_ARR_8, num_points);
  }

  static parlay::sequence<pargeo::point<dim>> RESOURCES_1000() {
    const char* test_file = "../resources/2d-UniformInSphere-1k.pbbs";
    //int check_dim = readDimensionFromFile(test_file);
    int check_dim = pargeo::pointIO::readHeader(test_file);
    if (check_dim != dim) throw std::runtime_error("Invalid input file!");

    return pargeo::pointIO::readPointsFromFile<pargeo::point<dim>>(test_file);
  }

  static Tree CONSTRUCT_RESOURCES_1000() {
    const char* test_file = "../resources/2d-UniformInSphere-1k.pbbs";
    //int check_dim = readDimensionFromFile(test_file);
    int check_dim = pargeo::pointIO::readHeader(test_file);
    if (check_dim != dim) throw std::runtime_error("Invalid input file!");

    parlay::sequence<pargeo::point<dim>> points =
      pargeo::pointIO::readPointsFromFile<pargeo::point<dim>>(test_file);
    return Tree(points);
  }

  static Tree CONSTRUCT_RESOURCES_1000_EVEN() {
    const char* test_file = "../resources/2d-UniformInSphere-1k.pbbs";
    //int check_dim = readDimensionFromFile(test_file);
    int check_dim = pargeo::pointIO::readHeader(test_file);
    if (check_dim != dim) throw std::runtime_error("Invalid input file!");

    parlay::sequence<pargeo::point<dim>> points =
      pargeo::pointIO::readPointsFromFile<pargeo::point<dim>>(test_file);
    const auto to_insert = KEEP_EVEN(points);

    Tree ret(10);
    ret.insert(to_insert.cut(0, to_insert.size()));
    return ret;
  }
};

#endif  // TEST_SHARED_BASICSTRUCTURE2D_H