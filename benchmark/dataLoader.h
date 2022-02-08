#pragma once

#include "dataset/uniform.h"
#include "pargeo/pointIO.h"
#include "parlay/sequence.h"

size_t N = 1000000;

/* Synthetic data sets */

template<int dim>
parlay::sequence<pargeo::point<dim>> inSphere(size_t n) {
  return pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 0);}

template<int dim>
parlay::sequence<pargeo::fpoint<dim>> inSphere_float(size_t n) {
  return pargeo::uniformInPolyPoints<dim, pargeo::fpoint<dim>>(n, 0);}

template<int dim>
parlay::sequence<pargeo::point<dim>> onSphere(size_t n) {
  return pargeo::uniformOnPolyPoints<dim, pargeo::point<dim>>(n, 0, 0.05);}

template<int dim>
parlay::sequence<pargeo::fpoint<dim>> onSphere_float(size_t n) {
  return pargeo::uniformOnPolyPoints<dim, pargeo::fpoint<dim>>(n, 0, 0.05);}

template<int dim>
parlay::sequence<pargeo::point<dim>> inCube(size_t n) {
  return pargeo::uniformInPolyPoints<dim, pargeo::point<dim>>(n, 1);}

template<int dim>
parlay::sequence<pargeo::fpoint<dim>> inCube_float(size_t n) {
  return pargeo::uniformInPolyPoints<dim, pargeo::fpoint<dim>>(n, 1);}

template<int dim>
parlay::sequence<pargeo::point<dim>> onCube(size_t n) {
  return pargeo::uniformOnPolyPoints<dim, pargeo::point<dim>>(n, 1, 0.05);}

template<int dim>
parlay::sequence<pargeo::fpoint<dim>> onCube_float(size_t n) {
  return pargeo::uniformOnPolyPoints<dim, pargeo::fpoint<dim>>(n, 1, 0.05);}

/* Graphics 3d data sets */

parlay::sequence<pargeo::point<3>> bunny_3d_36k() {
  try {
    auto filePath = "datasets/bunny.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: bunny data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> bunny_3d_36k_float() {
  try {
    auto filePath = "datasets/bunny.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: bunny data set is not available\n";
    return inSphere_float<3>(N);
  }
}

parlay::sequence<pargeo::point<3>> armadillo_3d_173k() {
  try {
    auto filePath = "datasets/armadillo.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: armadillo data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> armadillo_3d_173k_float() {
  try {
    auto filePath = "datasets/armadillo.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: armadillo data set is not available\n";
    return inSphere_float<3>(N);
  }
}

parlay::sequence<pargeo::point<3>> dragon_3d_438k() {
  try {
    auto filePath = "datasets/dragon.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: dragon data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> dragon_3d_438k_float() {
  try {
    auto filePath = "datasets/dragon.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: dragon data set is not available\n";
    return inSphere_float<3>(N);
  }
}


parlay::sequence<pargeo::point<3>> buddha_3d_544k() {
  try {
    auto filePath = "datasets/buddha.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: buddha data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> buddha_3d_544k_float() {
  try {
    auto filePath = "datasets/buddha.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: buddha data set is not available\n";
    return inSphere_float<3>(N);
  }
}

parlay::sequence<pargeo::point<3>> thaiStatue_3d_3_5m() {
  try {
    auto filePath = "datasets/thai-statue.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: thaiStatue data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> thaiStatue_3d_3_5m_float() {
  try {
    auto filePath = "datasets/thai-statue.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: thaiStatue data set is not available\n";
    return inSphere_float<3>(N);
  }
}

parlay::sequence<pargeo::point<3>> asianDragon_3d_3_6m() {
  try {
    auto filePath = "datasets/asian-dragon.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: asianDragon data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> asianDragon_3d_3_6m_float() {
  try {
    auto filePath = "datasets/asian-dragon.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: asianDragon data set is not available\n";
    return inSphere_float<3>(N);
  }
}

parlay::sequence<pargeo::point<3>> lucy_3d_3_14m() {
  try {
    auto filePath = "datasets/lucy.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::point<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::point<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: lucy data set is not available\n";
    return inSphere<3>(N);
  }
}

parlay::sequence<pargeo::fpoint<3>> lucy_3d_3_14m_float() {
  try {
    auto filePath = "datasets/lucy.csv";
    int dim = pargeo::pointIO::readHeader(filePath);
    parlay::sequence<pargeo::fpoint<3>> P =
      pargeo::pointIO::readPointsFromFile<pargeo::fpoint<3>>(filePath);
    return P;
  } catch (const std::runtime_error& error) {
    std::cout << "warning: lucy data set is not available\n";
    return inSphere_float<3>(N);
  }
}
