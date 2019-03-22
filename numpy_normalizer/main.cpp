#include <stdexcept>

#include "numpyNormalizer.hpp"

int			main(int argc, char** argv, char** env) {
  // ---------- Exploration of the cnpy library ----------

  cnpy::NpyArray	faces;
  cnpy::NpyArray	vertices;

  try {
    faces = cnpy::npy_load(NN::facesPath);
    vertices = cnpy::npy_load(NN::verticesPath);
  } catch(const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
  }

  // std::cout << "Faces:" << std::endl;
  // for(unsigned int i = 0; i < faces.shape.size(); ++i) {
  //   std::cout << "shape[" << i << "]: " << faces.shape[i] << std::endl;
  // }
  // std::cout << faces.num_vals << std::endl;
  // std::cout << faces.word_size << std::endl;

  // std::cout << "Vertices:" << std::endl;
  // for(unsigned int i = 0; i < vertices.shape.size(); ++i) {
  //   std::cout << "vertices[" << i << "]: " << vertices.shape[i] << std::endl;
  // }
  // std::cout << vertices.num_vals << std::endl;
  // std::cout << vertices.word_size << std::endl;

  // ---------- Example for one iteration of the final program ----------

  // int i = faces.data<int>()[0];
  // int j = faces.data<int>()[1];
  // int k = faces.data<int>()[2];

  // double i1 = vertices.data<double>()[i];
  // double i2 = vertices.data<double>()[i + 1];
  // double i3 = vertices.data<double>()[i + 2];

  // double j1 = vertices.data<double>()[j];
  // double j2 = vertices.data<double>()[j + 1];
  // double j3 = vertices.data<double>()[j + 2];

  // double k1 = vertices.data<double>()[k];
  // double k2 = vertices.data<double>()[k + 1];
  // double k3 = vertices.data<double>()[k + 2];
  
  // std::cout << faces.data<int>()[0] << ' ' << faces.data<int>()[1] << ' ' << faces.data<int>()[2] << std::endl;
  // std::cout << i1 << ' ' << i2 << ' ' << i3 << std::endl;
  // std::cout << j1 << ' ' << j2 << ' ' << j3 << std::endl;
  // std::cout << k1 << ' ' << k2 << ' ' << k3 << std::endl;

  // std::cout << "C-A" << std::endl;
  // std::cout << k1 - i1 << std::endl;
  // std::cout << k2 - i2 << std::endl;
  // std::cout << k3 - i3 << std::endl;

  // std::cout << "B-A" << std::endl;
  // std::cout << j1 - i1 << std::endl;
  // std::cout << j2 - i2 << std::endl;
  // std::cout << j3 - i3 << std::endl;

  // NN::t_vec3D<double, double> u(k1 - i1, k2 - i2, k3 - i3);
  // NN::t_vec3D<double, double> v(j1 - i1, j2 - i2, j3 - i3);

  // std::cout << u;
  // std::cout << v;
  
  // NN::t_vec3D<double, double> cross = u * v;;

  // cross.normalize();
  
  // std::cout << cross;

  // ---------- Concrete example of the final aim of the project ----------

  NN::Mesh<int, double> m(&faces, &vertices);

  m.computeNormals();
  m.normalizeNormals();
  m.dumpXVec(10);
  m.saveNormalized();
 
  // -------------------- Check if the created files contain the right information --------------------

  // cnpy::NpyArray	faces;
  // cnpy::NpyArray	vertices;

  // try {
  //   faces = cnpy::npy_load(NN::normalizedFacesPath);
  //   vertices = cnpy::npy_load(NN::normalizedVerticesPath);
  // } catch(const std::runtime_error& err) {
  //   std::cerr << err.what() << std::endl;
  // }

  // // std::cout << "Faces:" << std::endl;
  // for(unsigned int i = 0; i < faces.shape.size(); ++i) {
  //   std::cout << "shape[" << i << "]: " << faces.shape[i] << std::endl;
  // }
  // // std::cout << faces.num_vals << std::endl;
  // // std::cout << faces.word_size << std::endl;

  // // std::cout << "Vertices:" << std::endl;
  // for(unsigned int i = 0; i < vertices.shape.size(); ++i) {
  //   std::cout << "vertices[" << i << "]: " << vertices.shape[i] << std::endl;
  // }

  // double i = faces.data<double>()[0];
  // double j = faces.data<double>()[1];
  // double k = faces.data<double>()[2];

  // double l = vertices.data<double>()[0];
  // double m = vertices.data<double>()[1];
  // double n = vertices.data<double>()[2];

  // std::cout << i << std::endl;
  // std::cout << j << std::endl;
  // std::cout << k << std::endl;
  // std::cout << l << std::endl;
  // std::cout << m << std::endl;
  // std::cout << n << std::endl;

  return (0);
}
