#pragma once

# include "cnpy.h"

# include <iostream>
# include <memory>
# include <array>
# include <cmath>

/*
** This Namespace embeds all the C++ statements of the program.
*/
namespace NN {

  /*
  ** Entries that point to the files we have to compute.
  */
  const std::string facesPath = "./resources/bunny_faces.npy";
  const std::string verticesPath = "./resources/bunny_vertices.npy";

  /*
  ** Path to filenames where normalized values will be stored.
  */
  const std::string normalizedFacesPath = "./results/face_normals.npy";
  const std::string normalizedVerticesPath = "./results/vertex_normals.npy";

  /*
  ** This variable correponds to the number of steps we have to do in order
  ** to jump between each vector contained within the array (since we are in a R³)
  ** space given by cnpy::NpyArray;;data().
  */
  const unsigned int step = 3;

  /*
  ** Variable especially created for the creation of npy files after performing
  ** the normalization of normals for both faces and vertices. The normalized normals
  ** of the face, and normalized normals of vertices will be of type double, hence the value 8.
  */
  const size_t word_size = 8;

  /*
  ** Enum class that will allow us to differentiate to which the actual vector belongs to.
  */
  enum class Component {FACE, VERTEX};

  /*
  ** This templated struture will reproduce the behavior of vectors according to our needs
  ** for this challenge.
  **
  ** It is intended that this structure is composed of two templated parameters.
  ** Indeed, according to the needs of the challenge and by analyzing the constitution
  ** of the mesh, we can deduce that the mesh adopts the "Face-Vertex" behavior. Therefore,
  ** one set of values is composed of integers and the other is composed of floating point numbers.
  **
  ** C++: This structure follow the rule of three.
  ** T will usually represent all kind (type) of integer values.
  ** U will usually represent all kind (type) of floating point values.
  */
  template<typename T, typename U>
  struct	s_vector3D {
  public:
    /*
    ** These coordinates will represent the value retrieved from cnpy::NpyArray;;data().
    ** These values have been gathered by aggregating the initial values from
    ** cnpy's pointer per stack of 3 (R³ environment).
    **
    ** This action has been taken, although it might give an
    ** overhead to the initial aim of the challenge, for a better readability
    ** instead of playing with pointers and offsets.
    */
    T		_x; // Corresponds to each multiple of 3.
    T		_y; // Corresponds to each multiple of 3 + 1.
    T		_z; // Corresponds to each multiple of 3 + 2.
    /*
    ** Normal vector components.
    **
    ** To compute these normals, we will retrieve the vertices that are assigned
    ** to the face of the mesh.
    **
    ** Therefore, according to the given vector,
    ** we will either perform an assignment (if it is related to a face),
    ** or an aggregation (if it is related to a vertex).
    */
    U		_xNormal;
    U		_yNormal;
    U		_zNormal;
    /*
    ** Normalized normal vector components.
    */
    U		_xNormalized;
    U		_yNormalized;
    U		_zNormalized;
  public:
    /*
    ** Default Constructor.
    */
    constexpr s_vector3D() noexcept
    : _x(0), _y(0), _z(0),
      _xNormal(0), _yNormal(0), _zNormal(0),
      _xNormalized(0), _yNormalized(0), _zNormalized(0)
    {}

    constexpr s_vector3D(const T& x, const T& y, const T& z)
      : _x(x), _y(y), _z(z),
	_xNormal(0), _yNormal(0), _zNormal(0),
	_xNormalized(0), _yNormalized(0), _zNormalized(0)
    {}

    /*
    ** Constructor used when performing cross-product.
    */
    constexpr s_vector3D(const T& x, const T& y, const T& z,
			 const T& xn, const T& yn, const T& zn)
      : _x(x), _y(y), _z(z),
	_xNormal(xn), _yNormal(yn), _zNormal(zn),
	_xNormalized(0), _yNormalized(0), _zNormalized(0)
    {}

    /*
    ** Copy constructor.
    */
    constexpr s_vector3D(const struct s_vector3D& v) noexcept
      : _x(v._x), _y(v._y), _z(v._z),
	_xNormal(v._xNormal), _yNormal(v._yNormal), _zNormal(v._zNormal),
	_xNormalized(v._xNormalized), _yNormalized(v._yNormalized), _zNormalized(v._zNormalized)
    {}

    ~s_vector3D() = default;

    /*
    ** This method will perform an assignment on each vector normal's components.
    ** Since the normal of a face is the same on the entire face.
    */
    constexpr void addNormal(const U& x, const U& y, const U& z) noexcept {
      _xNormal = x;
      _yNormal = y;
      _zNormal = z;
    }

    /*
    ** This method will perform an aggregation on each vector normal's components.
    ** Since a vertex can be shared with several faces, hence the aggregation.
    */
    constexpr void aggregateNormal(const struct s_vector3D& n) noexcept {
      _xNormal += n._xNormal;
      _yNormal += n._yNormal;
      _zNormal += n._zNormal;
    }

    /*
    ** Normalize the vector normal.
    */
    constexpr void	normalize() noexcept {
      U			norm;

      // Compute the norm (magnitude) of the vector.
      norm = sqrt((_xNormal * _xNormal + _yNormal * _yNormal + _zNormal * _zNormal));

      // Normalize each component of the vector.
      _xNormalized = _xNormal / norm;
      _yNormalized = _yNormal / norm;
      _zNormalized = _zNormal / norm;
    }

    /*
    ** Copy assignment.
    */
    s_vector3D<T, U>& operator=(const struct s_vector3D&);

    /*
    ** Arithmetic operators.
    */
    s_vector3D<T, U> operator-(const struct s_vector3D&);

    s_vector3D<T, U> operator*(const struct s_vector3D&); // cross-product

  };

  /*
  ** Typical vector subtraction.
  */
  template<typename T, typename U>
  s_vector3D<T, U> s_vector3D<T, U>::operator-(const struct s_vector3D& v) {
    return (s_vector3D((_x - v._x),
		       (_y - v._y),
		       (_z - v._z)));
  };

  /*
  ** Cross-product between two vectors.
  */
  template<typename T, typename U>
  s_vector3D<T, U> s_vector3D<T, U>::operator*(const struct s_vector3D& v) {
    return (s_vector3D(0, 0, 0,
		       ((_y * v._z) - (v._y * _z)),
		       ((_x * v._z) - (v._x * _z)) * -1,
		       ((_x * v._y) - (v._x * _y))));
  };

  /*
  ** Definition of the copy assignment operator.
  */
  template<typename T, typename U>
  s_vector3D<T, U>& s_vector3D<T, U>::operator=(const struct s_vector3D& v) {
    _x = v._x;
    _y = v._y;
    _z = v._z;
    _xNormal = v._xNormal;
    _yNormal = v._yNormal;
    _zNormal = v._zNormal;
    _xNormalized = v._xNormalized;
    _yNormalized = v._yNormalized;
    _zNormalized = v._zNormalized;
    return (*this);
  };

  /*
  ** Definition of a pre-formatted output for our designed structure vector.
  */
  template<typename T, typename U>
  std::ostream&	operator<<(std::ostream& os, const struct s_vector3D<T, U>& v) {
    os << "{ x: " << v._x << " y: " << v._y << " z: " << v._z << " }" << std::endl;
    os << "{ xNormal: " << v._xNormal << " yNormal: " << v._yNormal
       << " zNormal: " << v._zNormal << " }" << std::endl;
    os << "{ xNormalized: " << v._xNormalized << " yNormalized: " << v._yNormalized
       << " zNormalized: " << v._zNormalized << " }" <<std::endl;
    return (os);
  }

  /*
  ** Define an alias for a better readability and understanding of the structure we created.
  */
  template<typename T, typename U>
  using t_vec3D = struct s_vector3D<T, U>;


  /*
  ** Templated Class that embeds one component that initially constitutes a mesh.
  **
  ** typename T: Corresponds to the type assigned with cnpy::NpyArray::data<T>()
  ** typename U: Corresponds to the type given for t_vec3D.
  **
  */
  template<typename T, typename U>
  class				MeshComponent {
  private:
    // Contains the "packed" values of NpyArray::data()
    std::vector<t_vec3D<T, U>>	_arrayVec;
    // Identify the type of component we are working on.
    Component			_c;
    // Initial size of the component according to NpyArray::shape[0]
    size_t			_size;
    // Array that contains the Normalized normal of the component.
    std::unique_ptr<U[]>	_UPtr;
  public:

    /*
    ** Default Constructor.
    ** Avoid the invocation withtout specifying the components of a mesh.
    */
    MeshComponent() = delete;

    MeshComponent(const T* ptr, size_t s, Component c)
      : _arrayVec((s)), _c(c), _size(s),
	_UPtr(std::make_unique<U[]>(s * NN::step)) {
      unsigned int jmp = 0;

      /*
      ** Reassign the values contained in the raw pointer given by cnpy,
      ** and packing them in t_vec3D by respecting their offset.
      */
      for (unsigned int i = 0; i < _size; ++i) {
	_arrayVec[i]._x = ptr[jmp];
	_arrayVec[i]._y = ptr[jmp + 1];
	_arrayVec[i]._z = ptr[jmp + 2];
	jmp += NN::step;
      }
    }

    ~MeshComponent() = default;

    /*
    ** Copy constructor.
    */
    MeshComponent(const MeshComponent&) = delete;

    /*
    ** Copy assignment.
    */
    MeshComponent<T, U>& operator=(const MeshComponent&) = delete;

    /*
    ** Subscript operator.
    */
    t_vec3D<T, U>& operator[](unsigned int);

    /*
    ** Return the "name" of the component as an enum.
    */
    inline Component	getComponent() const noexcept {
      return (_c);
    }

    /*
    ** Return the size of component (shape).
    */
    inline size_t	getSize() const noexcept {
      return (_size);
    }

    /*
    ** Assign the final normalized value of a vector at the corresponding
    ** step in the UPtr. This UPtr will be used when we will write the values
    ** contained in, in the npy file.
    */
    void addNormalizedToUPtr(const U& xn, const U& yn, const U& zn,
			     unsigned int step) noexcept {
      _UPtr[step] = xn;
      _UPtr[step + 1] = yn;
      _UPtr[step + 2] = zn;
    }

    /*
    ** Save each computed set of normalized normals of each mesh's component in an npy file.
    */
    void save() const noexcept {
      cnpy::npy_save((_c == NN::Component::FACE ? NN::normalizedFacesPath : NN::normalizedVerticesPath),
      		     _UPtr.get(),
      		     {_size, NN::step});
    }

  };

  /*
  ** Implemented for a better accessiblity.
  */
  template<typename T, typename U>
  t_vec3D<T, U>&	MeshComponent<T, U>::operator[](unsigned int index) {
    return (_arrayVec[index]);
  }

  /*
  ** Templated Class that represents the integrity, in terms of components, of a mesh.
  ** Thus, it embeds the essential correctly formatted components for the challenge,
  ** on top of extra methods to satisfy the requirements of the challenge.
  **
  ** This class acts as a "Face-vertex mesh".
  */
  template<typename T, typename U>
  class			Mesh {
  private:
    // Attribute corresponding to the face component.
    MeshComponent<T, U>	_face;
    // Attribute corresponding to the vertex component.
    MeshComponent<U, U>	_vertex;
  public:
    /*
    ** Default Constructor.
    ** Avoid the invocation withtout specifying the components of a mesh.
    */    
    Mesh() = delete;

    Mesh(cnpy::NpyArray* f, cnpy::NpyArray* v)
      : _face(f->data<T>(), f->shape[0], Component::FACE),
	_vertex(v->data<U>(), v->shape[0], Component::VERTEX)
    {}

    ~Mesh() = default;

    /*
    ** Copy constructor.
    */
    Mesh(const Mesh&) = delete;

    /*
    ** Copy assignment.
    */
    Mesh<T, U>& operator=(const Mesh&) = delete;

    /*
    ** Iterate through all the faces to compute the normals of
    ** both faces and vertices.
    */
    void computeNormals() noexcept {
      for (unsigned int i = 0; i < _face.getSize(); ++i) {
	// Retrieve the references of "face(vertex1, vertex2, vertex3)".
	auto& actualFace = _face[i];
	auto& v1 = _vertex[actualFace._x];
	auto& v2 = _vertex[actualFace._y];
	auto& v3 = _vertex[actualFace._z];

	// Determine vector U and vector V to compute the cross-product.
	auto uVec = v3 - v1;
	auto vVec = v2 - v1;

	// Cross-product between U and V.
	auto cp = uVec * vVec;

	// Assign or aggregate the normals based on the component the vector is related to.
	actualFace.addNormal(cp._xNormal, cp._yNormal, cp._zNormal);
	v1.aggregateNormal(cp);
	v2.aggregateNormal(cp);
	v3.aggregateNormal(cp);
      }
    }

    /*
    ** Iterate through all the faces to compute the normalized normals of
    ** both faces and vertices.
    */
    void normalizeNormals() noexcept {
      unsigned int step = 0;

      for (unsigned int i = 0; i < _face.getSize(); ++i) {
	// Retrieve the references of "face(vertex1, vertex2, vertex3)".
	auto& actualFace = _face[i];
	auto& v1 = _vertex[actualFace._x];
	auto& v2 = _vertex[actualFace._y];
	auto& v3 = _vertex[actualFace._z];

	// Normalize the normals.
	actualFace.normalize();
	v1.normalize();
	v2.normalize();
	v3.normalize();

	// Set at the right position each normalized vector (either face or vertex) in the UPtr.
	_face.addNormalizedToUPtr(actualFace._xNormalized, actualFace._yNormalized,
				  actualFace._zNormalized, step);
	_vertex.addNormalizedToUPtr(v1._xNormalized, v1._yNormalized,
				    v1._zNormalized, static_cast<unsigned int>(actualFace._x));
	_vertex.addNormalizedToUPtr(v2._xNormalized, v2._yNormalized,
				    v2._zNormalized, static_cast<unsigned int>(actualFace._y));
	_vertex.addNormalizedToUPtr(v3._xNormalized, v3._yNormalized,
				    v3._zNormalized, static_cast<unsigned int>(actualFace._z));

	// Follow the same semantics of npy to store values.
	step += NN::step;
      }
    }

    /*
    ** Call each component of the mesh to save the state of its set of normalized normals
    ** in a compliant way with the cnpy library.
    */
    void	saveNormalized() const noexcept {
      _face.save();
      _vertex.save();
    }

    /*
    ** Print with a "readable" format the initial vector's components/normals/normalized.
    **
    ** A set of faces/vertices have been chosen to verify, thanks to WolframAlpha,
    ** if their vector (normalized normals) has a length of 1. That's the case.
    */
    void dumpXVec(unsigned int nb) {
      for (unsigned int i = 0; i < nb; ++i) {
	auto& actualFace = _face[i];
	auto& v1 = _vertex[actualFace._x];
	auto& v2 = _vertex[actualFace._y];
	auto& v3 = _vertex[actualFace._z];

	std::cout << "---" << i << "---" << std::endl;

	std::cout << "-- Face --" << std::endl;
	std::cout << actualFace;
	std::cout << "-- --" << std::endl;

	std::cout << "-- Vertices--" << std::endl;
	std::cout << v1;
	std::cout << '-' << std::endl;
	std::cout << v2;
	std::cout << '-' << std::endl;
	std::cout << v3;
	std::cout << "-- --" << std::endl;

	std::cout << "--- ---" << std::endl;
      }
    }

  };

};
