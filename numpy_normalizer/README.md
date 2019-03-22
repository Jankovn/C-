## numpyNormalizer

### Disclaimer
This program is especially developed to compute _face vertex_ meshes. [Here](https://en.wikipedia.org/wiki/Polygon_mesh#Face-vertex_meshes) is a short description how _face vertex_ meshes are composed.

### Introduction

This program is able to normalize normals, i.e, it creates unit vectors, each set of our mesh that are represented as _numpy array_.

This program has been developed with the [cnpy](https://github.com/rogersce/cnpy) library.

Each file is composed of a set of short descriptions to understand how the program has been developed.

### Environment of work
Linux (ArchLinux)

### Compilation

Without debugging : ```g++ -o nn main.cpp -I./include/ -L/usr/local/lib/ -lcnpy -lz --std=c++14```

With debugging : ```g++ -g -o nn main.cpp -I./include/ -L/usr/local/lib/ -lcnpy -lz --std=c++14```

The debugging line is used with valgrind.

##### Dynamic library linking
In order to avoid a dynamic library linking during the execution of the program, remember to use this line before using the program:

```export LD_LIBRARY_PATH=/usr/local/lib```

### Additional links

- [cnpy library](https://github.com/rogersce/cnpy)
- [Face Vertex mesh](https://en.wikipedia.org/wiki/Polygon_mesh#Face-vertex_meshes)