# Source
 A C++ implementation in the [Boost Graph Library](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/index.html) of two algorithms for path
 coloring and path list-coloring plane graphs.

## Directories:

### path_coloring
 This is the main directory of the project and contains header files for the
 three algorithms implemented. It also includes some headers for helper functions
 and data structures shared between the different algorithms. 

### examples
 This directory contains examples displaying how to call each algorithm
 implemented.

### unit_tests
 This directory contains a unit testing suite for the headers in the
 path_coloring directory.

### visualization
 This directory contains a single header used to generate tikz drawings of
 boost graphs.

## Requirements
 The [Boost Libraries](http://www.boost.org/) and a compiler that supports the
 C++11 standard. The provided makefiles use [GCC](https://gcc.gnu.org/) by
 default; other compilers may be specified by changing the appropriate makefile
 variable.
