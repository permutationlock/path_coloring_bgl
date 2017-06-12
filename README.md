## Path Coloring Algorithms on Plane Graphs
 An implementation of three algorithms for path coloring and path list coloring
 plane graphs.

## Abstract
 A path coloring of a graph partitions its vertex set into color classes such
 that each class induces a disjoint union of paths. In this project we consider
 implementing several algorithms to compute path colorings of graphs embedded in
 the plane.

 We present two implementations of the algorithm to path $3$-color plane graphs
 based on a proof by Poh's in 1990. First we describe a naive algorithm that
 directly follows Poh's procedure, and then we give a modified algorithm
 that runs in linear time.

 The independent results of Hartman and Skrekovski describe a procedure that takes
 a plane graph $G$, and a list of $3$ colors for each vertex in $G$, and
 computes a path coloring of $G$ such that each vertex recieves a color from its
 list. We provide a linear time algorithm based on Hartman and Skrekovski's
 proofs.

 A C++ implementation is provided for all three algorithms, utilizing
 the boost graph library.

## Requirements
 The boost libraries and a compiler that supports the C++11 standard. The
 provided makefile uses GCC by default.

## Tests

 * To build the tests from unit_tests directory
  ```
 make
 ```

 * Execute test suite
 ```
 ./unit_tests
 ```
 
 * Cleaning build
 ```
 make clean
 ```
