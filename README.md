# Path Coloring Algorithms on Plane Graphs
 A C++ implementation in the Boost Graph Library of two algorithms for path
 coloring and path list coloring plane graphs.

## Abstract
 A path coloring of a graph partitions its vertex set into color classes such
 that each class induces a disjoint union of paths. In this project we consider
 implementing several algorithms to compute path colorings of graphs embedded in
 the plane.

 We present two algorithms to path *3*-color plane graphs based on a proof by Poh
 in 1990. First we describe a naive algorithm that directly follows Poh's
 procedure, then we give a modified algorithm that runs in linear time.

 Independent results of Hartman and Skrekovski describe a procedure that takes a
 plane graph *G* and a list of *3* colors for each vertex in *G*, and computes a path
 coloring of *G* such that each vertex recieves a color from its list. We provide
 a linear time algorithm based on Hartman and Skrekovski's proofs.

 A C++ implementation is provided for all three algorithms, utilizing the Boost
 Graph Library.

## Requirements
 The boost libraries and a compiler that supports the C++11 standard. The
 provided makefiles use GCC by default; other compilers may be specified by
 changing the appropriate variable.

