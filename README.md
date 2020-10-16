# Path Coloring Algorithms for Plane Graphs
 A C++ implementation in the [Boost Graph Library](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/index.html) of two algorithms for path
 coloring and path list-coloring plane graphs. The main parts of the project are
 the [header library and documentation](https://github.com/permutationlock/path_coloring_bgl/tree/master/src/include/path_coloring),
 a [research paper](https://github.com/permutationlock/path_coloring_bgl/blob/master/doc/writeup/final_paper.pdf),
 and [presentation slides](https://github.com/permutationlock/path_coloring_bgl/blob/master/doc/slides/presentation.pdf)
 from a talk given on the research.

## Abstract
 A path coloring of a graph partitions its vertex set into color classes such
 that each class induces a disjoint union of paths. In slightly less technical
 terms, a set of vertices S in a graph G induces a disjoint union of paths if
 the subgraph of G formed by the vertices in S and
 all the edges between vertices in S is a collection of (disjoint) paths. In this project we
 implement several algorithms to compute path colorings of graphs embedded in
 the plane.

 We present two algorithms to path color plane graphs with *3* colors, based on
 a proof by Poh in 1990. First we describe a naive algorithm that directly
 follows Poh's procedure, then we give a modified algorithm that runs in linear
 time.

 Independent results of Hartman and Skrekovski describe a procedure that takes a
 plane graph and a list of *3* colors for each vertex, and computes a path
 coloring such that each vertex recieves a color from its list. We provide
 a linear time algorithm based on Hartman and Skrekovski's proofs.

 A C++ implementation is provided for all three algorithms, utilizing the Boost
 Graph Library. Instructions are given on how to use the implementation
 to construct colorings for plane graphs represented by Boost data
 structures.
