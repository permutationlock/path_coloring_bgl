## Path coloring of graphs
 Implementation of path coloring algorithms for planar and toroid graphs in boost graph library.

## Requirements
 The boost libraries and a compiler that supports the C++11 standard.

## Tests

 * To build the tests
 ```
 make unit_tests
 ```

 * Execute test suite
 ```
 ./unit_tests
 ```
 
 * Cleaning build
 ```
 make clean
 ```

## To-do list
 
 * Write algorithm for embedding toroidal graphs and locating a non-contractible chordless cycle.
 	This is required to apply the Hartman algorithm to arbitrary toroidal graphs.
 ```
 Input: A planar graph G.
 ```
 ```
 Output: A triangulation of G with vertices arranged in clockwise order for
 	some embedding of G.
 ```
 
 * Write a dot language parser. This will allow us to thoroughly test on arbitrary graphs
 	and provide a nice interface into this algorithm package. Will likely be accomplished with flex and bison
 	and follow the existing dot language grammar.
 	
 * Other algorithms such as list path coloring, and variations on previously stated algorithms will
 	be future work.
 	
 * Consider properties of graphs correctly colored by Poh algorithms. We know planar implies it may be
   3-path-colored by Poh, does being 3-path-colorable by Poh imply planarity (note: almost certainly not)?
   Does being 3-path-colorable by Poh imply anything else?