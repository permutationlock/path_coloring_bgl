## Path coloring of graphs
 Implementation of path coloring algorithms for planar and toroid graphs.

## Requirements
 A compiler that supports the C++11 standard.

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

 * Write Poh algorithm for planar graph path coloring. This algorithm is required as a base for all future work.
 ```
 Input: Triangulated plane graph G with vertices arranged in clockwise order
 	(other embedding details not required), a cycle C of G 2-colored such that each
 	color class induces a single cycle.
 ```
 ```
 Output: A coloring of G such that the interior face of cycle C is 3-colored such
 	that each color class induces a disjoint union of paths.
 ```
 
 * Write planar graph embedding algorithm. This algorithm is required
 for toroidal graph coloring.
 ```
 Input: A planar graph G.
 ```
 ```
 Output: The graph G with vertices arranged in clockwise order for
 	some embedding in the plane.
 ```
 
 * Write planar graph triangulation algorithm. This algorithm is required
 for toroidal graph coloring. Combined with the embedding algorithm this allows us to apply
 the Poh algorithm to arbitrary planar graphs.
 ```
 Input: A planar graph G with vertices arranged in clockwise order for some embedding in the plane.
 ```
 ```
 Output: A weakly triangulated version of G with vertices still in clockwise order.
 ```
 
 * Write Hartman algorithm for toroidal graph path coloring. This algorithms builds upon Poh and
 	is one of the main goals of the project.
 ```
 Input: A toroidal graph G with vertices arranged in clockwise order for some embedding,
 	a non-contractible chordless cycle C of G.
 ```
 ```
 Output: A 3-coloring of G such that two color classes induce a disjoint union of paths,
 	and the third color class induces a disjoint union of paths and the cycle C.
 ```
 
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