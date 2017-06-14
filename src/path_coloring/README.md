## Path 3-Coloring Plane Graphs

 These files implement an algorithm for path coloring a plane graph based on a
 proof by Poh in 1990.

### poh_color_bfs.hpp

 This file implements Poh's algorithm using breadth first search to locate paths
 to color.
 

 ```c++
 template<
        typename graph_t, typename planar_embedding_t, typename color_map_t,
        typename vertex_iterator_t, typename color_t
    >
 void poh_color_bfs(
         const graph_t & graph, const planar_embedding_t & planar_embedding,
         color_map_t & color_map, vertex_iterator_t p_begin,
         vertex_iterator_t p_end, vertex_iterator_t q_begin,
         vertex_iterator_t q_end, color_t c_0, color_t c_1, color_t c_2
     )
 ```

#### Types

 *vertex_t* is *boost::graph_traits<graph_t>::vertex_descriptor*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [Assignable](http://www.sgi.com/tech/stl/EqualityComparable.html), [Assignable](http://www.boost.org/doc/libs/1_64_0/libs/utility/Assignable.html) | None |
 | *planar_embedding_t* | [PlanarEmbedding](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/PlanarEmbedding.html) | None |
 | *color_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_t* |
 | *mark_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be an integer type |
 | *parent_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *vertex_t* |

#### Input Requirements

 - *planar_embedding* must be a valid planar embedding of *graph*
 - the colors *c_0*, *c_1*, and *c_2* are distinct
 - *color_map* must not initially assign any vertex of *graph* any of *c_0*,
   *c_1*, or *c_2*
 - *mark_map* must initially assign each vertex a value of *0*
 - the iterator pairs *p_begin*, *p_end* and *q_begin*, *q_end* must be for a
   range of vertices *p_0,..,p_n* and *q_0,...,q_m* in *graph* such that *p_0...p_n*
   and *q_0...q_m* are induced paths in *graph* and such that *p_0...p_nq_0...q_m*
   is a cycle in *graph*

#### Output

 The *color_map* has been assigned such that it represents a valid path
 3-coloring of the vertices of subgraph of *graph* bounded by the cycle
 *p_0...p_nq_0...q_m* using the three colors *c_0*, *c_1*, and *c_2*.

#### Time Complexity
 
 We take the input size to be the number of vertices *n* in *graph*. A call to
 *poh_path_color* makes at most *O(n^2)* reads and writes to the provided property map
 structures.
 
 Property map lookups are constant time given a structure such as
 [boost::adjacency_list](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/adjacency_list.html).
 
