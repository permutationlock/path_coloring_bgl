# Path Coloring

 We implement an algorithm to path *3*-color a triangulated simple plane graph,
 based on a proof by Poh in 1990.

## poh_color_bfs.hpp

 This file implements the Poh algorithm using a breadth first search to locate
 paths to color.

 ```c++
template<
        typename graph_t, typename planar_embedding_t, typename color_map_t,
        typename mark_map_t, typename parent_map_t, typename vertex_iterator_t,
        typename color_t
    >
void poh_color_bfs(
        const graph_t & graph, const planar_embedding_t & planar_embedding,
        mark_map_t & mark_map, parent_map_t & parent_map,
        color_map_t & color_map, vertex_iterator_t p_begin,
        vertex_iterator_t p_end, vertex_iterator_t q_begin,
        vertex_iterator_t q_end, color_t c_0, color_t c_1, color_t c_2
    );
 ```

### Types

 *vertex_t* is *boost::graph_traits<graph_t>::vertex_descriptor*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [EqualityComparable](http://www.sgi.com/tech/stl/EqualityComparable.html), [Assignable](http://www.boost.org/doc/libs/1_64_0/libs/utility/Assignable.html) | None |
 | *planar_embedding_t* | [PlanarEmbedding](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/PlanarEmbedding.html) | None |
 | *color_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_t* |
 | *mark_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be an integer type |
 | *parent_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *vertex_t* |
 | *vertex_iterator_t* | [Input Iterator](http://www.cplusplus.com/reference/iterator/InputIterator/) | *value_type* must be *vertex_t* |

### Input Requirements

 - *graph* is triangulated and has no loops or parallel edges;
 - *planar_embedding* must be a valid planar embedding of *graph*;
 - the colors *c_0*, *c_1*, and *c_2* are distinct;
 - *color_map* must not initially assign any vertex of *graph* the value *c_0*,
   *c_1*, or *c_2*;
 - *mark_map* must initially assign each vertex a value of *0*;
 - the iterator pairs *p_begin*, *p_end* and *q_begin*, *q_end* must each be a
   range of vertices *p_0,..,p_n* and *q_0,...,q_m* in *graph* such that *p_0...p_n*
   and *q_0...q_m* are induced paths in *graph*, and such that *p_0...p_nq_0...q_m*
   is a cycle in *graph*.

### Output

 The *color_map* has been assigned such that it represents a valid path
 *3*-coloring of the vertices of the subgraph bounded by the cycle
 *p_0...p_nq_0...q_m* using the three colors *c_0*, *c_1*, and *c_2*.

### Time Complexity
 
 A call to *poh_color_bfs* on input *graph* with *n* vertices makes at most:
 
 - *O(n^2)* reads and writes to the provided property map structures;
 - *O(n^2)* comparisons between integer, *vertex_t*, and *color_t* type variables;
 - *O(n^2)* calls to *push* and *pop* on a *std::queue<vertex_t>*.
 
 Property map lookups are often constant time, and thus the algorithm often runs
 in *O(n^2)* time. For example: *O(n^2)* time may be achieved using
 [boost::adjacency_list](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/adjacency_list.html)
 for *graph_t* and
 [boost::iterator_property_map](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/iterator_property_map.html)
 for all property maps. See *src/examples/* for full example code.

## poh_color.hpp

 This file implements Poh's algorithm by walking along a colored path to find
 the next path to color.

 ```c++
template<
        typename graph_t, typename planar_embedding_t, typename color_map_t,
        typename mark_map_t, typename neighbor_range_map_t,
        typename vertex_iterator_t, typename color_t
    >
void poh_color(
        const graph_t & graph, const planar_embedding_t & planar_embedding,
        color_map_t & color_map, neighbor_range_map_t & neighbor_range_map,
        mark_map_t & mark_map, vertex_iterator_t p_begin,
        vertex_iterator_t p_end, vertex_iterator_t q_begin,
        vertex_iterator_t q_end, color_t c_0, color_t c_1, color_t c_2
    );
 ```

### Types

 *vertex_t* is *boost::graph_traits<graph_t>::vertex_descriptor*
 *edge_iterator_t* is *boost::property_traits<embedding_t>::value_type::const_iterator*
 *neighbor_range_t* is *std::pair<edge_iterator_t,edge_iterator_t>*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [EqualityComparable](http://www.sgi.com/tech/stl/EqualityComparable.html), [Assignable](http://www.boost.org/doc/libs/1_64_0/libs/utility/Assignable.html) | None |
 | *planar_embedding_t* | [PlanarEmbedding](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/PlanarEmbedding.html) | None |
 | *color_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_t* |
 | *mark_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be an integer type |
 | *neighbor_range_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *neighbor_range_t* |
 | *vertex_iterator_t* | [Input Iterator](http://www.cplusplus.com/reference/iterator/InputIterator/) | *value_type* must be *vertex_t* |

### Input Requirements

 - *graph* is triangulated and has no loops or parallel edges;
 - *planar_embedding* must be a valid planar embedding of *graph*;
 - the colors *c_0*, *c_1*, and *c_2* are distinct;
 - *color_map* must not initially assign any vertex of *graph* the value *c_0*,
   *c_1*, or *c_2*;
 - *mark_map* must initially assign each vertex a value of *0*;
 - the iterator pairs *p_begin*, *p_end* and *q_begin*, *q_end* must each be a
   range of vertices *p_0,..,p_n* and *q_0,...,q_m* in *graph* such that *p_0...p_n*
   and *q_0...q_m* are induced paths in *graph*, and such that *p_0...p_nq_0...q_m*
   is a cycle in *graph*.

### Output

 The *color_map* has been assigned such that it represents a valid path
 *3*-coloring of the vertices of subgraph of *graph* bounded by the cycle
 *p_0...p_nq_0...q_m* using the three colors *c_0*, *c_1*, and *c_2*.

### Time Complexity
 
 A call to *poh_color* on input *graph* with *n* vertices makes at most:
 
 - *O(n)* reads and writes to the provided property map structures;
 - *O(n)* comparisons between integer, *vertex_t*, and *color_t* type variables;
 
 Property map lookups are often constant time, and thus the algorithm often runs
 in *O(n)* time. For example: *O(n)* time may be achieved using
 [boost::adjacency_list](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/adjacency_list.html)
 for *graph_t* and
 [boost::iterator_property_map](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/iterator_property_map.html)
 for all property maps. See *src/examples/* for full example code.

# Path List-Coloring

We implement an algorithm to path list-color a plane graph given color lists of
size *3* or more for each vertex, based on proofs by Hartman and Skrekovski.

## hartman_skrekovski_color.hpp

 This file implements the Hartman-Skrekovski algorithm.

 ```c++
template<
        typename graph_t, typename planar_embedding_t, typename color_map_t,
        typename mark_map_t, typename neighbor_range_map_t,
        typename vertex_iterator_t, typename color_t
    >
void poh_color(
        const graph_t & graph, const planar_embedding_t & planar_embedding,
        color_map_t & color_map, neighbor_range_map_t & neighbor_range_map,
        mark_map_t & mark_map, vertex_iterator_t p_begin,
        vertex_iterator_t p_end, vertex_iterator_t q_begin,
        vertex_iterator_t q_end, color_t c_0, color_t c_1, color_t c_2
    );
 ```

### Types

 *vertex_t* is *boost::graph_traits<graph_t>::vertex_descriptor*
 
 *edge_iterator_t* is *boost::property_traits<embedding_t>::value_type::const_iterator*
 
 *neighbor_range_t* is *std::pair<edge_iterator_t,edge_iterator_t>*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [EqualityComparable](http://www.sgi.com/tech/stl/EqualityComparable.html), [Assignable](http://www.boost.org/doc/libs/1_64_0/libs/utility/Assignable.html) | None |
 | *planar_embedding_t* | [PlanarEmbedding](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/PlanarEmbedding.html) | None |
 | *color_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_t* |
 | *mark_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be an integer type |
 | *neighbor_range_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *neighbor_range_t* |
 | *vertex_iterator_t* | [Input Iterator](http://www.cplusplus.com/reference/iterator/InputIterator/) | *value_type* must be *vertex_t* |

### Input Requirements

 - *graph* is triangulated and has no loops or parallel edges;
 - *planar_embedding* must be a valid planar embedding of *graph*;
 - the colors *c_0*, *c_1*, and *c_2* are distinct;
 - *color_map* must not initially assign any vertex of *graph* the value *c_0*,
   *c_1*, or *c_2*;
 - *mark_map* must initially assign each vertex a value of *0*;
 - the iterator pairs *p_begin*, *p_end* and *q_begin*, *q_end* must each be a
   range of vertices *p_0,..,p_n* and *q_0,...,q_m* in *graph* such that *p_0...p_n*
   and *q_0...q_m* are induced paths in *graph*, and such that *p_0...p_nq_0...q_m*
   is a cycle in *graph*.

### Output

 The *color_map* has been assigned such that it represents a valid path
 3-coloring of the vertices of the subgraph bounded by the cycle
 *p_0...p_nq_0...q_m* using the three colors *c_0*, *c_1*, and *c_2*.

### Time Complexity
 
 A call to *poh_color* on input *graph* with *n* vertices makes at most:
 
 - *O(n)* reads and writes to the provided property map structures;
 - *O(n)* comparisons between integer, *vertex_t*, and *color_t* type variables;
 
 Property map lookups are often constant time, and thus the algorithm often runs
 in *O(n)* time. For example: *O(n)* time may be achieved by using
 [boost::adjacency_list](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/adjacency_list.html)
 for *graph_t* and
 [boost::iterator_property_map](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/iterator_property_map.html)
 for all property maps. See *src/examples/* for full example code.
 
