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

### Type definitions

 | Type | Definition |
 | ---- | --- |
 | *vertex_t* | *boost::graph_traits<graph_t>::vertex_descriptor*
 
### Template requirements

 *key_type* for all property maps must be *vertex_t*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [EqualityComparable](http://en.cppreference.com/w/cpp/concept/EqualityComparable), [CopyAssignable](http://en.cppreference.com/w/cpp/concept/CopyAssignable) | None |
 | *planar_embedding_t* | [PlanarEmbedding](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/PlanarEmbedding.html) | None |
 | *color_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_t* |
 | *mark_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be an integer type |
 | *parent_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *vertex_t* |
 | *vertex_iterator_t* | [Input Iterator](http://en.cppreference.com/w/cpp/concept/InputIterator) | *value_type* must be *vertex_t* |

### Input Requirements

 - *graph* is triangulated and has no loops or parallel edges;
 - *planar_embedding* must be ordered as a valid planar embedding of *graph*;
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

### Type Definitions

 | Type | Definition |
 | ---- | --- |
 | *vertex_t* | *boost::graph_traits<graph_t>::vertex_descriptor*
 | *neighbor_iterator_t* | *boost::property_traits<planar_embedding_t>::value_type::const_iterator*
 | *neighbor_range_t* | *std::pair<neighbor_iterator_t,neighbor_iterator_t>*

### Template Requirements

 *key_type* for all property maps must be *vertex_t*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [EqualityComparable](http://en.cppreference.com/w/cpp/concept/EqualityComparable), [CopyAssignable](http://en.cppreference.com/w/cpp/concept/CopyAssignable) | None |
 | *planar_embedding_t* | [PlanarEmbedding](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/PlanarEmbedding.html) | None |
 | *color_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_t* |
 | *mark_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be an integer type |
 | *neighbor_range_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *neighbor_range_t* |
 | *vertex_iterator_t* | [Input Iterator](http://en.cppreference.com/w/cpp/concept/InputIterator) | *value_type* must be *vertex_t* |

### Input Requirements

 - *graph* is triangulated and has no loops or parallel edges;
 - *planar_embedding* must be ordered as a valid planar embedding of *graph*;
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
         typename graph_t, typename augmented_embedding_t,
         typename color_list_map_t, typename neighbor_range_map_t,
         typename face_location_map_t, typename face_iterator_t
     >
 void hartman_skrekovski_color(
         const graph_t & graph,
         const augmented_embedding_t & augmented_embedding,
         const color_list_map_t & color_list_map,
         neighbor_range_map_t & neighbor_range_map,
         face_location_map_t & face_location_map,
         face_iterator_t face_begin, face_iterator_t face_end
     );
 ```

### Type Definitions

 | Type | Definition |
 | ---- | --- |
 | *vertex_t* | *boost::graph_traits<graph_t>::vertex_descriptor*|
 | *neighbor_iterator_t* | *boost::property_traits<augmented_embedding_t>::value_type::const_iterator* |
 | *neighbor_range_t* | *std::pair<neighbor_iterator_t,neighbor_iterator_t>*

### Template Requirements

 *key_type* for all property maps must be *vertex_t*
 
 | Type | Concept | Additional Requirements |
 | ---- | ------- | ----- |
 | *graph_t* | [VertexAndEdgeListGraph](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/VertexAndEdgeListGraph.html) | None |
 | *color_t* | [EqualityComparable](http://en.cppreference.com/w/cpp/concept/EqualityComparable), [CopyAssignable](http://en.cppreference.com/w/cpp/concept/CopyAssignable) | None |
 | *augmented_embedding_t* | [AugmentedEmbedding](https://github.com/permutationlock/path_coloring_bgl/tree/master/src/path_coloring#augmented-embeddings) | None |
 | *color_list_t* | [SequenceContainer](http://en.cppreference.com/w/cpp/concept/SequenceContainer) | *value_type* must be *color_t* |
 | *color_list_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *color_list_t* |
 | *face_location_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be a signed integer type |
 | *neighbor_range_map_t* | [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html) | *value_type* must be *neighbor_range_t* |
 | *vertex_iterator_t* | [Input Iterator](http://en.cppreference.com/w/cpp/concept/InputIterator) | *value_type* must be *vertex_t* |

### Input Requirements

 - *graph* is triangulated and has no loops or parallel edges;
 - *augmented_embedding* must be ordered as a valid planar embedding of *graph*;
 - *color_list_map_t* must assign each vertex on the outer cycle a list of 2 or
   more colors, and each vertex interior to the cycle a list of at 3 or more
   colors;
 - *face_location_map* must each initially assign each vertex a
   value of *0*;
 - the iterator pair *face_begin*, *face_end* must be for a range of vertices
   *v_0,..,v_n* such that *v_0...v_n* is a cycle in *graph*.

### Output

 For each vertex in the subgraph bounded by the cycle *v_0...v_n* the
 *color_list_map* has been modified such that all but one color remains in the
 vertex's list. Additionally, if we consider each vertex to be colored with the
 remaining color in its list, the coloring is a path coloring of the subgraph
 bounded by *v_0...v_n*.

### Time Complexity
 
 A call to *hartman_skrekovski_color* on input *graph* with *n* vertices makes
 at most:
 
 - *O(n)* reads and writes to the provided property map structures;
 - *O(n)* comparisons between integer, *vertex_t*, and *color_t* type variables;
 
 Property map lookups are often constant time, and thus the algorithm often runs
 in *O(n)* time. For example: *O(n)* time may be achieved by using
 [boost::adjacency_list](http://www.boost.org/doc/libs/1_64_0/libs/graph/doc/adjacency_list.html)
 for *graph_t* and
 [boost::iterator_property_map](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/iterator_property_map.html)
 for all property maps. See *src/examples/* for full example code.

# Augmented Embeddings

 An augmented embedding assigns each vertex an augmented adjacency list, ordered
 according to a valid planar embedding. An augmented adjacency list differs from
 a standard adjacency list in that given a vertex *v*, the augmented adjacency
 list entry for a neighbor *u* also provides a reference to the entry for *v* in
 *u*'s list.

## AugmentedEmbedding Concept
 
 The AugmentedEmbedding concept refines [LvaluePropertyMap](http://www.boost.org/doc/libs/1_64_0/libs/property_map/doc/LvaluePropertyMap.html),
 placing additional restrictions on the *value_type*.
 
### Types

 | Type | Description |
 | --- | --- |
 | *augmented_embedding_t* | a type that models the [AugmentedEmbedding](https://github.com/permutationlock/path_coloring_bgl/tree/master/src/path_coloring#augmented-embeddings) concept |
 | *node_t* | *boost::property_traits<augmented_embedding_t>::value_type::value_type* |
 | *iterator_t* | *boost::property_traits<augmented_embedding_t>::value_type::iterator* |
 | *graph_t* | the type of the underlying graph |
 | *vertex_t* | *boost::graph_traits<graph_t>::vertex_descriptor* |
 
 | Type | Additional Restrictions |
 | --- | --- |
 | *node_t* | A *node_t n* object must have members *vertex_t n.vertex* and *iterator_t n.iterator*  |
 | *iterator_t* | Must model the [BidirectionalIterator](http://www.cplusplus.com/reference/iterator/BidirectionalIterator/) concept |

### Notation

 | Object | Type |
 | --- | --- |
 | *u*, *v* | objects of the type *vertex_t* |
 | *augmented_embedding* | an object of type *augmented_embedding_t* |
 | *n* | *node_t* |

### Description of Specification

 The object *augmented_embedding* will assign a range of objects of type
 *node_t* to each vertex *v* in the underlying graph. There will be exactly one
 node in this range for each neighbor of *v* in the underlying graph. We will
 call this range of nodes the augmented adjacency list for *v*.
 
 The type *node_t* will represent a neighboring vertex *u* in the augmented
 adjacency list for a vertex *v*. The type *iterator_t* will be an iterator for
 the range of *node_t* objects for a vertex *v*.
 
 For a vertex *v* each node *n* in the range *augmented_embedding[v].begin()*,
 *augmented_embedding[v].end()* will have *n.vertex* be a neighboring vertex *u*
 and *n.iterator* be the unique iterator in the range *augmented_embedding[u].begin()*,
 *augmented_embedding[u].end()* such that *n.iterator->vertex* is equal to *v*.

### Valid Expressions
 
 | Expression | Return Type | Description |
 | --- | --- | --- |
 | *n.vertex | *vertex_t &* | Access the vertex member for a node in the augmented adjacency list |
 | *n.iterator | *iterator_t* | Access the iterator member for a node in the augmented adjacency list |
 | *augmented_embedding[v].begin()* | *iterator_t* | Returns an iterator to the beginning of the range of nodes for the neighbors of vertex *v* |
 | *augmented_embedding[v].end()* | *iterator_t* | Returns an iterator to the end of the range of nodes for the neighbors of the vertex *v* |
 | *augmented_embedding[v].push_back(n)* | *void* | Adds a node *n* to the end of the sequence of nodes for the vertex *v* |

### Example Definition

 ```c++
 typedef boost::adjacency_list<
         boost::setS,
         boost::vecS,
         boost::undirectedS,
         boost::property<boost::vertex_index_t, std::size_t>,
         boost::property<boost::edge_index_t, std::size_t>
     > graph_t;
 
 typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
 
 struct node_t {
         vertex_t vertex;
         typename std::vector<node_t>::iterator iterator;
     };
 
 typedef typename std::vector<std::list<node_t>> embedding_storage_t;
 
 typedef boost::iterator_property_map<
         typename embedding_storage_t::iterator,
         typename boost::property_map<graph_t, vertex_index_t>::const_type
     > augmented_embedding_t;
 ```

 See *src/examples* for full example code.

