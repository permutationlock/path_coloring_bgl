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
