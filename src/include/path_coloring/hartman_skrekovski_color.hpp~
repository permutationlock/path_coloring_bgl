/*
 * hartman_skrekovski_color.hpp
 * Author: Aven Bross
 * 
 * Implementation of Hartman-Skrekovski path 3-choosing algorithm.
 */

#ifndef __HARTMAN_SKREKOVSKI_COLOR_HPP
#define __HARTMAN_SKREKOVSKI_COLOR_HPP

// STL headers
#include <vector>
#include <utility>
#include <algorithm>

// Basic graph headers
#include <boost/graph/graph_traits.hpp>

// Local project headers
#include "disjoint_set.hpp"
#include "incidence_list_helpers.hpp"

namespace {
    // The three states a vertex may have
    const int INTERIOR_MARK = 0;

    /* 
     * set_face_location
     *
     * inputs: a vertex v, the new face location for v, a read-write-able vertex
     *     property map storing the face location of each vertex, a disjoint set
     *     structure to store face location sets.
     * 
     * outputs: sets the face location of v to the given face location, if the
     *     given location doesn't exist, that is, the integer new_face_location
     *     is negative, a new face location is created in the disjoint set
     *     structure.
     */
    template<typename vertex_t, typename face_location_map_t>
    inline int set_face_location(
            vertex_t v, int new_face_location,
            face_location_map_t & face_location_map,
            disjoint_set_t & face_location_sets
        )
    {
        if(!face_location_sets.exists(new_face_location)) {
            face_location_map[v] = face_location_sets.make_next();
        }
        else {
            face_location_map[v] = face_location_sets.find(new_face_location);
        }
        return face_location_map[v];
    }
    
    /*
     * hartman_skrekovski_color
     *
     * assumptions: There is a cycle C=v_0v_1...v_k in the weakly triangulated
     *     being colored. Suppose x, y, and p are in C such that p is between x
     *     and y clockwise around C. Suppose all vertices other than x, y, and p
     *     are uncolored. Suppose vertices x, y, and p have lists of size 1 or
     *     more, other vertices in C have lists of size 2 or more, and vertices
     *     in the subgraph bounded by C have lists of size 3 or more.
     *     Additionally, if p is colored some color c, assume vertices between x
     *     and p clockwise, including x, don't have the color c in their lists.
     * 
     * inputs: A weakly triangulated planar graph with vertex indices
     *     (predefined boost property), a valid augmented planar embedding of
     *     the graph (modeling the boost AugmentedEmbedding concept defined in
     *     documentation), a read-write-able vertex property map assigning each
     *     vertex an integer to track face location, a disjoint set structure to
     *     store face location sets, a read-write-able vertex property map
     *     assigning each vertex an integer to track state, a read-able vertex
     *     property map assigning a pair of iterators from its augmented
     *     embedding list to track the range of valid neighbors, a read-write-
     *     able vertex property map to which the coloring will be assigned, the
     *     vertices x, y, and p on the outer cycle of the current subgraph, and
     *     the integers for the current face location marks.
     * 
     * outputs: The coloring vertex property will contain a valid path list-
     *     coloring of the subgraph bounded by C such that x, y, and p each
     *     receive at most one same color neighbor. Moreover, if x=y=p then
     *     x,y,p receive no same color neighbors.
     */

    template<
            typename graph_t, typename augmented_embedding_t,
            typename face_location_map_t, typename neighbor_range_map_t,
            typename color_list_map_t, typename vertex_t
                = typename boost::graph_traits<graph_t>::vertex_descriptor
        >
    void hartman_skrekovski_color_recursive(
            const graph_t & graph,
            const augmented_embedding_t & augmented_embedding,
            face_location_map_t & face_location_map,
            disjoint_set_t & face_location_sets, 
            neighbor_range_map_t & neighbor_range_map,
            color_list_map_t & color_list_map,
            vertex_t x, vertex_t y, vertex_t p,
            int before_p, int before_y, int before_x
        )
    {
        // If p isn't colored yet we are starting a new path
        if(color_list_map[p].size() > 1) {
            // Color p the first color in its list
            color_list_map[p].erase(
                    ++color_list_map[p].begin(),
                    color_list_map[p].end()
                );
        }
        
        // Grab p's the color from p's singleton list
        auto p_color = color_list_map[p].front();
        
        // Track the vertices that will become x and y in the first cycle
        vertex_t new_x = x, new_y = y;
        
        // Iterate through p's adjacency list
        auto neighbor_iter = neighbor_range_map[p].first;
        do {
            // Wrap adjacency list
            if(neighbor_iter == augmented_embedding[p].end())
                neighbor_iter = augmented_embedding[p].begin();
            
            // Grab neighbor n, and the iterator to p in n's adjacency list
            vertex_t n = neighbor_iter -> vertex;
            auto back_iter = neighbor_iter -> iterator;
            int n_location = face_location_map[n];
            
            // Look for p_color in n's color list
            auto color_iter = std::find(
                    color_list_map[n].begin(),
                    color_list_map[n].end(), p_color
                );
            
            // The case n is not in C
            if(face_location_map[n] == INTERIOR_MARK) {
                // Note that n is between x and p on the first cycle
                before_p = set_face_location(
                        n, before_p, face_location_map, face_location_sets
                    );
                
                // Initialize n's neighbor range to start at the edge pn
                initialize_neighbor_range(
                        n, back_iter, neighbor_range_map, augmented_embedding
                    );
                
                // Remove the edge pn from n's neighbor range
                remove_first_neighbor(
                        n, neighbor_range_map, augmented_embedding
                    );
                
                // Remove the path color from n's list
                if(color_iter != color_list_map[n].end()) {
                    color_list_map[n].erase(color_iter);
                    color_iter = color_list_map[n].end();
                }
            }
            // The case n is immediately counterclockwise to p in C
            else if(neighbor_iter == neighbor_range_map[p].first) {
                // If p has a single neighbor the current subgraph is K_2
                if(neighbor_iter == neighbor_range_map[p].second) {
                    // If is neither x nor y, ensure it doesn't receive p_color
                    if(n != x && n != y) {
                        // Remove the path color from n's list
                        if(color_iter != color_list_map[n].end()) {
                            color_list_map[n].erase(color_iter);
                            color_iter = color_list_map[n].end();
                        }
                    }
                    
                    // Color n with any color remaining in its list
                    if(color_list_map[n].size() > 1) {
                        color_list_map[n].erase(
                                ++color_list_map[n].begin(),
                                color_list_map[n].end()
                            );
                    }
                    
                    break;
                }
                // If n is y and y has the path color, swap p and y so the
                //     path runs clockwise
                else if(n == y && color_iter != color_list_map[n].end()) {
                    if(color_list_map[n].size() > 1) {
                        // Remove all colors other than the path color
                        color_list_map[n].erase(
                                color_list_map[n].begin(),
                                color_iter
                            );
                        color_list_map[n].erase(
                                ++color_list_map[n].begin(),
                                color_list_map[n].end()
                            );
                    }
                    
                    hartman_skrekovski_color_recursive(
                            graph, augmented_embedding,
                            face_location_map, face_location_sets,
                            neighbor_range_map, color_list_map,
                            y, p, y,
                            -1, -1, before_y
                        );
                    
                    break;
                }
                else {
                    // If we are removing x, setup new_x as n
                    if(x == p) {
                        new_x = n;
                        
                        // Note that n is between x and p on the first cycle
                        before_p = set_face_location(
                                n, before_p, face_location_map,
                                face_location_sets
                            );
                    }
                    
                    // Remove path color from n's list
                    if(color_iter != color_list_map[n].end()) {
                        color_list_map[n].erase(color_iter);
                        color_iter = color_list_map[n].end();
                    }
                    
                    // Remove the edge pn from n's adjacency list
                    remove_last_neighbor(
                            n, neighbor_range_map, augmented_embedding
                        );
                }
            }
            else {
                // Divide the neighbor range of n at the edge pn
                auto n_ranges = split_neighbor_range(
                        n, back_iter, neighbor_range_map, augmented_embedding
                    );
                
                // p has been removed from the first cycle so we may remove
                //     those neighbors from p's neighbor range
                neighbor_range_map[p].first = neighbor_iter;
                
                // The case we are removing y
                if(p == y) {
                    if(color_iter != color_list_map[n].end()) {
                        color_list_map[n].erase(color_iter);
                        color_iter = color_list_map[n].end();
                    }
                    
                    // If p=x=y then we are removing a single vertex path
                    if(p == x) {
                        // Color the subgraph bounded by the first cycle
                        neighbor_range_map[n] = n_ranges.second;
                        hartman_skrekovski_color_recursive(
                                graph, augmented_embedding,
                                face_location_map, face_location_sets,
                                neighbor_range_map, color_list_map,
                                new_x, n, new_x,
                                -1, before_p, before_y
                            );
                        
                        // If the edge pn is a chord, we must color the rest
                        if(back_iter != n_ranges.first.first) {
                            // Color the subgraph bounded by the second cycle
                            neighbor_range_map[n] = n_ranges.first;
                            hartman_skrekovski_color_recursive(
                                    graph, augmented_embedding,
                                    face_location_map, face_location_sets,
                                    neighbor_range_map, color_list_map,
                                    n, p, p,
                                    -1, -1, before_y
                                );
                        }
                        
                        break;
                    }
                    // The case the edge pn is in the cycle, i.e. not a chord
                    else if(neighbor_iter == neighbor_range_map[p].second) {
                        // Color the subgraph bounded by the first cycle
                        neighbor_range_map[n] = n_ranges.second;
                        hartman_skrekovski_color_recursive(
                                graph, augmented_embedding,
                                face_location_map, face_location_sets,
                                neighbor_range_map, color_list_map,
                                new_x, n, new_x,
                                -1, before_p, before_x
                            );
                        
                        break;
                    }
                    // Otherwise we condition on the location of n as usual
                    
                    // Remember n will be our new y vertex for the first cycle
                    new_y = n;
                }
                
                // The case n is in C[p,y] 
                if(face_location_sets.compare(n_location, before_y) || n == y) {
                    vertex_t new_p = n;
                    
                    // The case n needs to be colored
                    if(color_iter != color_list_map[n].end()
                        && color_list_map[n].size() > 1)
                    {
                        // Remove all colors other than the path color
                        color_list_map[n].erase(
                                color_list_map[n].begin(),
                                color_iter
                            );
                        color_list_map[n].erase(
                                ++color_list_map[n].begin(),
                                color_list_map[n].end()
                            );
                    }
                    // The case n will not be added to the path
                    else if(color_list_map[n].size() > 1
                        || color_list_map[n].front() != p_color)
                    {
                        // Start new path at the new x on the first cycle
                        new_p = new_x;
                        
                        // Combine the before_p and before_y segments
                        before_y = face_location_sets.take_union(
                                before_p, before_y
                            );
                        
                        // Note before_p is gone
                        before_p = -1;
                    }
                    
                    /*
                     * This is the most complicated case.
                     * There are two sub-cases to consider:
                     * 
                     * Case 1: p_color in L[n]
                     *   In this case |L[n]|=1. In the first call we will
                     *   continue coloring the path from n. In the second
                     *   call, we will immediately find the edge pn and,
                     *   since it the first edge and in the cycle, we will
                     *   immediately remove the two vertex path np. Thus
                     *   no vertex in the path will recieve a new same color
                     *   neighbor from the second call.
                     *
                     * Case 2: p_color is not in L[n]
                     *   In this case the first call will start a new path
                     *   from the vertex new_x. By removing p_color from all
                     *   neighbors of vertices in the path, we ensure no
                     *   vertices in the path, including n, will recieve any
                     *   new same color neighbors in the first call. In the
                     *   second call we will be continuing the path from p.
                     *   Since p_color is not in L[n] and x=y=n, the vertex
                     *   n will not recieve any new same color neighbors in
                     *   the second call.
                     
                     */
                    
                    // Color the subgraph bounded by the first cycle
                    neighbor_range_map[n] = n_ranges.second;
                    hartman_skrekovski_color_recursive(
                            graph, augmented_embedding,
                            face_location_map, face_location_sets,
                            neighbor_range_map, color_list_map,
                            new_x, new_y, new_p,
                            before_p, before_y, before_x
                        );
                    
                    // If the edge pn is a chord, handle the rest
                    if(back_iter != n_ranges.first.first) {
                        // Color the subgraph bounded by the second cycle
                        neighbor_range_map[n] = n_ranges.first;
                        hartman_skrekovski_color_recursive(
                                graph, augmented_embedding,
                                face_location_map, face_location_sets,
                                neighbor_range_map, color_list_map,
                                n, n, p,
                                -1, before_y, -1
                            );
                    }
                }
                // The case n is in C[y,x], n != y
                else if(face_location_sets.compare(n_location, before_x)) {
                    if(color_iter != color_list_map[n].end()) {
                        color_list_map[n].erase(color_iter);
                        color_iter = color_list_map[n].end();
                    }
                    
                    // Color the subgraph bounded by the first cycle
                    neighbor_range_map[n] = n_ranges.second;
                    hartman_skrekovski_color_recursive(
                            graph, augmented_embedding,
                            face_location_map, face_location_sets,
                            neighbor_range_map, color_list_map,
                            new_x, n, new_x,
                            -1, before_p, before_x
                        );
                    
                    // Color the subgraph bounded by the second cycle
                    neighbor_range_map[n] = n_ranges.first;
                    hartman_skrekovski_color_recursive(
                            graph, augmented_embedding,
                            face_location_map, face_location_sets,
                            neighbor_range_map, color_list_map,
                            n, y, p,
                            -1, before_y, before_x
                        );
                }
                // The case n is in C[x,p]
                else {
                    // Color the subgraph bounded b the second cycle
                    neighbor_range_map[n] = n_ranges.first;
                    hartman_skrekovski_color_recursive(
                            graph, augmented_embedding,
                            face_location_map, face_location_sets,
                            neighbor_range_map, color_list_map,
                            x, y, p,
                            before_p, before_y, before_x
                        );
                    
                    // Color the subgraph bounded by the first cycle
                    neighbor_range_map[n] = n_ranges.second;
                    hartman_skrekovski_color_recursive(
                            graph, augmented_embedding,
                            face_location_map, face_location_sets,
                            neighbor_range_map, color_list_map,
                            n, n, n,
                            -1, before_p, -1
                        );
                }
                
                break;
            }
        } while(neighbor_iter++ != neighbor_range_map[p].second);
    }
}


/*
 * hartman_skrekovski_color
 * 
 * inputs: A weakly triangulated planar graph with vertex indices (predefined
 *     boost property), a valid augmented planar embedding of the graph
 *     (modeling the boost AugmentedEmbedding concept defined in documentation),
 *     a read-able vertex property map assigning a range of colors to each
 *     vertex (each vertex must recieve a range of at least 3 colors if
 *     interior, and at least two colors if on the outer face), a read-write-
 *     able vertex property map to which the coloring will be assigned, and
 *     finally a pair of iterators for the list of vertices on the outer face of
 *     the graph in clockwise order.
 * 
 * outputs: The coloring will be a valid assignment of colors from the input
 *     color lists such that each color class induces a disjoint union of paths.
 */
 
template<
        typename graph_t,
        typename augmented_embedding_t,
        typename color_list_map_t,
        typename neighbor_range_map_t,
        typename face_location_map_t,
        typename face_iterator_t
    >
void hartman_skrekovski_color(
        const graph_t & graph,
        const augmented_embedding_t & augmented_embedding,
        face_iterator_t face_begin, face_iterator_t face_end,
        neighbor_range_map_t & neighbor_range_map,
        face_location_map_t & face_location_map,
        color_list_map_t & color_list_map
    )
{
    // Type definitions
    typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
    
    // Setup face location sets (we will have only the region before_y)
    disjoint_set_t face_location_sets;
    
    // Make a set for 0, the location given to interior vertices
    face_location_sets.make_next();
    
    // Intitially all vertices will be
    int before_y = -1;
    
    // Initialize vertices on outer face
    for(auto face_iter = face_begin; face_iter != face_end; ++face_iter)
    {
        auto next = face_iter;
        if(++next == face_end) next = face_begin;
        
        vertex_t l = *face_iter;
        vertex_t v = *next;
        
        before_y = set_face_location(
                v, before_y, face_location_map, face_location_sets
            );
        
        auto back_iter = find_neighbor_iterator(
                v, l, augmented_embedding, graph
            );
        initialize_neighbor_range(
                v, back_iter, neighbor_range_map, augmented_embedding
            );
    }
    
    vertex_t x = *face_begin, y = *(--face_end);
    
    // Construct the path choosing from given color lists
    hartman_skrekovski_color_recursive(
            graph, augmented_embedding,
            face_location_map, face_location_sets,
            neighbor_range_map, color_list_map,
            x, y, x,
            -1, before_y, -1
        );
}

template<
        typename graph_t,
        typename augmented_embedding_t,
        typename color_list_map_t,
        typename face_iterator_t
    >
void hartman_skrekovski_color(
        const graph_t & graph,
        const augmented_embedding_t & augmented_embedding,
        face_iterator_t face_begin, face_iterator_t face_end,
        color_list_map_t & color_list_map
    )
{
    // Vertex property map to store vertex marks
    typedef boost::iterator_property_map<
            std::vector<int>::iterator,
            typename boost::property_map<
                    graph_t, boost::vertex_index_t
                >::const_type
        > integer_property_map_t;
    
    // Vertex property map type for the neighbor ranges of planar_embedding_t
    typedef typename boost::property_traits<augmented_embedding_t>::value_type
            ::const_iterator embedding_iterator_t;
    typedef typename std::vector<
            std::pair<embedding_iterator_t, embedding_iterator_t>
        > neighbor_range_storage_t;
    typedef boost::iterator_property_map<
            typename neighbor_range_storage_t::iterator,
            typename boost::property_map<
                    graph_t, boost::vertex_index_t
                >::const_type
        > neighbor_range_map_t;
    
    // Construct a vertex property for neighbor ranges
    neighbor_range_storage_t neighbor_range_storage(boost::num_vertices(graph));
    neighbor_range_map_t neighbor_range_map(
            neighbor_range_storage.begin(),
            boost::get(boost::vertex_index, graph)
        );
    
    // Construct a vertex property map to store face location marks
    std::vector<int> face_location_storage(boost::num_vertices(graph));
    integer_property_map_t face_location_map(
            face_location_storage.begin(),
            boost::get(boost::vertex_index, graph)
        );
    
    // Call Hartman-Skrekovski with the given cycle and 
    hartman_skrekovski_color(
            graph,
            augmented_embedding,
            face_begin, face_end,
            neighbor_range_map,
            face_location_map,
            color_list_map
        );
}

#endif
