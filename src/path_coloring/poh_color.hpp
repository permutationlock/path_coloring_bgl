/*
 * poh_color.hpp
 * Author: Aven Bross
 * 
 * Implementation of Poh path 3-coloring algorithm for triangulated plane graphs
 * that traces a chordless path along the inside of the outer face.
 */

#ifndef __POH_COLOR_HPP
#define __POH_COLOR_HPP

// STL headers
#include <vector>
#include <queue>
#include <stdexcept>
#include <utility>

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Local project headers
#include "incidence_list_helpers.hpp"


/* poh_color_recursive
 * 
 * assumptions: There are two disjoint colored paths P=p_0...p_n and Q=q_0...q_m
 *     such that p_0...p_nq_m...q_0 is a cycle, and all vertices on the interior
 *     of the subgrapph bounded by this cycle with at least one neighbor in P
 *     have been marked with an identifying integer.
 *
 * inputs: A weakly triangulated planar graph with vertex indices (predfined
 *     boost property), a valid planar embedding of the graph modeling the boost
 *     PlanarEmbedding concept, a read-write-able vertex property map to store
 *     vertex colors, a read-write-able vertex property map to store integer
 *     marks for BFS, a read-write-able vertex property map to store start and
 *     stop iterators in each incidence list, the vertex u such that p_0uq_0
 *     is a triangle, the integer mark for vertices with neighbors in P, an
 *     unsigned integer count such that all vertex marks assigned so far are
 *     less than count, and finally the color for the path T to be constructed,
 *     the color of the path P, and the color of the path Q.
 *
 * output: The coloring vertex property will contain a valid path 3-coloring of
 *     the subgraph bounded by p_0...p_nq_m...q_0 such that no interior vertex
 *     shares a color with a neighbor in P or Q.
 */

namespace {
    template<
            typename graph_t, typename planar_embedding_t, typename color_map_t,
            typename mark_map_t, typename neighbor_range_map_t,
            typename color_t, typename vertex_t
                = typename boost::graph_traits<graph_t>::vertex_descriptor
        >
    void poh_color_recursive(
            const graph_t & graph, const planar_embedding_t & planar_embedding,
            color_map_t & color_map, mark_map_t & mark_map,
            neighbor_range_map_t & neighbor_range_map, vertex_t u,
            std::size_t face_mark, std::size_t count,
            color_t t_color, color_t p_color, color_t q_color
        )
    {
        // Let t_i be the last vertex of T, and w the first interior vertex
        vertex_t t_i = u, w = u;
    
        // Initialize remaining variables for the loop constructing T
        auto edge_iter = neighbor_range_map[u].first;
        std::size_t below_t_mark = count++;
        color_map[u] = t_color;
    
        // If t_0 is the only remaining interior vertex we are are done
        if(neighbor_range_map[u].first == neighbor_range_map[u].second){
            return;
        }
    
        // Construct the path T one vertex at a time, with current vertex t_i
        while(edge_iter++ != neighbor_range_map[t_i].second)
        {
            // If we hit the end of the incidence list, wrap to the start
            if(edge_iter == planar_embedding[t_i].end())
                edge_iter = planar_embedding[t_i].begin();
        
            vertex_t n = get_incident_vertex(t_i, *edge_iter, graph);
        
            // Note if the vertex n may be added to T
            bool continue_path = (mark_map[n] == face_mark
                && color_map[n] != p_color && color_map[n] != q_color);
        
            // If we are adding n to T or at the end, color any vertices above T
            if(continue_path || color_map[n] == p_color) {
                if(continue_path) {
                    color_map[n] = t_color;
                }
            
                vertex_t l = n;
            
                // Loop through neighbors of p that lie above T
                while(edge_iter++ != neighbor_range_map[t_i].second) {
                    // If we hit the end of the list, wrap to the start
                    if(edge_iter == planar_embedding[t_i].end())
                        edge_iter = planar_embedding[t_i].begin();
                
                    vertex_t v = get_incident_vertex(t_i, *edge_iter, graph);
                
                    // If v is in P and l is uncolored, make recursive call
                    if(color_map[v] == p_color)
                    {
                        if(color_map[l] != p_color && color_map[l] != q_color
                            && color_map[l] != t_color)
                        {
                            poh_color_recursive(
                                    graph, planar_embedding, color_map,
                                    mark_map, neighbor_range_map, l, face_mark,
                                    count, q_color, p_color, t_color
                                );
                        }
                    }
                
                    l = v;
                }
            
                // If we are adding n to T, setup next loop with t_i = n
                if(continue_path) {
                    auto n_back_iter = find_edge_iterator_restricted(
                            n, t_i, neighbor_range_map[n].first,
                            neighbor_range_map[n].second, planar_embedding,
                            graph
                        );
            
                    neighbor_range_map[n].first = n_back_iter;
                    edge_iter = n_back_iter;
                    t_i = n;
                }
                // Otherwise we have completed T and we are done
                else {
                    return;
                }
            }
            // If n is in Q we have an edge T to Q and may split along it
            else if(color_map[n] == q_color) {
                // Color the left cycle if it has uncolored vertices
                if(w != u) {
                    poh_color_recursive(
                            graph, planar_embedding, color_map, mark_map,
                            neighbor_range_map, w, below_t_mark, count,
                            p_color, t_color, q_color
                        );
                
                    u = t_i;
                    w = t_i;
                }
            }
            // If n is interior and has not been visited, mark and initialize it
            else if(mark_map[n] != below_t_mark) {
                auto back_iter = find_edge_iterator(
                        n, t_i, planar_embedding, graph
                    );
            
                // If this is the first interior vertex hit, save it as w
                if(w == u) {
                    w = n;
                }
            
                // Increment back_iter counterclockwise, wrapping in list 
                if(++back_iter == planar_embedding[n].end()) {
                    back_iter = planar_embedding[n].begin();
                }
            
                initialize_neighbor_range(
                        n, back_iter, neighbor_range_map, planar_embedding
                    );
                mark_map[n] = below_t_mark;
            }
        }
    }
}


/* poh_color
 * 
 * inputs: A weakly triangulated planar graph with vertex indices (predfined
 *     boost property), a valid planar embedding of the graph modeling the boost
 *     PlanarEmbedding concept, a read-write-able vertex property map to store
 *     vertex colors, pairs of bidirectional iterators for lists of vertices
 *     of two disjoint colored paths P=p_0...p_n and Q=q_0...q_m such that
 *     p_0...p_nq_m...q_0 is a cycle, and three colors for the path 3-coloring.
 *
 * output: The coloring vertex property will contain a valid path 3-coloring of
 *     the subgraph bounded by p_0...p_nq_m...q_0 such that vertices in P are
 *     recieve the first color, vertices in Q the second color, and no interior
 *     vertex shares a color with a neighbor in P or Q.
 */

template<
        typename graph_t, typename planar_embedding_t, typename color_map_t,
        typename vertex_iterator_t, typename color_t
    >
void poh_color(
        const graph_t & graph, const planar_embedding_t & planar_embedding,
        color_map_t & color_map, vertex_iterator_t p_begin,
        vertex_iterator_t p_end, vertex_iterator_t q_begin,
        vertex_iterator_t q_end, color_t c_0, color_t c_1, color_t c_2
    )
{
    // Type definitions
    typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
    typedef typename boost::property_traits<planar_embedding_t>::value_type
        ::const_iterator edge_iterator_t;
    
    // Construct a vertex property for marking vertices
    std::vector<std::size_t> mark_storage(boost::num_vertices(graph));
    boost::iterator_property_map<
            std::vector<std::size_t>::iterator,
            typename boost::property_map<graph_t, boost::vertex_index_t>
                ::const_type
        > mark_map(
            mark_storage.begin(), boost::get(boost::vertex_index, graph)
        );
    
    // Construct a vertex property for start/stop positions in incidence lists
    std::vector<std::pair<edge_iterator_t, edge_iterator_t>>
        neighbor_range_storage(boost::num_vertices(graph));
    boost::iterator_property_map<
            typename std::vector<std::pair<edge_iterator_t, edge_iterator_t>>
                ::iterator,
            typename boost::property_map<graph_t, boost::vertex_index_t>
                ::const_type
        > neighbor_range_map(
                neighbor_range_storage.begin(),
                boost::get(boost::vertex_index, graph)
            );
    
    // Intitialize neighbor ranges for vertices in the path P
    vertex_t l = *q_begin;
    for(vertex_iterator_t p_iter = p_begin; p_iter != p_end; ++p_iter) {
        vertex_t v = *p_iter;
        
        auto back_iter = find_edge_iterator(v, l, planar_embedding, graph);
        initialize_neighbor_range(
                v, back_iter, neighbor_range_map, planar_embedding
            );
        
        l = v;
    }
    
    // Color the path Q
    for(vertex_iterator_t q_iter = q_begin; q_iter != q_end; ++q_iter) {
        vertex_t v = *q_iter;
        
        color_map[v] = c_1;
    }
    
    // Construct the path 3-coloring
    poh_color_recursive(
            graph, planar_embedding, color_map, mark_map,
            neighbor_range_map, *p_begin, 1, 2,
            c_0, c_2, c_1
        );
}

#endif
