/*
 * poh_color_bfs.hpp
 * Author: Aven Bross
 * 
 * Implementation of Poh path 3-coloring algorithm for triangulated plane graphs
 * that uses a breadth first search to find chordless paths.
 */

#ifndef __POH_COLOR_BFS_HPP
#define __POH_COLOR_BFS_HPP

// STL headers
#include <vector>
#include <queue>
#include <stdexcept>
#include <utility>
#include <string>

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Local project headers
#include "incidence_list_helpers.hpp"


/* poh_color_bfs_recursive
 * 
 * inputs: A weakly triangulated planar graph with vertex indices (predfined
 *     boost property), a valid planar embedding of the graph modeling the boost
 *     PlanarEmbedding concept, a read-write-able vertex property map to store
 *     vertex colors, a read-write-able vertex property map to store integer
 *     marks for BFS, a read-write-able vertex property map to store the parent
 *     vertex for backtracking after the BFS, the first and last vertex of two
 *     disjoint colored paths P=p_0...p_n and Q=q_0...q_m such that p_0...p_nq_m
 *     ...q_0 is a cycle, an unsigned int count such that all vertex marks
 *     assigned so far are less than count, and finally the third color that has
 *     not been used on the paths P and Q.
 *
 * output: The coloring vertex property will contain a valid path 3-coloring of
 *     the subgraph bounded by p_0...p_nq_m...q_0 such that no interior vertex
 *     shares a color with a neighbor in P or Q.
 */

namespace {
    template<
            typename graph_t,
            typename planar_embedding_t,
            typename color_map_t,
            typename mark_map_t,
            typename parent_map_t,
            typename color_t,
            typename vertex_t
                = typename boost::graph_traits<graph_t>::vertex_descriptor
        >
    void poh_color_bfs_recursive(
            const graph_t & graph,
            const planar_embedding_t & planar_embedding,
            color_map_t & color_map,
            mark_map_t & mark_map,
            parent_map_t & parent_map,
            vertex_t p_0, vertex_t p_n, vertex_t q_0,
            vertex_t q_m, std::size_t count, color_t new_color
        )
    {
        color_t p_color = color_map[p_0], q_color = color_map[q_0];
        vertex_t t_0 = p_0, t_l = p_n;
    
        // Remove triangles from the end until we find an interior vertex t_l
        do {
            if(p_0 == p_n && q_0 == q_m) return;
        
            // Find the edge between p_n and q_m in q_m's incidence list
            auto edge_iter = find_edge_iterator(
                    q_m, p_n, planar_embedding, graph
                );
        
            // Find the neighbor t_l counterclockwise from p_n around q_m
            if(edge_iter == planar_embedding[q_m].end()) {
                throw std::runtime_error("No edge between p_n and q_m).");
            }
            else if(++edge_iter == planar_embedding[q_m].end()) {
                t_l = get_incident_vertex(
                        q_m, *planar_embedding[q_m].begin(), graph
                    );
            }
            else {
                t_l = get_incident_vertex(q_m, *edge_iter, graph);
            }
        
            // If t_l is in P or Q we have found a colored face and remove it
            if(color_map[t_l] == p_color) {
                p_n = t_l;
            }
            else if(color_map[t_l] == q_color) {
                q_m = t_l;
            }
        } while(color_map[t_l] == p_color || color_map[t_l] == q_color);
    
        vertex_t current_vertex, p_i = p_0, q_j = q_0;
    
        // Perform a BFS from t_l to find and color a path T between P and Q
        {
            std::queue<vertex_t> bfs_queue;
            std::size_t bfs_mark = count++;
            bfs_queue.push(t_l);
            mark_map[t_l] = bfs_mark;
            parent_map[t_l] = t_l;
    
            while(t_0 == p_0 && !bfs_queue.empty()) {
                current_vertex = bfs_queue.front();
                bfs_queue.pop();
        
                auto edge_iter = planar_embedding[current_vertex].begin();
                vertex_t last_neighbor = get_incident_vertex(
                        current_vertex, *edge_iter, graph
                    );
            
                // Loop through the neighbors of current_vertex
                do {
                    // If we hit the end of the list, wrap to the start
                    if(++edge_iter == planar_embedding[current_vertex].end())
                        edge_iter = planar_embedding[current_vertex].begin();
                
                    vertex_t neighbor = get_incident_vertex(
                            current_vertex, *edge_iter, graph
                        );
                
                    std::size_t mark = mark_map[neighbor];
                    color_t color = color_map[neighbor];
                    color_t last_color = color_map[last_neighbor];
                
                    // If we hit an unmarked vertex, add it to the queue
                    if(
                            mark != bfs_mark && color != p_color &&
                            color != q_color
                        )
                    {
                        parent_map[neighbor] = current_vertex;
                        mark_map[neighbor] = bfs_mark;
                        bfs_queue.push(neighbor);
                    }
                    // If we find an edge P to Q, we have completed T
                    else if(color == q_color && last_color == p_color) {
                        t_0 = current_vertex;
                        p_i = last_neighbor;
                        q_j = neighbor;
                    
                        break;
                    }
                
                    last_neighbor = neighbor;
                } while(edge_iter != planar_embedding[current_vertex].begin());
            }
        
            if(t_0 == p_0) {
                throw std::runtime_error(
                        "BFS failed to find edge between paths."
                    );
            }
        }
        
        // If the edge p_iq_j is a chord we must color the rest
        if(p_i != p_0 || q_j != q_0) {
            poh_color_bfs_recursive(
                    graph, planar_embedding, color_map,
                    mark_map, parent_map, p_0, p_i, q_0, q_j,
                    count, new_color
                );
        }
    
        color_map[current_vertex] = new_color;
    
        // Backtrack through BFS tree from t_0 to color the path T
        while(parent_map[current_vertex] != current_vertex) {
            current_vertex = parent_map[current_vertex];
            color_map[current_vertex] = new_color;
        }
    
        // Color the subgraph bounded by p_i...p_n and T
        poh_color_bfs_recursive(
                graph, planar_embedding, color_map, mark_map, parent_map,
                p_i, p_n, t_0, t_l, count, q_color
            );
    
        // Color the subgraph bounded by T and q_m...q_j
        poh_color_bfs_recursive(
                graph, planar_embedding, color_map, mark_map, parent_map,
                q_m, q_j, t_l, t_0, count, p_color
            );
    }
}


/* poh_color_bfs
 * 
 * inputs: A weakly triangulated planar graph with vertex indices (predfined
 *     boost property), a valid planar embedding of the graph modeling the boost
 *     PlanarEmbedding concept, a read-write-able vertex property map to store
 *     vertex colors, pairs of bidirectional iterators for lists of vertices
 *     of two disjoint colored paths P=p_0...p_n and Q=q_0...q_m such that
 *     p_0...p_nq_m...q_0 is a cycle, and three colors for the path 3-coloring
 *     (no vertex should be assigned any of these three colors before calling
 *     this function).
 *
 * output: The coloring vertex property will contain a valid path 3-coloring of
 *     the subgraph bounded by p_0...p_nq_m...q_0 such that vertices in P are
 *     recieve the first color, vertices in Q the second color, and no interior
 *     vertex shares a color with a neighbor in P or Q.
 */

template<
        typename graph_t,
        typename planar_embedding_t,
        typename color_map_t,
        typename mark_map_t,
        typename parent_map_t,
        typename vertex_iterator_t,
        typename color_t
    >
void poh_color_bfs(
        const graph_t & graph,
        const planar_embedding_t & planar_embedding,
        color_map_t & color_map,
        mark_map_t & mark_map,
        parent_map_t & parent_map,
        vertex_iterator_t p_begin, vertex_iterator_t p_end,
        vertex_iterator_t q_begin, vertex_iterator_t q_end,
        color_t c_0, color_t c_1, color_t c_2
    )
{ 
    // Color the path P
    for(vertex_iterator_t p_iter = p_begin; p_iter != p_end; ++p_iter) {
        color_map[*p_iter] = c_0;
    }
    
    // Color the path Q
    for(vertex_iterator_t q_iter = q_begin; q_iter != q_end; ++q_iter) {
        color_map[*q_iter] = c_1;
    }
    
    // Construct the path 3-coloring
    poh_color_bfs_recursive(
            graph, planar_embedding, color_map, mark_map, parent_map,
            *p_begin, *(--p_end), *q_begin, *(--q_end), 1, c_2
        );
}

#endif
