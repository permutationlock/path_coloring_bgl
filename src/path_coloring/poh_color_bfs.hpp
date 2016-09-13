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

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Local project headers
#include "incidence_list_helpers.hpp"

template<
		typename graph_t, typename planar_embedding_t, typename color_map_t,
		typename mark_map_t, typename parent_map_t, typename color_t,
		typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor
	>
static void poh_color_bfs_recursive(
		const graph_t & graph, const planar_embedding_t & planar_embedding,
		color_map_t & color_map, mark_map_t & mark_map,
		parent_map_t & parent_map, vertex_t p_0, vertex_t p_n, vertex_t q_0,
		vertex_t q_m, std::size_t & count, std::size_t p_mark,
		std::size_t q_mark, color_t new_color
	)
{
	vertex_t t_0 = p_0, t_1 = p_n;
	
	do {
		if(p_0 == p_n && q_0 == q_m) return;
		
		auto edge_iter = find_neighbor_iterator(
				q_m, p_n, planar_embedding, graph
			);
			
		if(edge_iter == planar_embedding[q_m].end()) {
			throw std::runtime_error("No edge between p_n and q_m).");
		}
		else if(++edge_iter == planar_embedding[q_m].end()) {
			t_1 = get_incident_vertex(
					q_m, *planar_embedding[q_m].begin(), graph
				);
		}
		else {
			t_1 = get_incident_vertex(q_m, *edge_iter, graph);
		}
		
		if(mark_map[t_1] == p_mark) {
			p_n = t_1;
		}
		else if(mark_map[t_1] == q_mark) {
			q_m = t_1;
		}
	} while(mark_map[t_1] == p_mark || mark_map[t_1] == q_mark);
	
	vertex_t current_vertex;
	
	{
		std::queue<vertex_t> bfs_queue;
		std::size_t bfs_mark = count++;
		bfs_queue.push(t_1);
		mark_map[t_1] = bfs_mark;
		parent_map[t_1] = t_1;
	
		while(t_0 == p_0 && !bfs_queue.empty()) {
			current_vertex = bfs_queue.front();
			bfs_queue.pop();
		
			auto edge_iter = planar_embedding[current_vertex].begin();
			vertex_t last_neighbor = get_incident_vertex(
					current_vertex, *edge_iter, graph
				);
			
			do {
				if(++edge_iter == planar_embedding[current_vertex].end())
					edge_iter = planar_embedding[current_vertex].begin();
				
				vertex_t neighbor = get_incident_vertex(
						current_vertex, *edge_iter, graph
					);
				
				std::size_t mark = mark_map[neighbor];
				std::size_t last_mark = mark_map[last_neighbor];
				
				if(mark != bfs_mark && mark != p_mark && mark != q_mark) {
					parent_map[neighbor] = current_vertex;
					mark_map[neighbor] = bfs_mark;
					bfs_queue.push(neighbor);
				}
				else if(mark == q_mark && last_mark == p_mark) {
					t_0 = current_vertex;
					vertex_t p_i = last_neighbor;
					vertex_t q_j = neighbor;
					
					if(p_i != p_0 || q_j != q_0) {
						poh_color_bfs_recursive(
								graph, planar_embedding, color_map, mark_map,
								parent_map, p_0, p_i, q_0, q_j,
								count, p_mark, q_mark, new_color
							);
				
						p_0 = p_i;
						q_0 = q_j;
					}
				
					break;
				}
			
				last_neighbor = neighbor;
			} while(edge_iter != planar_embedding[current_vertex].begin());
		}
		
		if(t_0 == p_0) {
			throw std::runtime_error("BFS failed to find edge between paths.");
		}
	}
	
	std::size_t new_mark = count++;
	color_map[current_vertex] = new_color;
	mark_map[current_vertex] = new_mark;
	
	while(parent_map[current_vertex] != current_vertex) {
		current_vertex = parent_map[current_vertex];
		color_map[current_vertex] = new_color;
		mark_map[current_vertex] = new_mark;
	}
	
	poh_color_bfs_recursive(
			graph, planar_embedding, color_map, mark_map, parent_map,
			p_0, p_n, t_0, t_1,
			count, p_mark, new_mark, color_map[q_0]
		);
	
	poh_color_bfs_recursive(
			graph, planar_embedding, color_map, mark_map, parent_map,
			t_0, t_1, q_0, q_m,
			count, new_mark, q_mark, color_map[p_0]
		);
}

template<
		typename graph_t, typename planar_embedding_t,
		typename color_map_t, typename vertex_iterator_t, typename color_t
	>
void poh_color_bfs(
		const graph_t & graph, const planar_embedding_t & planar_embedding,
		color_map_t & color_map, vertex_iterator_t p_begin,
		vertex_iterator_t p_end, vertex_iterator_t q_begin,
		vertex_iterator_t q_end, color_t c_0, color_t c_1, color_t c_2
	)
{
	std::vector<std::size_t> mark_storage(boost::num_vertices(graph));
	typename boost::iterator_property_map<
			std::vector<std::size_t>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> mark_map(
				mark_storage.begin(), boost::get(boost::vertex_index, graph)
			);
	
	typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
	
	std::vector<vertex_t> parent_storage(boost::num_vertices(graph));
	typename boost::iterator_property_map<
			typename std::vector<vertex_t>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> parent_map(
				parent_storage.begin(), boost::get(boost::vertex_index, graph)
			);
	
	for(vertex_iterator_t p_iter = p_begin; p_iter != p_end; ++p_iter) {
		mark_map[*p_iter] = 1;
		color_map[*p_iter] = c_0;
	}
	
	for(vertex_iterator_t q_iter = q_begin; q_iter != q_end; ++q_iter) {
		mark_map[*q_iter] = 2;
		color_map[*q_iter] = c_1;
	}
	
	std::size_t count = 3;
	
	poh_color_bfs_recursive(
			graph, planar_embedding, color_map, mark_map, parent_map,
			*p_begin, *(--p_end), *q_begin, *(--q_end),
			count, 1, 2, c_2
		);
}

#endif
