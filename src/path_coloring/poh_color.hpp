/*
 * poh_color.hpp
 * Author: Aven Bross
 * 
 * Implementation of Poh path 3-coloring algorithm for triangulated plane graphs.
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

template<
		typename index_graph, typename planar_embedding, typename color_map,
		typename mark_map, typename vertex_map, typename color_type,
		typename vertex_descriptor = typename boost::graph_traits<index_graph>::vertex_descriptor
	>
static void poh_color_recursive(
		const index_graph & graph, const planar_embedding & embedding,
		color_map & coloring, mark_map & vertex_marks, vertex_map & parent_map,
		vertex_descriptor p_0, vertex_descriptor p_1, vertex_descriptor t_0,
		vertex_descriptor q_0, vertex_descriptor q_1,
		std::size_t & count, std::size_t p_mark, std::size_t q_mark, color_type new_color
	)
{
	if(t_0 == p_0) {
		if(p_0 == p_1 && q_0 == q_1) return;
		
		do {
			auto edge_iter = find_neighbor_iterator(p_0, q_0, embedding, graph);
			
			if(edge_iter == embedding[p_0].end()) {
				throw std::runtime_error("Invalid embedding (no t_0).");
			}
			else if(++edge_iter == embedding[p_0].end()) {
				t_0 = get_incident_vertex(p_0, *embedding[p_0].begin(), graph);
			}
			else {
				t_0 = get_incident_vertex(p_0, *edge_iter, graph);
			}
			
			if(vertex_marks[t_0] == p_mark) {
				p_0 = t_0;
			}
			else if(vertex_marks[t_0] == q_mark) {
				q_0 = t_0;
			}
		} while(vertex_marks[t_0] == p_mark || vertex_marks[t_0] == q_mark);
	}
	
	if(p_0 == p_1 && q_0 == q_1) return;
	
	vertex_descriptor t_1;
	
	do {
		auto edge_iter = find_neighbor_iterator(q_1, p_1, embedding, graph);
			
		if(edge_iter == embedding[q_1].end()) {
			throw std::runtime_error("Invalid embedding (no t_1).");
		}
		else if(++edge_iter == embedding[q_1].end()) {
			t_1 = get_incident_vertex(q_1, *embedding[q_1].begin(), graph);
		}
		else {
			t_1 = get_incident_vertex(q_1, *edge_iter, graph);
		}
		
		if(vertex_marks[t_1] == p_mark) {
			p_1 = t_1;
		}
		else if(vertex_marks[t_1] == q_mark) {
			q_1 = t_1;
		}
	} while(vertex_marks[t_1] == p_mark || vertex_marks[t_1] == q_mark);
	
	if(p_0 == p_1 && q_0 == q_1) return;
	
	vertex_descriptor current_vertex;
	
	{
		std::queue<vertex_descriptor> bfs_queue;
		std::size_t bfs_mark = count++;
		bfs_queue.push(t_1);
		vertex_marks[t_1] = bfs_mark;
		parent_map[t_1] = t_1;
	
		while(!bfs_queue.empty()) {
			current_vertex = bfs_queue.front();
			bfs_queue.pop();
		
			if(current_vertex == t_0) break;
		
			auto edge_iter = embedding[current_vertex].begin();
			vertex_descriptor last_neighbor = get_incident_vertex(current_vertex, *edge_iter, graph);
			
			do {
				if(++edge_iter == embedding[current_vertex].end())
					edge_iter = embedding[current_vertex].begin();
				
				vertex_descriptor neighbor = get_incident_vertex(current_vertex, *edge_iter, graph);
			
				std::size_t mark = vertex_marks[neighbor];
				std::size_t last_mark = vertex_marks[last_neighbor];
				if(mark != bfs_mark && mark != p_mark && mark != q_mark) {
					parent_map[neighbor] = current_vertex;
					vertex_marks[neighbor] = bfs_mark;
					bfs_queue.push(neighbor);
				}
				else if(mark == q_mark && last_mark == p_mark) {
					poh_color_recursive(
							graph, embedding, coloring, vertex_marks, parent_map,
							p_0, last_neighbor, t_0, q_0, neighbor,
							count, p_mark, q_mark, new_color
						);
				
					p_0 = last_neighbor;
					q_0 = neighbor;
					t_0 = current_vertex;
				
					break;
				}
			
				last_neighbor = neighbor;
			} while(edge_iter != embedding[current_vertex].begin());
		
			if(current_vertex == t_0) break;
		}
	}
	
	if(current_vertex != t_0) throw std::runtime_error("Invalid embedding (no splitting path).");
	
	std::size_t new_mark = count++;
	coloring[current_vertex] = new_color;
	vertex_marks[current_vertex] = new_mark;
	
	while(parent_map[current_vertex] != current_vertex) {
		current_vertex = parent_map[current_vertex];
		coloring[current_vertex] = new_color;
		vertex_marks[current_vertex] = new_mark;
	}
	
	poh_color_recursive(
			graph, embedding, coloring, vertex_marks, parent_map,
			p_0, p_1, p_0, t_0, t_1,
			count, p_mark, new_mark, coloring[q_0]
		);
	
	poh_color_recursive(
			graph, embedding, coloring, vertex_marks, parent_map,
			t_0, t_1, t_0, q_0, q_1,
			count, new_mark, q_mark, coloring[p_0]
		);
}

template<
		typename index_graph, typename planar_embedding,
		typename color_map, typename vertex_iterator, typename color_type
	>
void poh_color(
		const index_graph & graph, const planar_embedding & embedding,
		color_map & coloring, vertex_iterator p_begin, vertex_iterator p_end,
		vertex_iterator q_begin, vertex_iterator q_end,
		color_type c_0, color_type c_1, color_type c_2
	)
{
	std::vector<std::size_t> mark_storage(boost::num_vertices(graph));
	typename boost::iterator_property_map<
			std::vector<std::size_t>::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> vertex_marks(mark_storage.begin(), boost::get(boost::vertex_index, graph));
	
	typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	
	std::vector<vertex_descriptor> parent_storage(boost::num_vertices(graph));
	typename boost::iterator_property_map<
			typename std::vector<vertex_descriptor>::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> parent_map(parent_storage.begin(), boost::get(boost::vertex_index, graph));
	
	for(vertex_iterator p_iter = p_begin; p_iter != p_end; ++p_iter) {
		vertex_marks[*p_iter] = 1;
		coloring[*p_iter] = c_0;
	}
	
	for(vertex_iterator q_iter = q_begin; q_iter != q_end; ++q_iter) {
		vertex_marks[*q_iter] = 2;
		coloring[*q_iter] = c_1;
	}
	
	std::size_t count = 3;
	
	poh_color_recursive(
			graph, embedding, coloring, vertex_marks, parent_map,
			*p_begin, *(--p_end), *p_begin, *q_begin, *(--q_end),
			count, 1, 2, c_2
		);
}

#endif
