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

template<
		typename graph_t, typename planar_embedding_t, typename color_map_t,
		typename mark_map_t, typename neighbor_range_map_t, typename color_t,
		typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor
	>
static void poh_color_recursive(
		const graph_t & graph, const planar_embedding_t & planar_embedding,
		color_map_t & color_map, mark_map_t & mark_map,
		neighbor_range_map_t & neighbor_range_map, vertex_t p_0,
		std::size_t p_mark, std::size_t q_mark, std::size_t count,
		color_t p_color, color_t t_color, color_t q_color
	)
{
	if(neighbor_range_map[p_0].first == neighbor_range_map[p_0].second)
		return;
	
	vertex_t p = p_0, t_0 = p_0;
	auto p_edge_iter = neighbor_range_map[p].first;
	std::size_t below_p_mark = count++;
	
	while(p_edge_iter++ != neighbor_range_map[p].second)
	{
		if(p_edge_iter == planar_embedding[p].end())
			p_edge_iter = planar_embedding[p].begin();
		
		vertex_t n = get_incident_vertex(p, *p_edge_iter, graph);
		
		if(mark_map[n] == p_mark) {
			p_edge_iter = neighbor_range_map[n].first;
			p = n;
		}
		else if(mark_map[n] == q_mark) {
			if(t_0 != p_0) {
				if(p_edge_iter != neighbor_range_map[p].second) {
					auto temp = neighbor_range_map[p].first;
					neighbor_range_map[p].first = p_edge_iter;
			
					poh_color_recursive(
							graph, planar_embedding, color_map, mark_map,
							neighbor_range_map, p, p_mark, q_mark, count,
							p_color, t_color, q_color
						);
			
					neighbor_range_map[p].first = temp;
				}
				
				neighbor_range_map[p].second = p_edge_iter;
				remove_last_edge(p, neighbor_range_map, planar_embedding);
				
				break;
			}
			else {
				if(t_0 == p_0) t_0 = p;
				p_0 = p;
				neighbor_range_map[p].first = p_edge_iter;
			}
		}
		else if(mark_map[n] != below_p_mark) {
			if(t_0 == p_0) {
				t_0 = n;
				
				auto back_iter = find_neighbor_iterator(
						n, p, planar_embedding, graph
					);
				
				initialize(n, back_iter, neighbor_range_map, planar_embedding);
				remove_first_edge(n, neighbor_range_map, planar_embedding);
				
				neighbor_range_map[p].first = p_edge_iter;
			}
			
			mark_map[n] = below_p_mark;
		}
	}
	
	if(t_0 == p_0) return;
	
	vertex_t t = t_0;
	auto t_edge_iter = neighbor_range_map[t].first;
	std::size_t t_mark = count++;
	
	mark_map[t] = t_mark;
	color_map[t] = t_color;
	
	while(t_edge_iter++ != neighbor_range_map[t].second)
	{
		if(t_edge_iter == planar_embedding[t].end())
			t_edge_iter = planar_embedding[t].begin();
		
		vertex_t n = get_incident_vertex(t, *t_edge_iter, graph);
		
		if(mark_map[n] == below_p_mark) {
			auto n_back_iter = find_neighbor_iterator(
					n, t, planar_embedding, graph
				);
			
			initialize(n, n_back_iter, neighbor_range_map, planar_embedding);
			
			neighbor_range_map[t].second = t_edge_iter;
			
			t_edge_iter = neighbor_range_map[n].first;
			
			color_map[n] = t_color;
			mark_map[n] = t_mark;
			
			t = n;
		}
		else if(mark_map[n] == p_mark) {
			neighbor_range_map[t].second = t_edge_iter;
			remove_last_edge(t, neighbor_range_map, planar_embedding);
			
			break;
		}
	}
	
	poh_color_recursive(
			graph, planar_embedding, color_map, mark_map,
			neighbor_range_map, p_0, p_mark, t_mark, count,
			p_color, q_color, t_color
		);
	
	poh_color_recursive(
			graph, planar_embedding, color_map, mark_map,
			neighbor_range_map, t_0, t_mark, q_mark, count,
			t_color, p_color, q_color
		);
}

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
	typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename boost::property_traits<planar_embedding_t>::value_type
		::const_iterator edge_iterator_t;  
	typedef typename boost::iterator_property_map<
			std::vector<std::size_t>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> mark_map_t;
	typedef boost::iterator_property_map<
			typename std::vector<std::pair<edge_iterator_t, edge_iterator_t>>
				::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> neighbor_range_map_t;
	
	std::vector<std::size_t> mark_storage(boost::num_vertices(graph));
	mark_map_t mark_map(
			mark_storage.begin(), boost::get(boost::vertex_index, graph)
		);
	
	std::vector<std::pair<edge_iterator_t, edge_iterator_t>>
		neighbor_range_storage(boost::num_vertices(graph));
	neighbor_range_map_t neighbor_range_map(
				neighbor_range_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
	
	vertex_t l = *q_begin;
	
	for(vertex_iterator_t p_iter = p_begin; p_iter != p_end; ++p_iter) {
		vertex_t v = *p_iter;
		
		color_map[v] = c_0;
		mark_map[v] = 1;
		
		auto back_iter = find_neighbor_iterator(v, l, planar_embedding, graph);
		initialize(v, back_iter, neighbor_range_map, planar_embedding);
		
		l = v;
	}
	
	for(vertex_iterator_t q_iter = q_begin; q_iter != q_end; ++q_iter) {
		vertex_t v = *q_iter;
		
		color_map[v] = c_1;
		mark_map[v] = 2;
	}
	
	poh_color_recursive(
			graph, planar_embedding, color_map, mark_map,
			neighbor_range_map, *p_begin, 1, 2, 3,
			c_0, c_2, c_1
		);
}

#endif
