/*
 * hartman_skrekovski_choose.hpp
 * Author: Aven Bross
 * 
 * Implementation of Hartman-Skrekovski path 3-choosing algorithm.
 */

#ifndef __HARTMAN_SKREKOVSKI_CHOOSE_HPP
#define __HARTMAN_SKREKOVSKI_CHOOSE_HPP

// STL headers
#include <vector>
#include <utility>
#include <algorithm>

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Local project headers
#include "disjoint_set.hpp"
#include "incidence_list_helpers.hpp"


const static int INTERIOR = 0, ON_FACE = 1, COLORED = 2;


template<
		typename vertex_t, typename face_location_map_t
	>
static inline int set_face_location(
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

template<
		typename graph_t, typename planar_embedding_t,
		typename face_location_map_t, typename state_map_t,
		typename neighbor_range_map_t, typename color_list_map_t,
		typename color_map_t, typename vertex_t 
			= typename boost::graph_traits<graph_t>::vertex_descriptor
	>
static void hartman_skrekovski_choose_recursive(
		const graph_t & graph, const planar_embedding_t & planar_embedding,
		face_location_map_t & face_location_map,
		disjoint_set_t & face_location_sets, state_map_t & state_map,
		neighbor_range_map_t & neighbor_range_map,
		color_list_map_t & color_list_map, color_map_t & color_map,
		vertex_t x, vertex_t y, vertex_t p,
		int before_p, int before_y, int before_x
	)
{
	typedef typename boost::property_traits<planar_embedding_t>::value_type
		::const_iterator edge_iterator_t;
	typedef typename std::pair<edge_iterator_t, edge_iterator_t>
		neighbor_range_t; 
	
	auto neighbor_range = neighbor_range_map[p];
	
	if(state_map[p] != COLORED) {
		state_map[p] = COLORED;
		color_map[p] = color_list_map[p].front();
	}
	
	vertex_t new_x = x, new_y = y;
	auto edge_iter = neighbor_range.first;
	
	do {
		if(edge_iter == planar_embedding[p].end())
			edge_iter = planar_embedding[p].begin();
		
		vertex_t n = get_incident_vertex(p, *edge_iter, graph);
		int n_location = face_location_map[n];
		
		if(state_map[n] == INTERIOR) {
			auto back_iter = find_neighbor_iterator(
					n, p, planar_embedding, graph
				);
			
			state_map[n] = ON_FACE;
			
			before_p = set_face_location(
					n, before_p, face_location_map, face_location_sets
				);
			
			initialize(n, back_iter, neighbor_range_map, planar_embedding);
			
			color_list_map[n].remove(color_map[p]);
			
			remove_first_edge(n, neighbor_range_map, planar_embedding);
		}
		else if(edge_iter == neighbor_range.first) {
			if(edge_iter == neighbor_range.second) {
				if(n != x && n != y) {
					color_list_map[n].remove(color_map[p]);
				}
				
				if(state_map[n] != COLORED) {
					state_map[n] = COLORED;
					color_map[n] = color_list_map[n].front();
				}
			
				break;
			}
			else if(n == y) {
				if(state_map[y] != COLORED) {
					state_map[y] = COLORED;
					color_map[y] = color_map[p];
				}
				
				if(x == p) new_x = y;
				
				hartman_skrekovski_choose_recursive(
						graph, planar_embedding,
						face_location_map, face_location_sets,
						state_map, neighbor_range_map, 
						color_list_map, color_map,
						new_x, p, y,
						-1, -1, before_y
					);
				
				break;
			}
			else {
				if(x == p) {
					new_x = n;
					
					before_p = set_face_location(
							n, before_p, face_location_map, face_location_sets
						);
				}
				
				color_list_map[n].remove(color_map[p]);
		
				remove_last_edge(n, neighbor_range_map, planar_embedding);
			}
		}
		else {
			auto back_iter = find_neighbor_iterator_restricted(
					n, p, neighbor_range_map[n].first,
					neighbor_range_map[n].second, planar_embedding, graph
				);
			
			auto n_ranges = split_range(
					n, back_iter, neighbor_range_map, planar_embedding
				);
			
			neighbor_range_map[p] = neighbor_range_t(
					edge_iter, neighbor_range.second
				);
	
			if(p == y) {
				color_list_map[n].remove(color_map[p]);
		
				if(p == x) {
					neighbor_range_map[n] = n_ranges.second;
					
					hartman_skrekovski_choose_recursive(
							graph, planar_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							new_x, n, new_x,
							-1, before_p, before_y
						);
			
					if(back_iter != n_ranges.first.first) {
						neighbor_range_map[n] = n_ranges.first;
						
						hartman_skrekovski_choose_recursive(
								graph, planar_embedding,
								face_location_map, face_location_sets,
								state_map, neighbor_range_map, 
								color_list_map, color_map,
								n, p, p,
								-1, -1, before_y
							);
					}
					
					break;
				}
				else if(edge_iter == neighbor_range.second) {
					neighbor_range_map[n] = n_ranges.second;
					
					hartman_skrekovski_choose_recursive(
							graph, planar_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							new_x, n, new_x,
							-1, before_p, before_x
						);
					
					break;
				}
				
				new_y = n;
			}
			
			if(face_location_sets.compare(n_location, before_y) || n == y) {
				vertex_t new_p = n;
				
				if(std::find(
								color_list_map[n].begin(),
								color_list_map[n].begin(), color_map[p]
							) != color_list_map[n].end() 
						&& state_map[n] != COLORED
					)
				{
					color_map[n] = color_map[p];
					state_map[n] = COLORED;
				}
				else if(state_map[n] != COLORED || color_map[n] != color_map[p]){
					new_p = new_x;
					
					before_y = face_location_sets.take_union(
							before_p, before_y
						);
					
					before_p = -1;
				}
	
				neighbor_range_map[n] = n_ranges.second;
				
				hartman_skrekovski_choose_recursive(
						graph, planar_embedding,
						face_location_map, face_location_sets,
						state_map, neighbor_range_map, 
						color_list_map, color_map,
						new_x, new_y, new_p,
						before_p, before_y, before_x
					);
	
				if(back_iter != n_ranges.first.first) {
					neighbor_range_map[n] = n_ranges.first;
				
					hartman_skrekovski_choose_recursive(
							graph, planar_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							p, n, p,
							-1, before_y, -1
						);
				}
			}
			else if(face_location_sets.compare(n_location, before_x)) {
				color_list_map[n].remove(color_map[p]);
		
				neighbor_range_map[n] = n_ranges.second;
		
				hartman_skrekovski_choose_recursive(
						graph, planar_embedding,
						face_location_map, face_location_sets,
						state_map, neighbor_range_map, 
						color_list_map, color_map,
						new_x, n, new_x,
						-1, before_p, before_x
					);
		
				neighbor_range_map[n] = n_ranges.first;
		
				hartman_skrekovski_choose_recursive(
						graph, planar_embedding,
						face_location_map, face_location_sets,
						state_map, neighbor_range_map, 
						color_list_map, color_map,
						n, y, p,
						-1, before_y, before_x
					);
			}
			else {
				color_list_map[n].remove(color_map[p]);
		
				neighbor_range_map[n] = n_ranges.first;
		
				hartman_skrekovski_choose_recursive(
						graph, planar_embedding,
						face_location_map, face_location_sets,
						state_map, neighbor_range_map, 
						color_list_map, color_map,
						x, y, p,
						before_p, before_y, before_x
					);
				
				neighbor_range_map[n] = n_ranges.second;
				
				hartman_skrekovski_choose_recursive(
						graph, planar_embedding,
						face_location_map, face_location_sets,
						state_map, neighbor_range_map, 
						color_list_map, color_map,
						n, n, n,
						-1, before_p, -1
					);
			}
	
			break;
		}
	} while(edge_iter++ != neighbor_range.second);
}


/*
 * hartman_skrekovski_choose
 * inputs: A weakly triangulated planar graph with vertex indices (predefined
 * boost property), a valid planar embedding of the graph (modeling the boost
 * PlanarEmbedding concept), a read-able vertex property map assigning a range
 * of colors to each vertex (each vertex must recieve a range of at least 3
 * colors if interior, and at least two colors if on the outer face), a
 * read-write-able vertex property map to which the coloring will be assigned,
 * and finally a pair of iterators providing the outer face of the graph in
 * clockwise order.
 * 
 * outputs: The coloring will be a valid assignment of colors from the input
 * color lists such that each color class induces a disjoint union of paths.
 */
 
template<
		typename graph_t, typename planar_embedding_t,
		typename color_list_map_t, typename color_map_t,
		typename face_iterator_t
	>
void hartman_skrekovski_choose(
		const graph_t & graph,
		const planar_embedding_t & planar_embedding,
		const color_list_map_t & color_list_map,
		color_map_t & color_map,
		face_iterator_t face_begin,
		face_iterator_t face_end
	)
{
	typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename color_map_t::value_type color_t;
	typedef typename boost::graph_traits<graph_t>::vertex_iterator
		vertex_iterator_t;
	typedef typename boost::property_traits<planar_embedding_t>::value_type
		::const_iterator edge_iterator_t;
	
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
	
	std::vector<int> face_location_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			std::vector<int>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> face_location_map(
				face_location_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
		
	std::vector<int> state_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			std::vector<int>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> state_map(
				state_storage.begin(), boost::get(boost::vertex_index, graph)
			);
	
	std::vector<std::list<color_t>> color_list_storage(
			boost::num_vertices(graph)
		);
	boost::iterator_property_map<
			typename std::vector<std::list<color_t>>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> color_list_map_copy(
				color_list_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
		
	vertex_iterator_t vertex_iter, vertex_end;
	for(boost::tie(vertex_iter, vertex_end) = vertices(graph);
			vertex_iter != vertex_end; vertex_iter++
		)
	{
		vertex_t v = *vertex_iter;
		
		std::copy(
			color_list_map[v].begin(), color_list_map[v].end(),
			std::back_inserter(color_list_map_copy[v])
		);
	}
	
	disjoint_set_t face_location_sets;
	int before_y = -1;
	
	for(auto face_iter = face_begin; face_iter != face_end; ++face_iter)
	{
		auto next = face_iter;
		if(++next == face_end) next = face_begin;
		
		vertex_t l = *face_iter;
		vertex_t v = *next;
		
		state_map[v] = ON_FACE;
		
		before_y = set_face_location(
				v, before_y, face_location_map, face_location_sets
			);
		
		auto back_iter = find_neighbor_iterator(v, l, planar_embedding, graph);
		initialize(v, back_iter, neighbor_range_map, planar_embedding);
	}
	
	vertex_t x = *face_begin, y = *(--face_end);
	
	hartman_skrekovski_choose_recursive(
			graph, planar_embedding,
			face_location_map, face_location_sets,
			state_map, neighbor_range_map, 
			color_list_map_copy, color_map,
			x, y, x,
			-1, before_y, -1
		);
}

#endif
