/*
 * hartman_skrekovski_choose.hpp
 * Author: Aven Bross
 * 
 * Implementation of Hartman-Skrekovski path 3-choose algorithm.
 */

#ifndef __PATH_COLORING_HPP
#define __PATH_COLORING_HPP

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

template<typename vertex_descriptor, typename neighbor_range_map, typename planar_embedding>
static inline void remove_first_edge(vertex_descriptor v, neighbor_range_map & neighbor_ranges,
	const planar_embedding & embedding)
{
	if(++neighbor_ranges[v].first == embedding[v].end()) {
		neighbor_ranges[v].first = embedding[v].begin();
	}
}

template<typename vertex_descriptor, typename neighbor_range_map, typename planar_embedding>
static inline void remove_last_edge(vertex_descriptor v, neighbor_range_map & neighbor_ranges,
	const planar_embedding & embedding)
{
	if(neighbor_ranges[v].second == embedding[v].begin()) {
		neighbor_ranges[v].second = embedding[v].end();
	}
	--neighbor_ranges[v].second;
}

template<
		typename vertex_descriptor, typename edge_iterator, typename neighbor_range_map,
		typename neighbor_range = typename std::pair<edge_iterator, edge_iterator>,
		typename neighbor_range_pair = typename std::pair<neighbor_range, neighbor_range>
	>
static inline neighbor_range_pair split_range(vertex_descriptor v, edge_iterator mid_iter,
	neighbor_range_map & neighbor_ranges)
{
	return neighbor_range_pair(
			neighbor_range(neighbor_ranges[v].first, mid_iter),
			neighbor_range(mid_iter, neighbor_ranges[v].second)
		);
}

template<typename vertex_descriptor, typename face_location_map>
static inline int set_face_location(vertex_descriptor v, int new_face_location,
	face_location_map & face_locations, disjoint_set & face_location_sets)
{
	if(!face_location_sets.exists(new_face_location)) {
		face_locations[v] = face_location_sets.make_next();
	}
	else {
		face_locations[v] = face_location_sets.find(new_face_location);
	}
	return face_locations[v];
}

template<
		typename vertex_descriptor, typename face_location_map,
		typename state_map, typename neighbor_range_map, typename planar_embedding,
		typename edge_iterator = typename boost::property_traits<planar_embedding>::value_type::const_iterator
	>
static inline int initialize(vertex_descriptor v, edge_iterator start_iter, int new_face_location,
	face_location_map & face_locations, disjoint_set & face_location_sets, state_map & states,
	neighbor_range_map & neighbor_ranges, const planar_embedding & embedding)
{
	states[v] = ON_FACE;
	
	neighbor_ranges[v].first = start_iter;
	neighbor_ranges[v].second = start_iter;
	remove_last_edge(v, neighbor_ranges, embedding);
	
	return set_face_location(v, new_face_location, face_locations, face_location_sets);
}

template<
		typename index_graph, typename planar_embedding,
		typename face_location_map, typename state_map,
		typename neighbor_range_map,
		typename color_list_map, typename color_map,
		typename vertex_descriptor = typename boost::graph_traits<index_graph>::vertex_descriptor
	>
static void hartman_skrekovski_choose_recursive(
		const index_graph & graph, const planar_embedding & embedding,
		face_location_map & face_locations, disjoint_set & face_location_sets,
		state_map & states, neighbor_range_map & neighbor_ranges,
		color_list_map & color_lists, color_map & coloring,
		vertex_descriptor x, vertex_descriptor y, vertex_descriptor p,
		int before_p, int before_y, int before_x
	)
{
	auto range = neighbor_ranges[p];
	
	if(states[p] != COLORED) {
		states[p] = COLORED;
		coloring[p] = color_lists[p].front();
	}
	
	vertex_descriptor new_x = x, new_y = y;
	auto edge_iter = range.first;
	
	do {
		if(edge_iter == embedding[p].end()) edge_iter = embedding[p].begin();
		
		vertex_descriptor n = get_incident_vertex(p, *edge_iter, graph);
		int n_location = face_locations[n];
		
		if(states[n] == INTERIOR) {
			auto back_iter = find_neighbor_iterator(n, p, embedding, graph);
			before_p = initialize(n, back_iter, before_p, face_locations,
				face_location_sets, states, neighbor_ranges, embedding);
			
			color_lists[n].remove(coloring[p]);
			
			remove_first_edge(n, neighbor_ranges, embedding);
		}
		else if(edge_iter == range.first) {
			if(edge_iter == range.second) {
				if(n != x && n != y) {
					color_lists[n].remove(coloring[p]);
				}
				
				if(states[n] != COLORED) {
					states[n] = COLORED;
					coloring[n] = color_lists[n].front();
				}
			
				break;
			}
			else if(n == y) {
				if(states[y] != COLORED) {
					states[y] = COLORED;
					coloring[y] = coloring[p];
				}
				
				if(x == p) new_x = y;
				
				hartman_skrekovski_choose_recursive(
						graph, embedding,
						face_locations, face_location_sets,
						states, neighbor_ranges, 
						color_lists, coloring,
						new_x, p, y,
						-1, -1, before_y
					);
				
				break;
			}
			else {
				if(x == p) {
					new_x = n;
					before_p = set_face_location(n, before_p, face_locations, face_location_sets);
				}
				
				color_lists[n].remove(coloring[p]);
		
				remove_last_edge(n, neighbor_ranges, embedding);
			}
		}
		else {
			auto back_iter = find_neighbor_iterator_restricted(n, p, neighbor_ranges[n].first,
				neighbor_ranges[n].second, embedding, graph);
			
			auto p_ranges = split_range(p, edge_iter, neighbor_ranges);
			auto n_ranges = split_range(n, back_iter, neighbor_ranges);
	
			if(p == y) {
				color_lists[n].remove(coloring[p]);
		
				if(p == x) {
					neighbor_ranges[n] = n_ranges.second;
					remove_first_edge(n, neighbor_ranges, embedding);
					
					hartman_skrekovski_choose_recursive(
							graph, embedding,
							face_locations, face_location_sets,
							states, neighbor_ranges, 
							color_lists, coloring,
							new_x, n, new_x,
							-1, before_p, before_y
						);
			
					if(back_iter != n_ranges.first.first) {
						neighbor_ranges[n] = n_ranges.first;
						neighbor_ranges[p] = p_ranges.second;
						
						hartman_skrekovski_choose_recursive(
								graph, embedding,
								face_locations, face_location_sets,
								states, neighbor_ranges, 
								color_lists, coloring,
								n, p, p,
								-1, -1, before_y
							);
					}
					
					break;
				}
				else if(edge_iter == range.second) {
					remove_first_edge(n, neighbor_ranges, embedding);
					
					hartman_skrekovski_choose_recursive(
							graph, embedding,
							face_locations, face_location_sets,
							states, neighbor_ranges, 
							color_lists, coloring,
							new_x, n, new_x,
							-1, before_p, before_x
						);
					
					break;
				}
				
				new_y = n;
			}
			
			if(face_location_sets.compare(n_location, before_y) || n == y) {
				vertex_descriptor new_p = n;
				
				if(std::find(color_lists[n].begin(), color_lists[n].begin(), coloring[p])
					!= color_lists[n].end() && states[n] != COLORED)
				{
					coloring[n] = coloring[p];
					states[n] = COLORED;
				}
				else if(states[n] != COLORED || coloring[n] != coloring[p]) {
					new_p = new_x;
					before_y = face_location_sets.take_union(before_p, before_y);
					before_p = -1;
				}
	
				neighbor_ranges[n] = n_ranges.second;
				remove_first_edge(n, neighbor_ranges, embedding);
				
				hartman_skrekovski_choose_recursive(
						graph, embedding,
						face_locations, face_location_sets,
						states, neighbor_ranges, 
						color_lists, coloring,
						new_x, new_y, new_p,
						before_p, before_y, before_x
					);
	
				if(back_iter != n_ranges.first.first) {
					neighbor_ranges[n] = n_ranges.first;
					neighbor_ranges[p] = p_ranges.second;
				
					hartman_skrekovski_choose_recursive(
							graph, embedding,
							face_locations, face_location_sets,
							states, neighbor_ranges, 
							color_lists, coloring,
							p, n, p,
							-1, before_y, -1
						);
				}
			}
			else if(face_location_sets.compare(n_location, before_x)) {
				color_lists[n].remove(coloring[p]);
		
				neighbor_ranges[n] = n_ranges.second;
				remove_first_edge(n, neighbor_ranges, embedding);
		
				hartman_skrekovski_choose_recursive(
						graph, embedding,
						face_locations, face_location_sets,
						states, neighbor_ranges, 
						color_lists, coloring,
						new_x, n, new_x,
						-1, before_p, before_x
					);
		
				neighbor_ranges[n] = n_ranges.first;
				neighbor_ranges[p] = p_ranges.second;
		
				hartman_skrekovski_choose_recursive(
						graph, embedding,
						face_locations, face_location_sets,
						states, neighbor_ranges, 
						color_lists, coloring,
						n, y, p,
						-1, before_y, before_x
					);
			}
			else {
				color_lists[n].remove(coloring[p]);
		
				neighbor_ranges[n] = n_ranges.first;
				neighbor_ranges[p] = p_ranges.second;
		
				hartman_skrekovski_choose_recursive(
						graph, embedding,
						face_locations, face_location_sets,
						states, neighbor_ranges, 
						color_lists, coloring,
						x, y, p,
						before_p, before_y, before_x
					);
				
				neighbor_ranges[n] = n_ranges.second;
				remove_first_edge(n, neighbor_ranges, embedding);
				
				hartman_skrekovski_choose_recursive(
						graph, embedding,
						face_locations, face_location_sets,
						states, neighbor_ranges, 
						color_lists, coloring,
						n, n, n,
						-1, before_p, -1
					);
			}
	
			break;
		}
	} while(edge_iter++ != range.second);
}


/*
 * hartman_skrekovski_choose
 * inputs: A weakly triangulated planar graph with vertex indices (predefined boost property),
 *         a valid planar embedding of the graph (modeling the boost PlanarEmbedding concept),
 *         a read-write vertex property map assigning a range of colors to each vertex (each vertex must
 *         recieve a range of at least 3 colors if interior, and at least two colors if on the outer face),
 *         a read-write vertex property map to which the coloring will be assigned,
 *         and finally a pair of iterators providing the outer face of the graph in clockwise order.
 * 
 * outputs: The coloring will be a valid assignment of colors from the input color lists such that
 *          each color class induces a disjoint union of paths.
 */
 
template<
		typename index_graph, typename planar_embedding,
		typename color_list_map, typename color_map,
		typename face_iterator
	>
void hartman_skrekovski_choose(
		const index_graph & graph,
		const planar_embedding & embedding,
		const color_list_map & color_lists,
		color_map & coloring,
		face_iterator face_begin,
		face_iterator face_end
	)
{
	typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	typedef typename color_map::value_type color_type;
	typedef typename boost::graph_traits<index_graph>::vertex_iterator vertex_iterator;
	typedef typename boost::property_traits<planar_embedding>::value_type::const_iterator edge_iterator;
	
	std::vector<std::pair<edge_iterator, edge_iterator>> neighbor_range_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			typename std::vector<std::pair<edge_iterator, edge_iterator>>::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> neighbor_ranges(neighbor_range_storage.begin(), boost::get(boost::vertex_index, graph));
	
	std::vector<int> face_location_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			std::vector<int>::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> face_locations(face_location_storage.begin(), boost::get(boost::vertex_index, graph));
		
	std::vector<int> state_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			std::vector<int>::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> states(state_storage.begin(), boost::get(boost::vertex_index, graph));
	
	std::vector<std::list<color_type>> color_list_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			typename std::vector<std::list<color_type>>::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> color_lists_copy(color_list_storage.begin(), boost::get(boost::vertex_index, graph));
		
	vertex_iterator vertex_iter, vertex_end;
	for(boost::tie(vertex_iter, vertex_end) = vertices(graph); vertex_iter != vertex_end; vertex_iter++) {
		vertex_descriptor v = *vertex_iter;
		std::copy(color_lists[v].begin(), color_lists[v].end(), std::back_inserter(color_lists_copy[v]));
	}
	
	disjoint_set face_location_sets;
	int before_y = -1;
	
	for(auto face_iter = face_begin; face_iter != face_end; ++face_iter)
	{
		auto next = face_iter;
		if(++next == face_end) next = face_begin;
		
		vertex_descriptor l = *face_iter;
		vertex_descriptor v = *next;
		
		auto back_iter = find_neighbor_iterator(v, l, embedding, graph);
		before_y = initialize(v, back_iter, before_y, face_locations,
			face_location_sets, states, neighbor_ranges, embedding);
	}
	
	vertex_descriptor x = *face_begin, y = *(--face_end);
	
	hartman_skrekovski_choose_recursive(
			graph, embedding,
			face_locations, face_location_sets,
			states, neighbor_ranges, 
			color_lists_copy, coloring,
			x, y, x,
			-1, before_y, -1
		);
}

#endif
