/*
 * hartman_skrekovski_color.hpp
 * Author: Aven Bross
 * 
 * Implementation of Hartman-Skrekovski path 3-list-coloring algorithm.
 */

#ifndef __PATH_COLORING_HPP
#define __PATH_COLORING_HPP

#include <vector>
#include <utility>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include "disjoint_set.hpp"
#include "augmented_embedding.hpp"

/*
 * list_color_properties
 * This structure stores all properties of a vertex relative to the Hartman-Skrekovski path 3-list-coloring
 * implementation below. Firstly, a vertex has a state that proceeds from "interior" to "on the outer face,"
 * and finally to "colored." Each vertex also tracks the current range of valid neighbors in its adjacency
 * list. Finally, each vertex recieves a mark describing which section of the outer face it belongs to
 * (note that this mark is only valid while the vertex has the state "on the outer face").
 */

template<typename index_graph>
class list_color_properties {
	public:
		typedef typename augmented_embedding<index_graph>::vertex_descriptor vertex_descriptor;
		typedef typename std::pair<std::size_t, std::size_t> neighbor_range;
		
		list_color_properties() : state(INTERIOR) {}
		
		int initialize(std::size_t start_index, std::size_t num_neighbors,
			int new_face_location, disjoint_set & face_locations)
		{
			// Vertices are initialized as they are added to the outer face
			state = ON_FACE;
			set_face_location(new_face_location, face_locations);
			
			// Set up current neighbor range
			degree = num_neighbors;
			range.first = start_index;
			range.second = (start_index + degree - 1) % degree;
			
			return face_location;
		}
		void color() {
			state = COLORED;
		}
		int set_face_location(int new_face_location, disjoint_set & face_locations) {
			if(!face_locations.exists(new_face_location)) {
				face_location = face_locations.make_next();
			}
			else {
				face_location = face_locations.find(new_face_location);
			}
			return face_location;
		}
		int get_face_location() {
			return face_location;
		}
		void set_range(const neighbor_range & new_range) {
			range = new_range;
		}
		const neighbor_range & get_range() const {
			return range;
		}
		void remove_begin() {
			range.first = (range.first + 1) % degree;
		}
		void remove_end() {
			range.second = (range.second + degree - 1) % degree;
		}
		std::pair<neighbor_range, neighbor_range> split_range(std::size_t mid_index) const {
			neighbor_range second_range(mid_index, range.second);
			neighbor_range first_range(range.first, mid_index);

			return std::pair<neighbor_range, neighbor_range>(first_range, second_range);
		}
		bool interior() const {
			return state == INTERIOR;
		}
		bool on_face() const {
			return state == ON_FACE;
		}
		bool colored() const {
			return state == COLORED;
		}

	private:
		enum vertex_state { INTERIOR, ON_FACE, COLORED };
	
		int face_location;	// Face location union find data
		vertex_state state;
		neighbor_range range;	// Current incidence range in embedding
		std::size_t degree;
};

/*
 * hartman_path_list_color
 * inputs: A weakly triangulated planar graph with vertex and edge indices (predefined boost properties),
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
		typename vertex_iterator
	>
void hartman_path_list_color(
		const index_graph & graph,
		const planar_embedding & embedding,
		color_list_map & color_list,
		color_map & coloring,
		vertex_iterator face_begin,
		vertex_iterator face_end
	)
{
	typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	
	// Construct our special embedded adjacency list with backwards lookup
	augmented_embedding<index_graph> plane_graph(graph, embedding);
	
	// Define the type for our property map
	typedef boost::iterator_property_map<
			typename std::vector<list_color_properties<index_graph> >::iterator,
			typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
		> list_color_property_map;
	
	// Make property map to store all properties
	std::vector<list_color_properties<index_graph> > list_color_property_storage(num_vertices(graph));
	list_color_property_map properties(list_color_property_storage.begin(),
		boost::get(boost::vertex_index, graph));
	
	// Initialize face locations
	disjoint_set face_locations;
	int before_y = -1;
	
	// Initialize all vertices on the given outer face
	for(auto vertex_iter = face_begin; vertex_iter != face_end; ++vertex_iter)
	{
		auto next = vertex_iter;
		if(++next == face_end) next = face_begin;
		
		vertex_descriptor last_vertex = *vertex_iter;
		vertex_descriptor vertex = *next;
		
		for(std::size_t i = 0; i < plane_graph[vertex].size(); ++i)
		{
			vertex_descriptor neighbor = plane_graph[vertex][i].neighbor;
			
			if(neighbor == last_vertex) {
				before_y = properties[vertex].initialize(
						i, plane_graph[vertex].size(),
						before_y, face_locations
					);
			}
		}
	}
	
	// Star with x and y as the beginning and end of the outer face list given
	vertex_descriptor x = *face_begin, y = *(--face_end);
	
	// Recursively list color the graph
	hartman_path_list_color_recursive(
			plane_graph,
			properties, face_locations,
			color_list, coloring,
			x, y, x,
			-1, before_y, -1
		);
}

/*
 * hartman_path_list_color_recursive
 * inputs: An augmented embedding representing a weakly triangulated plane graph,
 *         a read-write vertex property map assigning a list_color_properties object to each vertex
 *         (the properties for the outer face of the current subgraph must already be assigned),
 *         a disjoint set structure storing integer identifiers for segments of the outer face,
 *         a read-write vertex property map assigning a range of colors to each vertex (each vertex must
 *         recieve a range of at least 3 colors if interior, at least two colors if on the outer face,
 *         and at least one color if they are x or y),
 *         a read-write vertex property map to which the coloring will be assigned,
 *         vertices x, y, and p on the outer face with p between x and y clockwise,
 *         and finally the current identifiers for the outer face segments of vertices clockwise between
 *         x and p (before_p), p and y (before_y), and y and x (before_x) respectively.
 *
 * outputs: Colors from the provided lists will be assigned to the vertices of the subgraph bounded by
 *          the given outer face such that each color class induces a disjoint union of paths. Furthermore,
 *          if x = y, then no neighbors of x (= y) will recieve the same color as x (= y). If x != y,
 *          then x and y will have at most one neighbor sharing a color in the resulting coloring.
 */

template<
		typename index_graph, typename list_color_property_map,
		typename color_list_map, typename color_map,
		typename vertex_descriptor = typename augmented_embedding<index_graph>::vertex_descriptor
	>
void hartman_path_list_color_recursive(
		const augmented_embedding<index_graph> & plane_graph,
		list_color_property_map & properties, disjoint_set & face_locations,
		color_list_map & color_list, color_map & coloring,
		vertex_descriptor x, vertex_descriptor y, vertex_descriptor p,
		int before_p, int before_y, int before_x
	)
{
	auto range(properties[p].get_range());
	
	if(!properties[p].colored()) {
		p = x;
		properties[p].color();
		coloring[p] = color_list[p].front();
	}
	
	vertex_descriptor new_x = x, new_y = y;
	std::size_t current_index = range.first;
	
	do {
		current_index %= plane_graph[p].size();
		
		vertex_descriptor n = plane_graph[p][current_index].neighbor;
		std::size_t back_index = plane_graph[p][current_index].back_index;
		int n_location = properties[n].get_face_location();
		
		if(properties[n].interior()) {
			before_p = properties[n].initialize(
					back_index, plane_graph[n].size(),
					before_p, face_locations
				);
			
			color_list[n].erase(
					std::remove(color_list[n].begin(), color_list[n].end(),
					coloring[p]), color_list[n].end()
				);
			
			properties[n].remove_begin();
		}
		else if(current_index == range.first) {
			if(current_index == range.second) {
				if(n != y) {
					color_list[n].erase(
							std::remove(color_list[n].begin(), color_list[n].end(),
							coloring[p]), color_list[n].end()
						);
				}
				if(!properties[n].colored()) {
					properties[n].color();
					coloring[n] = color_list[n].front();
				}
			
				break;
			}
			else if(n == y) {
				coloring[y] = coloring[p];
				properties[y].color();
				
				if(x == p) new_x = y;
				
				hartman_path_list_color_recursive(
						plane_graph,
						properties, face_locations,
						color_list, coloring,
						new_x, p, y,
						-1, -1, before_y
					);
				
				break;
			}
			else {
				if(x == p) {
					new_x = n;
					before_p = properties[n].set_face_location(before_p, face_locations);
				}
		
				color_list[n].erase(
						std::remove(color_list[n].begin(), color_list[n].end(),
						coloring[p]), color_list[n].end()
					);
		
				properties[n].remove_end();
			}
		}
		else {
			auto p_ranges = properties[p].split_range(current_index);
			auto n_ranges = properties[n].split_range(back_index);
	
			if(p == y) {
				color_list[n].erase(
						std::remove(color_list[n].begin(), color_list[n].end(),
						coloring[p]), color_list[n].end()
					);
		
				if(p == x) {
					properties[n].set_range(n_ranges.second);
					properties[n].remove_begin();
			
					hartman_path_list_color_recursive(
							plane_graph,
							properties, face_locations,
							color_list, coloring,
							new_x, n, new_x,
							-1, before_p, before_y
						);
			
					if(back_index != n_ranges.first.first) {
						properties[n].set_range(n_ranges.first);
						properties[p].set_range(p_ranges.second);
			
						hartman_path_list_color_recursive(
								plane_graph,
								properties, face_locations,
								color_list, coloring,
								p, p, p,
								-1, before_y, -1
							);
					}
					
					break;
				}
				else if(current_index == range.second) {
					before_p = properties[n].set_face_location(before_p, face_locations);
					properties[n].remove_begin();
					
					hartman_path_list_color_recursive(
							plane_graph,
							properties, face_locations,
							color_list, coloring,
							new_x, n, new_x,
							-1, before_p, before_x
						);
					
					break;
				}
				else {
					new_y = n;
				}
			}
			
			if(face_locations.compare(n_location, before_y) || n == y) {
				vertex_descriptor new_p = n;
				
				if(p != y && std::find(color_list[n].begin(), color_list[n].begin(), coloring[p])
					!= color_list[n].end() && !properties[n].colored())
				{
					coloring[n] = coloring[p];
					properties[n].color();
				}
				else if(!properties[n].colored() || coloring[n] != coloring[p]) {
					new_p = new_x;
					before_y = face_locations.take_union(before_p, before_y);
					before_p = -1;
					
					color_list[n].erase(
							std::remove(color_list[n].begin(), color_list[n].end(),
							coloring[p]), color_list[n].end()
						);
				}
	
				properties[n].set_range(n_ranges.second);
				properties[n].remove_begin();
	
				hartman_path_list_color_recursive(
						plane_graph,
						properties, face_locations,
						color_list, coloring,
						new_x, new_y, new_p,
						before_p, before_y, before_x
					);
	
				if(back_index != n_ranges.first.first) {
					properties[n].set_range(n_ranges.first);
					properties[p].set_range(p_ranges.second);
			
					hartman_path_list_color_recursive(
							plane_graph,
							properties, face_locations,
							color_list, coloring,
							p, n, p,
							-1, before_y, -1
						);
				}
			}
			else if(face_locations.compare(n_location, before_x)) {
				color_list[n].erase(
						std::remove(color_list[n].begin(), color_list[n].end(),
						coloring[p]), color_list[n].end()
					);
		
				properties[n].set_range(n_ranges.second);
				properties[n].remove_begin();
		
				hartman_path_list_color_recursive(
						plane_graph,
						properties, face_locations,
						color_list, coloring,
						new_x, n, new_x,
						-1, before_p, before_x
					);
		
				properties[n].set_range(n_ranges.first);
				properties[p].set_range(p_ranges.second);
		
				hartman_path_list_color_recursive(
						plane_graph,
						properties, face_locations,
						color_list, coloring,
						n, y, p,
						-1, before_y, before_x
					);
			}
			else {
				color_list[n].erase(
						std::remove(color_list[n].begin(), color_list[n].end(),
						coloring[p]), color_list[n].end()
					);
		
				properties[n].set_range(n_ranges.first);
				properties[p].set_range(p_ranges.second);
		
				hartman_path_list_color_recursive(
						plane_graph,
						properties, face_locations,
						color_list, coloring,
						x, y, p,
						before_p, before_y, before_x
					);
				
				properties[n].set_range(n_ranges.second);
				properties[n].remove_begin();
				
				hartman_path_list_color_recursive(
						plane_graph,
						properties, face_locations,
						color_list, coloring,
						n, n, n,
						-1, before_p, -1
					);
			}
	
			break;
		}
	} while(current_index++ != range.second);
}

#endif