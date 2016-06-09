/*
 * hartman_skrekovski_color.hpp
 * Author: Aven Bross
 * 
 * Implementation of Hartman-Skrekovski path 3-list-coloring algorithm.
 */

#ifndef __PATH_COLORING_HPP
#define __PATH_COLORING_HPP

#include <vector>
#include <map>
#include <stdexcept>
#include <utility>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include "disjoint_set.hpp"
#include "incidence_list_helpers.hpp"

namespace boost {
	 
	// Stores all vertex propreties for path 3-list-coloring algorithm
	template<typename vertex_descriptor, typename edge_iterator>
	class list_color_properties {
		public:
			typedef typename std::pair<edge_iterator, edge_iterator> edge_iterator_pair;

			list_color_properties() : state(INTERIOR) {}

			template<typename index_graph, typename planar_embedding>
			int initialize(vertex_descriptor vertex, vertex_descriptor adding_vertex,
				int new_face_location, const index_graph & graph, const planar_embedding & embedding,
				disjoint_set & face_location_sets)
			{
				// Vertices are initialized as they are added to the outer face
				state = ON_FACE;
				set_face_location(new_face_location, face_location_sets);
	
				// Set up incidence list
				range.first = embedding[vertex].begin();
				range.second = embedding[vertex].end();
	
				// Loop through all neighbors
				for(edge_iterator edge_iter = range.first; edge_iter != range.second; ++edge_iter) {
					vertex_descriptor neighbor = get_incident_vertex(vertex, *edge_iter, graph);
		
					if(neighbor == adding_vertex) {
						// Orient our incidence list to start at the adding vertex
						current_range.first = edge_iter;
						current_range.second = edge_iter;
			
						// Correction to make sure list wraps correctly
						if(current_range.second == range.first) {
							current_range.second = range.second;
						}
					}
		
					// Remember incidence list location of each neighbor
					neighbor_iterator[neighbor] = edge_iter;
				}
				
				return face_location;
			}
			void color() {
				state = COLORED;
			}
			int set_face_location(int new_face_location, disjoint_set & face_location_sets) {
				if(!face_location_sets.exists(new_face_location)) {
					face_location = face_location_sets.make_next();
				}
				else {
					face_location = face_location_sets.find(new_face_location);
				}
				return face_location;
			}
			int get_face_location() {
				return face_location;
			}
			void set_current_range(const edge_iterator_pair & new_range) {
				current_range = new_range;
			}
			edge_iterator_pair get_current_range() const {
				return current_range;
			}
			edge_iterator_pair get_range() const {
				return range;
			}
			void remove_begin() {
				++current_range.first;
				if(current_range.first == range.second && current_range.first != current_range.second) {
					current_range.first = range.first;
				}
			}
			void remove_end() {
				--current_range.second;
				if(current_range.second == range.first && current_range.first != current_range.second) {
					current_range.second = range.second;
				}
			}
			std::pair<edge_iterator_pair, edge_iterator_pair> split_range_at(vertex_descriptor neighbor) const {
				edge_iterator split_iter = neighbor_iterator.at(neighbor);
				edge_iterator_pair second_range(split_iter, current_range.second);
				edge_iterator_pair first_range(current_range.first, ++split_iter);
	
				return std::pair<edge_iterator_pair, edge_iterator_pair>(first_range, second_range);
			}
			bool single_neighbor() const {
				edge_iterator temp_first = current_range.first;
				return ++temp_first == current_range.second;
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
			edge_iterator_pair range;	// Incidence range from embedding
			edge_iterator_pair current_range;	// Current incidence range in embedding
			std::map<vertex_descriptor, edge_iterator> neighbor_iterator;	// Reverse incidence list lookup
	};
	 
	 template<typename index_graph, typename planar_embedding, typename color_list_map, typename color_map,
	 	typename vertex_iterator>
	void hartman_path_list_color(const index_graph & graph, const planar_embedding & embedding,
		color_list_map & color_list, color_map & coloring, vertex_iterator face_begin,
		vertex_iterator face_end)
	{
		typedef typename property_traits<planar_embedding>::value_type::const_iterator edge_iterator;
		typedef typename graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
		
		// Define the type for our property map
		typedef iterator_property_map<
				typename std::vector<list_color_properties<vertex_descriptor, edge_iterator> >::iterator,
				typename property_map<index_graph, vertex_index_t>::type
			> list_color_property_map;
		
		// Make property map to store all properties
		std::vector<list_color_properties<vertex_descriptor, edge_iterator> >
			property_storage(num_vertices(graph));
		list_color_property_map properties(property_storage.begin(), get(vertex_index, graph));
		
		// Initialize face locations
		disjoint_set face_location_sets;
		int before_y = -1;
		
		// Initialize all vertices on the given outer face
		for(auto vertex_iter = face_begin; vertex_iter != face_end; ++vertex_iter)
		{
			// Grab next vertex on outer face
			auto next = vertex_iter;
			if(++next == face_end) next = face_begin;
			
			// Initialize properties of the next vertex on the face using the current face edge
			before_y = properties[*next].initialize(*next, *vertex_iter, before_y, graph, embedding,
				face_location_sets);
		}
		
		// Star with x and y as the beginning and end of the outer face list given
		vertex_descriptor x = *face_begin, y = *(--face_end);
		
		// Recursively list color the graph
		hartman_path_list_color_recursive(graph, embedding, color_list, properties, face_location_sets,
			coloring, x, properties[x].get_current_range(), y, properties[y].get_current_range(), x,
			-1, before_y, -1);
	}
	
	template<
		typename index_graph, typename planar_embedding, typename color_list_map,
		typename color_map, typename list_color_property_map,
		typename edge_iterator = typename property_traits<planar_embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<index_graph>::vertex_descriptor,
		typename edge_iterator_pair = typename std::pair<edge_iterator, edge_iterator>
	>
	void hartman_path_list_color_recursive(const index_graph & graph, const planar_embedding & embedding, 
		color_list_map & color_list, list_color_property_map & properties, disjoint_set & face_location_sets,
		color_map & coloring, vertex_descriptor x, const edge_iterator_pair & x_range, vertex_descriptor y,
		const edge_iterator_pair & y_range, vertex_descriptor p, int before_p, int before_y, int before_x)
	{
		// Update potentially conditional ranges for x and y
		properties[x].set_current_range(x_range);
		
		if(face_location_sets.exists(before_p)) {
			before_p = properties[x].set_face_location(before_p, face_location_sets);
		}
		
		if(x != y) {
			properties[y].set_current_range(y_range);
			before_x = properties[y].set_face_location(before_x, face_location_sets);
		}
		
		// Base Case 2: K_2
		if(properties[p].single_neighbor()) {
			// Color other vertex if not colored
			vertex_descriptor neighbor =
				get_incident_vertex(p, *properties[p].get_current_range().first, graph);
			
			// If there is no colored path, handle color both vertices from their remaining lists
			if(!properties[p].colored()) {
				properties[p].color();
				coloring[p] = color_list[p].front();
				
				if(!properties[neighbor].colored()) {
					properties[neighbor].color();
					coloring[neighbor] = color_list[neighbor].front();
				}
			}
			// Else if the last vertex is not part of the path, ensure it recieves a different color
			else if(!properties[neighbor].colored()) {
				color_list[neighbor].erase(
						std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), coloring[p]),
						color_list[neighbor].end()
					);
				
				properties[neighbor].color();
				coloring[neighbor] = color_list[neighbor].front();
			}
			
			return;
		}
		
		// If colored path doesn't exist, we create one
		if(!properties[p].colored()) {
			// Begin coloring path with first color in x's list
			auto path_color = color_list[x].front();
			
			// Must track current path end and next path vertex starting with x
			vertex_descriptor path_end = x, next_vertex = x;
			coloring[x] = path_color;
			properties[x].color();
			
			// Reassign p to our new path start
			p = x;
			
			if(x != y) {
				do {
					// Update the path end
					path_end = next_vertex;
				
					// Look counterclockwise through our current range of interior neighbors
					edge_iterator_pair range = properties[path_end].get_range();
					edge_iterator_pair current_range = properties[path_end].get_current_range();
					for(auto edge_iter = current_range.first; edge_iter != current_range.second; ++edge_iter) {
						if(edge_iter == range.second) edge_iter = range.first;
						vertex_descriptor neighbor = get_incident_vertex(path_end, *edge_iter, graph);
					
						// Check if vertex is on the outer face between x and y
						if(properties[neighbor].on_face() && (neighbor == y ||
							face_location_sets.compare(properties[neighbor].get_face_location(), before_y))) {
							// Check if it may be colored path_color
							for(auto color : color_list[neighbor]) {
								if(color == path_color) {
									// Color and append the vertex to the path
									next_vertex = neighbor;
									coloring[next_vertex] = path_color;
									properties[next_vertex].color();
									
									break;
								}
							}
						}
						else if(neighbor == y && properties[y].colored() && coloring[y] == path_color) {
							next_vertex = y;
						}
						
						// If we have found a new end vertex, stop looking and see if we made a lobe
						if(path_end != next_vertex) {
							// Grab current range of new path vertex
							edge_iterator_pair n_current_range = properties[next_vertex].get_current_range();
					
							// Check if we took a chord
							vertex_descriptor prev_on_face =
								get_incident_vertex(next_vertex, *(n_current_range.first), graph);
							if(!properties[prev_on_face].colored()) {
								// Split range at last path vertex
								auto nv_ranges = properties[next_vertex].split_range_at(path_end);
								auto pe_ranges = properties[path_end].split_range_at(next_vertex);
							
								// Color lobe
								hartman_path_list_color_recursive(graph, embedding, color_list, properties,
									face_location_sets, coloring, next_vertex, nv_ranges.first,
									path_end, pe_ranges.second, next_vertex, -1, -1, before_y);
						
								properties[next_vertex].set_current_range(nv_ranges.second);
								properties[path_end].set_current_range(pe_ranges.first);
							}
							
							break;
						}
					}
				} while(path_end != next_vertex && next_vertex != y);
			}
			
			hartman_path_list_color_recursive(graph, embedding, color_list, properties, face_location_sets,
				coloring, x, properties[x].get_current_range(), y,
				properties[y].get_current_range(), p, -1, before_y, before_x);
			
			return;
		}
		
		// New parameters for coloring the graph after removing p
		vertex_descriptor new_x = x, new_y = y;
		
		// Remeber iterator to the last edge in incidence range
		edge_iterator pre_end = --(properties[p].get_current_range().second);
		
		// Iterate counterclockwise through interior neighbors of p
		edge_iterator_pair range = properties[p].get_range();
		edge_iterator_pair current_range = properties[p].get_current_range();
		for(auto edge_iter = current_range.first; edge_iter != current_range.second; ++edge_iter) {
			if(edge_iter == range.second) edge_iter = range.first;
			vertex_descriptor neighbor = get_incident_vertex(p, *edge_iter, graph);
			
			// Remove p's color from the list of all adjacent neighbors
			color_list[neighbor].erase(
					std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), coloring[p]),
					color_list[neighbor].end()
				);
			
			if(properties[neighbor].interior()) {
				// Add neighbor to the outer face with incidence list starting at p
				before_p = properties[neighbor].initialize(neighbor, p, before_p, graph, embedding,
					face_location_sets);
				
				// Remove p from incidence list
				properties[neighbor].remove_begin();
			}
			else {
				// Vertex prior to p on the outer face
				if(edge_iter == current_range.first) {
					// Reassign x or y as needed
					if(x == p) {
						if(y == p) {
							before_x = properties[neighbor].set_face_location(before_x, face_location_sets);
							new_y = neighbor;
						}
						else {
							before_p = properties[neighbor].set_face_location(before_p, face_location_sets);
							new_x = neighbor;
						}
					}
					
					// Remove p from incidence list
					properties[neighbor].remove_end();
				}
				else {
					int neighbor_face_location = properties[neighbor].get_face_location();
					before_p = properties[neighbor].set_face_location(before_p, face_location_sets);
					
					// Reassign x or y as needed
					if(y == p) {
						if(x == p) {
							new_x = neighbor;
						}
						else {
							new_y = neighbor;
						}
					}
					
					if(edge_iter == pre_end) {
						// Remove p from incidence list
						properties[neighbor].remove_begin();
						vertex_descriptor new_p = neighbor;
						
						// Correction for weird case where x = y = p
						if(x == p && y == p) {
							before_x = before_p;
							before_p = -1;
						}
						// If we are removing the last path vertex, we must merge before_p and before_y
						else if(!properties[neighbor].colored() || coloring[neighbor] != coloring[p]) {
							before_y = face_location_sets.take_union(before_p, before_y);
							before_p = -1;
							new_p = new_x;
						}
						
						// Color graph with p removed
						hartman_path_list_color_recursive(graph, embedding, color_list, properties,
							face_location_sets, coloring, new_x, properties[new_x].get_current_range(), new_y,
							properties[new_y].get_current_range(), new_p, before_p, before_y, before_x);
					}
					else {
						// Split ranges along edge
						auto n_ranges = properties[neighbor].split_range_at(p);
						auto p_ranges = properties[p].split_range_at(neighbor);
					
						// Special case for x = y = p
						if(x == p && y == p) {
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
							properties[p].set_current_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, new_x,
								properties[new_x].get_current_range(), new_y,
								properties[new_y].get_current_range(), new_x, -1, before_y, before_p);
							
							// Reassign ranges for remaining "right" subgraph
							properties[neighbor].set_current_range(n_ranges.first);
							
							// Color remaining "right" subgraph
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, p, p_ranges.second, p, p_ranges.second, p,
								-1, before_y, -1);
						}
						// Case 1: C[x,p) Cutvertex between x and p inclusive
						else if(face_location_sets.compare(neighbor_face_location, before_p)) {
							// Set ranges for "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							properties[neighbor].set_current_range(n_ranges.first);
							
							// Color "right" subgraph
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, x,
								properties[x].get_current_range(), y, properties[y].get_current_range(), p,
								before_p, before_y, before_x);
							
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
							properties[p].set_current_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, neighbor,
								properties[neighbor].get_current_range(), neighbor,
								properties[neighbor].get_current_range(), neighbor, -1, before_p, -1);
						}
						// Case 2: C(p,y] Cutvertex between p and y inclusive
						else if(face_location_sets.compare(neighbor_face_location, before_y)) {
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
							properties[p].set_current_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Intersect face sections we are joining
							before_y = face_location_sets.take_union(before_p, before_y);
							before_p = -1;
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, new_x, properties[new_x].get_current_range(),
								new_y, properties[new_y].get_current_range(), new_x, -1, before_y, before_x);
								
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							
							// Color remaining "right" subgraph
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, neighbor, n_ranges.first, neighbor,
								n_ranges.first, neighbor, -1, before_y, -1);
						}
						//Case 3: C(y,x) Cutvertex between y and x exclusive
						else {
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							
							// Color remaining "right" subgraph
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, neighbor, n_ranges.first, y,
								properties[y].get_current_range(), p, -1, before_y, before_x);
							
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
							properties[p].set_current_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Reassign p
							vertex_descriptor new_p = new_x;
							if(neighbor == y && properties[p].colored() && coloring[p] == coloring[neighbor]) {
								new_p = y;
							}
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, new_x, properties[new_x].get_current_range(),
								neighbor, properties[neighbor].get_current_range(), new_p, -1, before_p,
								before_x);
						}
					}
					
					break;
				}
			}
		}
	}
}

#endif