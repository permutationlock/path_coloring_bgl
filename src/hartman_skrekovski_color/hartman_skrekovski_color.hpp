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
#include "augmented_embedding.hpp"

namespace boost {
	 
	// Stores all vertex propreties for path 3-list-coloring algorithm
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
			bool single_neighbor() const {
				return range.first == range.second;
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
		typedef typename graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
		
		// Construct our special embedded adjacency list with backwards lookup
		augmented_embedding<index_graph> plane_graph(graph, embedding);
		
		// Define the type for our property map
		typedef iterator_property_map<
				typename std::vector<list_color_properties<index_graph> >::iterator,
				typename property_map<index_graph, vertex_index_t>::const_type
			> list_color_property_map;
		
		// Make property map to store all properties
		std::vector<list_color_properties<index_graph> > list_color_property_storage(num_vertices(graph));
		list_color_property_map properties(list_color_property_storage.begin(), get(vertex_index, graph));
		
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
				    std::cout << "initializing " << vertex << " from " << last_vertex << "\n";
					// Initialize properties of the next vertex on the face using the current face edge
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
				x, properties[x].get_range(),
				y, properties[y].get_range(),
				x,
				-1, before_y, -1
			);
	}
	
	template<
			typename index_graph, typename list_color_property_map,
			typename color_list_map, typename color_map,
			typename neighbor_range,
			typename vertex_descriptor = typename augmented_embedding<index_graph>::vertex_descriptor
		>
	void hartman_path_list_color_recursive(
			const augmented_embedding<index_graph> & plane_graph,
			list_color_property_map & properties, disjoint_set & face_locations,
			color_list_map & color_list, color_map & coloring,
			vertex_descriptor x, const neighbor_range & x_range,
			vertex_descriptor y, const neighbor_range & y_range,
			vertex_descriptor p,
			int before_p, int before_y, int before_x
		)
	{
	    std::cout << "x = " << x << ", y = " << y << ", p = " << p << "\n";
		// Update potentially conditional ranges for x and y
		properties[x].set_range(x_range);
		if(x != y) properties[y].set_range(y_range);
		
		std::cout << "x_range = (" << plane_graph[x][properties[x].get_range().first].neighbor
		    << ", " << plane_graph[x][properties[x].get_range().second].neighbor << ")\n";
	    std::cout << "y_range = (" << plane_graph[y][properties[y].get_range().first].neighbor
		    << ", " << plane_graph[y][properties[y].get_range().second].neighbor << ")\n";
		
		// If there is a before_p section, add x to it
		if(face_locations.exists(before_p)) {
			before_p = properties[x].set_face_location(before_p, face_locations);
		}
		
		// Base Case 2: K_2
		if(properties[p].single_neighbor()) {
			vertex_descriptor neighbor = plane_graph[p][properties[p].get_range().first].neighbor;
			
			// If there is no colored path, color both vertices from their remaining lists
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
		    std::cout << "coloring path\n";
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
					neighbor_range range(properties[path_end].get_range());
					std::size_t current_index = range.first;
					do {
					    current_index %= plane_graph[path_end].size();
						vertex_descriptor neighbor = plane_graph[path_end][current_index].neighbor;
						std::size_t back_index = plane_graph[path_end][current_index].back_index;
					
					    std::cout << "looking at " << neighbor << "\n";
					    
						// Check if vertex is on the outer face between x and y
						if(properties[neighbor].on_face() && (neighbor == y ||
							face_locations.compare(properties[neighbor].get_face_location(), before_y)))
						{
							// Check if it may be colored path_color
							for(auto color : color_list[neighbor]) {
								if(color == path_color) {
								    std::cout << "coloring " << neighbor << "\n";
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
							neighbor_range n_range(properties[next_vertex].get_range());
					
							// Check if we took a chord
							vertex_descriptor prev_on_face = plane_graph[next_vertex][n_range.first].neighbor;
							if(!properties[prev_on_face].colored()) {
								// Split range at last path vertex
								auto pe_ranges = properties[path_end].split_range(current_index);
								auto nv_ranges = properties[next_vertex].split_range(back_index);
								
								// Color lobe
								hartman_path_list_color_recursive(
										plane_graph,
										properties, face_locations,
										color_list, coloring,
										next_vertex, nv_ranges.first,
										path_end, pe_ranges.second,
										next_vertex,
										-1, -1, before_y
									);
						
								properties[next_vertex].set_range(nv_ranges.second);
								properties[path_end].set_range(pe_ranges.first);
							}
							
							break;
						}
					} while(current_index++ != range.second);
				} while(path_end != next_vertex && next_vertex != y);
			}
			
			hartman_path_list_color_recursive(
					plane_graph,
					properties, face_locations,
					color_list, coloring,
					x, properties[x].get_range(),
					y, properties[y].get_range(),
					p,
					-1, before_y, before_x
				);
			
			return;
		}
		
		// New parameters for coloring the graph after removing p
		vertex_descriptor new_x = x, new_y = y;
		
		std::cout << "removing " << p << "\n";
		
		// Iterate counterclockwise through interior neighbors of p
		neighbor_range range(properties[p].get_range());
		std::size_t current_index = range.first;
		do {
		    current_index %= plane_graph[p].size();
			vertex_descriptor neighbor = plane_graph[p][current_index].neighbor;
			std::size_t back_index = plane_graph[p][current_index].back_index;
			
			std::cout << "looking at " << neighbor << "\n";
			
			// Remove p's color from the list of all adjacent neighbors
			color_list[neighbor].erase(
					std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), coloring[p]),
					color_list[neighbor].end()
				);
			
			if(properties[neighbor].interior()) {
			    std::cout << "interior\n";
				// Add neighbor to the outer face with incidence list starting at p
				before_p = properties[neighbor].initialize(
				        back_index, plane_graph[neighbor].size(),
				        before_p, face_locations
			        );
				
				// Remove p from incidence list
				properties[neighbor].remove_begin();
			}
			else {
				// Vertex directly prior to p on the outer face
				if(current_index == range.first) {
				    std::cout << "begin\n";
					// Reassign x or y as needed
					if(x == p) {
						if(y == p) {
							before_x = properties[neighbor].set_face_location(before_x, face_locations);
							new_y = neighbor;
						}
						else {
							before_p = properties[neighbor].set_face_location(before_p, face_locations);
							new_x = neighbor;
						}
					}
					
					// Remove p from incidence list
					properties[neighbor].remove_end();
				}
				else {
					int neighbor_face_location = properties[neighbor].get_face_location();
					before_p = properties[neighbor].set_face_location(before_p, face_locations);
					
					// Reassign x or y as needed if they are being removed
					if(y == p) {
						if(x == p) {
							new_x = neighbor;
						}
						else {
							new_y = neighbor;
						}
					}
					
					// Vertex directly following p on the outer face
					if(current_index == range.second) {
					    std::cout << "end\n";
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
							before_y = face_locations.take_union(before_p, before_y);
							before_p = -1;
							new_p = new_x;
						}
						
						// Color graph with p removed
						hartman_path_list_color_recursive(
								plane_graph,
								properties, face_locations,
								color_list, coloring,
								new_x, properties[new_x].get_range(),
								new_y, properties[new_y].get_range(),
								new_p,
								before_p, before_y, before_x
							);
					}
					else {
						// Split ranges along edge
						auto p_ranges = properties[p].split_range(current_index);
						auto n_ranges = properties[neighbor].split_range(back_index);
					
						// Special case for x = y = p
						if(x == p && y == p) {
						    std::cout << "x = y\n";
							// Set ranges for "left" half to color
							properties[neighbor].set_range(n_ranges.second);
							properties[p].set_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									new_x, properties[new_x].get_range(),
									new_y, properties[new_y].get_range(),
									new_x,
									-1, before_y, before_p
								);
							
							// Reassign ranges for remaining "right" subgraph
							properties[neighbor].set_range(n_ranges.first);
							
							// Color remaining "right" subgraph
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									p, p_ranges.second,
									p, p_ranges.second,
									p,
									-1, before_y, -1
								);
						}
						//Case 1: C[y,x) Cutvertex between y and x, including y
						else if(neighbor == y || face_locations.compare(neighbor_face_location, before_x)) {
						    std::cout << "before_x\n";
						    
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_range(p_ranges.second);
							
							// Color remaining "right" subgraph
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									neighbor, n_ranges.first,
									y, properties[y].get_range(),
									p,
									-1, before_y, before_x
								);
							
							// Set ranges for "left" half to color
							properties[neighbor].set_range(n_ranges.second);
							properties[p].set_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Reassign p
							vertex_descriptor new_p = new_x;
							if(neighbor == y && properties[p].colored() && coloring[p] == coloring[neighbor]) {
								new_p = y;
							}
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									new_x, properties[new_x].get_range(),
									neighbor, properties[neighbor].get_range(),
									new_p,
									-1, before_p, before_x
								);
						}
						// Case 2: C[x,p) Cutvertex between x and p, including x
						else if(face_locations.compare(neighbor_face_location, before_p)) {
						    std::cout << "before_p\n";
							// Set ranges for "right" subgraph
							properties[p].set_range(p_ranges.second);
							properties[neighbor].set_range(n_ranges.first);
							
							// Color "right" subgraph
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									x, properties[x].get_range(),
									y, properties[y].get_range(),
									p,
									before_p, before_y, before_x
								);
							
							// Set ranges for "left" half to color
							properties[neighbor].set_range(n_ranges.second);
							properties[p].set_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list,  coloring,
									neighbor, properties[neighbor].get_range(),
									neighbor, properties[neighbor].get_range(),
									neighbor,
									-1, before_p, -1
								);
						}
						// Case 3: C(p,y) Cutvertex between p and y, exclusive
						else {
						    std::cout << "before_y\n";
							// Set ranges for "left" half to color
							properties[neighbor].set_range(n_ranges.second);
							properties[p].set_range(p_ranges.first);
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Intersect face sections we are joining
							before_y = face_locations.take_union(before_p, before_y);
							
							// Color "left" subgraph with p removed
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									new_x, properties[new_x].get_range(),
									new_y, properties[new_y].get_range(),
									new_x,
									-1, before_y, before_x
								);
								
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_range(p_ranges.second);
							
							// Color remaining "right" subgraph
							hartman_path_list_color_recursive(
									plane_graph,
									properties, face_locations,
									color_list, coloring,
									neighbor, n_ranges.first,
									neighbor, n_ranges.first,
									neighbor,
									-1, before_y, -1
								);
						}
					}
					
					break;
				}
			}
		} while(current_index++ != range.second);
	}
}

#endif