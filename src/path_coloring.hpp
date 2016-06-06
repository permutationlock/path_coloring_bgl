/*
 * path_coloring.h
 * Author: Aven Bross
 * 
 * Implementation of several algorithms for path coloring graphs.
 */

#ifndef __PATH_COLORING_HPP
#define __PATH_COLORING_HPP

#include <vector>
#include <list>
#include <unordered_map>
#include <queue>
#include <stdexcept>
#include <utility>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include "draw_tikz_graph.hpp"
#include "disjoint_set.hpp"

//#define PLC_SHOW_ANNOTATIONS

namespace boost {
	
	/*
	 * get_incident_vertex
	 *
	 * Preconditions: edge is incident to vertex in graph.
	 * 
	 * Postconditions: returns the endpoint of edge that is not vertex.
	 */
	template<typename Graph>
	typename graph_traits<Graph>::vertex_descriptor
		get_incident_vertex(typename graph_traits<Graph>::vertex_descriptor vertex,
		typename graph_traits<Graph>::edge_descriptor edge, const Graph & graph)
	{
		auto adj = source(edge, graph);
		if(adj == vertex)
			adj = target(edge, graph);
		return adj;
	}
	
	/*
	 * poh_path_color
	 *
	 * Preconditions: graph is a weakly triangulated planar graph the incident edges of each
	 * vertex arranged in counterclockwise order in the embedding. p and q are paths with
	 * p above q (in counterclockwise orientation scheme) such that joining the beginnings
	 * and ends of the paths forms a chordless cycle, and with p and q each colored a
	 * different color. new_color is a third color different from the color of p and q.
	 * 
	 * Postconditions: a coloring of g with the face enclosed by the cycle pq
	 * 3-colored such that each color class induces a disjoint union of paths.
	 */
	template<typename Graph, typename Embedding, typename Coloring, typename VertexIter,
		typename color_type = typename property_traits<Coloring>::value_type>
	void poh_path_color(const Graph & graph, const Embedding & embedding,
		VertexIter p_begin, VertexIter p_end, VertexIter q_begin, VertexIter q_end, Coloring & coloring,
		color_type new_color)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		
		// Base case: Paths are each single vertices, K_2
		if(p_end - p_begin == 1 && q_end - q_begin == 1)
		{
			return;
		}
		
		// Vertices at the start and end of our paths
		vertex_descriptor p_0 = *p_begin, p_1 = *(p_end - 1), q_0 = *q_begin,
			q_1 = *(q_end - 1);
			
		//std::cout << "p_0=" << p_0 << ", p_1=" << p_1 << ", q_0=" << q_0 << ", q_1=" << q_1 << "\n";
		
		// Triangulation vertices at the beginning and end of the paths
		vertex_descriptor t_0, t_1;
		
		// Case 1.1: Find t_0 and check if it hits a path
		{
			t_0 = q_0;
			
			// Iterate q_0's incident edges in counterclockwise order
			auto ordering = embedding[q_0];
			for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); ++edge_iter)
			{
				// Find p_0 for orientation
				if(get_incident_vertex(q_0, *edge_iter, graph) == p_0)
				{
					// t_0 is one back from p_0 in the adjacency list
					if(edge_iter == ordering.begin())
					{
						t_0 = get_incident_vertex(q_0, *(ordering.end() - 1), graph);
					}
					else
					{
						t_0 = get_incident_vertex(q_0, *(edge_iter - 1), graph);
					}
					break;
				}
			}
			
			if(t_0 == q_0) throw std::runtime_error("Invalid embedding (no t_0).");
			
			// Case 1.1.1: Triangle formed at start with path p
			if(p_end - p_begin > 1 && t_0 == *(p_begin + 1))
			{
				poh_path_color(graph, embedding, p_begin + 1, p_end, q_begin, q_end, coloring, new_color);
				return;
			}
			
			// Case 1.1.2: Triangle formed at start with path q
			if(q_end - q_begin > 1 && t_0 == *(q_begin + 1))
			{
				poh_path_color(graph, embedding, p_begin, p_end, q_begin + 1, q_end, coloring, new_color);
				return;
			}
		}
		
		// Case 1.2: Find t_1 and check if it hits a path
		{
			t_1 = p_1;
			
			// Iterate p_1's incident edges in counterclockwise order
			auto ordering = embedding[p_1];
			for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); ++edge_iter)
			{
				// Find q_1 for orientation
				if(get_incident_vertex(p_1, *edge_iter, graph) == q_1)
				{
					// t_1 is one back from q_1 in the adjacency list
					if(edge_iter == ordering.begin())
					{
						t_1 = get_incident_vertex(p_1, *(ordering.end() - 1), graph);
					}
					else
					{
						t_1 = get_incident_vertex(p_1, *(edge_iter - 1), graph);
					}
					break;
				}
			}
			
			if(t_1 == p_1) throw std::runtime_error("Invalid embedding (no t_1).");
			
			// Case 1.2.1: Triangle formed at end with path p
			if(p_end - p_begin > 1 && t_1 == *(p_end - 2))
			{
				//std::cout << "Triangle formed at end with path p.\n";
				poh_path_color(graph, embedding, p_begin, p_end - 1, q_begin, q_end,
					coloring, new_color);
				return;
			}
			
			// Case 1.2.1: Triangle formed at end with path q
			if(q_end - q_begin > 1 && t_1 == *(q_end - 2))
			{
				//std::cout << "Triangle formed at end with path q.\n";
				poh_path_color(graph, embedding, p_begin, p_end, q_begin, q_end - 1,
					coloring, new_color);
				return;
			}
		}
		
		// Case 2: Check for bridging edge between the two paths
		{
			std::unordered_map<vertex_descriptor, VertexIter> p_map;
			
			// Mark all vertices in path p
			for(auto p_iter = p_begin; p_iter != p_end; ++p_iter)
			{
				p_map[*p_iter] = p_iter;
			}
			
			// Check vertices in path q for neighbors in p
			for(auto q_iter = q_begin; q_iter != q_end; ++q_iter)
			{
				// Look at each neighbor
				auto ordering = embedding[*q_iter];
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); ++edge_iter)
				{
					vertex_descriptor neighbor = get_incident_vertex(*q_iter, *edge_iter, graph);
					
					// If a vertex has neighbor in p
					if(p_map.count(neighbor) != 0)
					{
						// Check if our edge is one of the outside edges
						if(*q_iter == q_0 && *p_map[neighbor] == p_0)
							continue;
						if(*q_iter == q_1 && *p_map[neighbor] == p_1)
							continue;
						
						// Recurse on left and right halves
						poh_path_color(graph, embedding, p_begin, p_map[neighbor] + 1,
							q_begin, q_iter + 1, coloring, new_color);
						poh_path_color(graph, embedding, p_map[neighbor], p_end,
							q_iter, q_end, coloring, new_color);
						return;
					}
				}
			}
		}
		
		// Case 3: Find splitting path between triangulation vertices
		{
			/* 
			 * We will perform a BFS from t_1 to find a chordless t_0t_1-path.
			 * We initially mark paths p and q to contain the search within the
			 * pq-cycle.
			 */
			
			std::unordered_map<vertex_descriptor, vertex_descriptor> parent_map;
			
			// Mark vertices in path p
			for(auto p_iter = p_begin; p_iter != p_end; ++p_iter)
			{
				parent_map[*p_iter] = *p_iter;
			}
			
			// Mark vertices in path q
			for(auto q_iter = q_begin; q_iter != q_end; ++q_iter)
			{
				parent_map[*q_iter] = *q_iter;
			}
			
			std::queue<vertex_descriptor> bfs_queue;
			
			// Start bfs
			parent_map[t_1] = t_1;
			bfs_queue.push(t_1);
			
			vertex_descriptor curr_vertex = t_1;
			
			while(!bfs_queue.empty())
			{
				curr_vertex = bfs_queue.front();
				bfs_queue.pop();
				
				// We are done as soon as we find t_0
				if(curr_vertex == t_0) break;
				
				// Add all unmarked neighbors of curr_vertex to bfs_queue
				auto ordering = embedding[curr_vertex];
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); ++edge_iter)
				{
					vertex_descriptor neighbor = get_incident_vertex(curr_vertex, *edge_iter, graph);
					
					if(parent_map.count(neighbor) == 0)
					{
						parent_map[neighbor] = curr_vertex;
						bfs_queue.push(neighbor);
					}
				}
			}
			
			if(curr_vertex != t_0) throw std::runtime_error("Invalid embedding (no splitting path).");
			
			// Backtrack and generate splitting path, coloring it new_color
			std::vector<vertex_descriptor> splitting_path;
			splitting_path.push_back(curr_vertex);
			coloring[curr_vertex] = new_color;
			
			while(parent_map[curr_vertex] != curr_vertex)
			{
				curr_vertex = parent_map[curr_vertex];
				splitting_path.push_back(curr_vertex);
				coloring[curr_vertex] = new_color;
			}
			
			// Recurse on top and bottom halves with the appropriate color
			poh_path_color(graph, embedding, p_begin, p_end, splitting_path.begin(),
				splitting_path.end(), coloring, coloring[q_0]);
			poh_path_color(graph, embedding, splitting_path.begin(), splitting_path.end(),
				q_begin, q_end, coloring, coloring[p_0]);
		}
	}
	
	/*
	 * hartman_path_list_color_biconnected
	 *
	 * Preconditions: graph is a weakly triangulated planar biconnected graph. Incidence lists of each
	 * vertex are arranged in counterclockwise order in the embedding. p and q are paths with
	 * p above q (in counterclockwise orientation scheme) such that joining the beginnings
	 * and ends of the paths forms a cycle. The first vertex of p is our vertex x
	 * and the first vertex of q is our vertex y. The property_map color_list maps each
	 * vertex to a list of potential colors. Vertices x and y have a list size >=1,
	 * vertices in p or q have list size >=2, and all other vertices have list size
	 * >= 3.
	 * 
	 * Postconditions: coloring maps each vertex of g contained in the pq-cycle,
	 * including vertices in p and q, to a color in their respective color_list
	 * and each color class induces a disjoint union of cycles. Furthermore, x and
	 * y have no more than one neighbor assigned the same color.
	 */
	 
	// Stores all vertex propreties for path 3-list-coloring algorithm
	template<typename vertex_descriptor, typename edge_iterator>
	class list_color_properties
	{
	public:
		typedef typename std::pair<edge_iterator, edge_iterator> edge_iterator_pair;
		enum vertex_state { INTERIOR, ON_FACE, COLORED };

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
			for(edge_iterator edge_iter = range.first; edge_iter != range.second; ++edge_iter)
			{
				vertex_descriptor neighbor = get_incident_vertex(vertex, *edge_iter, graph);
		
				if(neighbor == adding_vertex)
				{
					// Orient our incidence list to start at the adding vertex
					current_range.first = edge_iter;
					current_range.second = edge_iter;
			
					// Correction to make sure list wraps correctly
					if(current_range.second == range.first)
						current_range.second = range.second;
				}
		
				// Remember incidence list location of each neighbor
				neighbor_iterator[neighbor] = edge_iter;
			}
			
			return face_location;
		}
		void color()
		{
			state = COLORED;
		}
		vertex_state get_state() const
		{
			return state;
		}
		int set_face_location(int new_face_location, disjoint_set & face_location_sets)
		{
			if(!face_location_sets.exists(new_face_location)) face_location = face_location_sets.make_next();
			else face_location = face_location_sets.find(new_face_location);
			return face_location;
		}
		int get_face_location()
		{
			return face_location;
		}
		void set_current_range(const edge_iterator_pair & new_range)
		{
			current_range = new_range;
		}
		void remove_begin()
		{
			++current_range.first;
			if(current_range.first == range.second && current_range.first != current_range.second)
				current_range.first = range.first;
		}
		void remove_end()
		{
			--current_range.second;
			if(current_range.second == range.first && current_range.first != current_range.second)
				current_range.second = range.second;
		}
		std::pair<edge_iterator_pair, edge_iterator_pair> split_range_at(vertex_descriptor neighbor) const
		{
			edge_iterator split_iter = neighbor_iterator.at(neighbor);
	
			// Range after split edge
			edge_iterator_pair second_range(split_iter, current_range.second);
	
			// Range before split edge
			edge_iterator_pair first_range(current_range.first, ++split_iter);
	
			return std::pair<edge_iterator_pair, edge_iterator_pair>(first_range, second_range);
		}
		edge_iterator_pair get_range() const
		{
			return range;
		}
		edge_iterator_pair get_current_range() const
		{
			return current_range;
		}
		bool no_neighbors() const
		{
			return current_range.first == current_range.second;
		}
		bool single_neighbor() const
		{
			edge_iterator temp_first = current_range.first;
			return ++temp_first == current_range.second;
		}
		bool interior() const
		{
			return state == INTERIOR;
		}
		bool on_face() const
		{
			return state == ON_FACE;
		}
		bool colored() const
		{
			return state == COLORED;
		}

	private:
		// Face location union find data
		int face_location;
		
		vertex_state state;
		edge_iterator_pair range;	// Incidence range from embedding
		edge_iterator_pair current_range;	// Current incidence range in embedding

		// Maps each neighbor to its corresponding edge in embedding
		std::unordered_map<vertex_descriptor, edge_iterator> neighbor_iterator;
	};
	 
	 template<typename index_graph, typename planar_embedding, typename color_list_map, typename color_map,
	 	typename vertex_iterator>
	void path_list_color(const index_graph & graph, const planar_embedding & embedding,
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
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "INTIALIZE: ";
				print_incidence_range(*next, graph, embedding, properties);
			#endif
		}
		
		// Star with x and y as the beginning and end of the outer face list given
		vertex_descriptor x = *face_begin, y = *(--face_end);
		
		// Recursively list color the graph
		path_list_color_recursive(graph, embedding, color_list, properties, face_location_sets, coloring,
			x, properties[x].get_current_range(), y, properties[y].get_current_range(), x,
			-1, before_y, -1);
	}
	
	bool compare_face_locations(int first_location, int second_location, disjoint_set & face_location_sets)
	{
		if(!face_location_sets.exists(second_location)) return false;
		return (face_location_sets.find(first_location) == face_location_sets.find(second_location));
	}
	
	template<
		typename index_graph, typename planar_embedding, typename color_list_map,
		typename color_map, typename list_color_property_map,
		typename edge_iterator = typename property_traits<planar_embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<index_graph>::vertex_descriptor,
		typename edge_iterator_pair = typename std::pair<edge_iterator, edge_iterator>
	>
	void path_list_color_recursive(const index_graph & graph, const planar_embedding & embedding, 
		color_list_map & color_list, list_color_property_map & properties, disjoint_set & face_location_sets,
		color_map & coloring, vertex_descriptor x, edge_iterator_pair x_range, vertex_descriptor y,
		edge_iterator_pair y_range, vertex_descriptor p, int before_p, int before_y, int before_x)
	{
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "x = " << x << ", y = " << y << ", p = " << p << "\n";
			std::cout << "before_p = " << before_p << ", before_y = " << before_y << ", before_x = "
				<< before_x << "\n";
		#endif
		
		// Update potentially conditional ranges for x and y
		properties[x].set_current_range(x_range);
		
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\t";
			print_incidence_range(x, graph, embedding, properties);
		#endif
		
		properties[y].set_current_range(y_range);
		
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\t";
			print_incidence_range(y, graph, embedding, properties);
		#endif
		
		// Base Case 2: K_2
		if(properties[x].single_neighbor())
		{
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tbase_case K_2\n";
			#endif
			
			// Color x if not colored
			if(!properties[x].colored())
			{
				properties[x].color();
				coloring[x] = color_list[x].front();
			}
			
			// Color other vertex if not colored
			vertex_descriptor neighbor =
				get_incident_vertex(x, *properties[x].get_current_range().first, graph);
			
			if(!properties[neighbor].colored())
			{
				// If neighbor is not a not a special vertex, ensure it recieves a different color
				if(x == y)
				{
					std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), coloring[x]);
				}
				properties[neighbor].color();
				coloring[neighbor] = color_list[neighbor].front();
			}
			return;
		}
		
		// If colored path doesn't exist, we create one
		if(!properties[p].colored())
		{
			// Begin coloring path with first color in x's list
			auto path_color = color_list[x].front();
			
			// Must track current path end and next path vertex starting with x
			vertex_descriptor path_end = x, next_vertex = x;
			coloring[x] = path_color;
			properties[x].color();
			
			// Reassign p to our new path start
			p = x;
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tcoloring path with path_color = " << path_color << "\n";
			#endif
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\t\t\tcolored vertex = " << x << "\n";
			#endif
			
			if(x != y)
			{
				do
				{
					// Update the path end
					path_end = next_vertex;
				
					// Look counterclockwise through our current range of interior neighbors
					edge_iterator_pair range = properties[path_end].get_range();
					edge_iterator_pair current_range = properties[path_end].get_current_range();
					for(auto edge_iter = current_range.first; edge_iter != current_range.second; ++edge_iter)
					{
						if(edge_iter == range.second) edge_iter = range.first;
						vertex_descriptor neighbor = get_incident_vertex(path_end, *edge_iter, graph);
					
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\tchecking state of vertex " << neighbor << "\n";
						#endif
					
						// Check if vertex is on the outer face between x and y
						if(neighbor == y || (properties[neighbor].on_face() &&
							compare_face_locations(properties[neighbor].get_face_location(),
								before_y, face_location_sets)))
						{
							// Check if it may be colored path_color
							for(auto color : color_list[neighbor])
							{
								if(color == path_color)
								{
									// Color and append the vertex to the path
									next_vertex = neighbor;
									coloring[next_vertex] = path_color;
									properties[next_vertex].color();
								
									#ifdef PLC_SHOW_ANNOTATIONS
										std::cout << "\t\t\tcolored vertex = " << next_vertex << "\n";
									#endif
								
									// Grab current range of new path vertex
									edge_iterator_pair n_current_range = properties[next_vertex].get_current_range();
								
									// Check if we took a chord
									vertex_descriptor prev_on_face =
										get_incident_vertex(next_vertex, *(n_current_range.first), graph);
									if(!properties[prev_on_face].colored())
									{
										// Split range at last path vertex
										auto nv_ranges = properties[next_vertex].split_range_at(path_end);
										auto pe_ranges = properties[path_end].split_range_at(next_vertex);
										
										#ifdef PLC_SHOW_ANNOTATIONS
											std::cout << "\t\t\t\tcoloring lobe\n";
										#endif
										
										// Color lobe
										path_list_color_recursive(graph, embedding, color_list, properties,
											face_location_sets, coloring, next_vertex, nv_ranges.first,
											path_end, pe_ranges.second, next_vertex, -1, -1, before_y);
									
										properties[next_vertex].set_current_range(nv_ranges.second);
									
										#ifdef PLC_SHOW_ANNOTATIONS
											std::cout << "\t\t\t\t";
											print_incidence_range(next_vertex, graph, embedding, properties);
										#endif
									
										properties[path_end].set_current_range(pe_ranges.first);
									
										#ifdef PLC_SHOW_ANNOTATIONS
											std::cout << "\t\t\t\t";
											print_incidence_range(path_end, graph, embedding, properties);
										#endif
									}
								
									break;
								}
							}
						}
					
						// If we have found a new end vertex, stop looking at neighbors
						if(path_end != next_vertex || neighbor == y) break;
					}
				}
				while(path_end != next_vertex && next_vertex);
			}
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tcoloring remaining graph\n";
			#endif
			
			path_list_color_recursive(graph, embedding, color_list, properties, face_location_sets,
				coloring, x, properties[x].get_current_range(), y,
				properties[y].get_current_range(), p, -1, before_y, before_x);
			
			return;
		}
		
		// New parameters for coloring the graph after removing p
		vertex_descriptor new_x = x, new_y = y;
		
		// Remeber iterator to the last edge in incidence range
		edge_iterator pre_end = --(properties[p].get_current_range().second);
		
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\tremoving p = " << p << "\n";
		#endif
		
		// Iterate counterclockwise through interior neighbors of p
		edge_iterator_pair range = properties[p].get_range();
		edge_iterator_pair current_range = properties[p].get_current_range();
		for(auto edge_iter = current_range.first; edge_iter != current_range.second; ++edge_iter)
		{
			if(edge_iter == range.second) edge_iter = range.first;
			vertex_descriptor neighbor = get_incident_vertex(p, *edge_iter, graph);
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\t\tlooking at vertex " << neighbor << "\n";
			#endif
			
			// Remove p's color from the list of all adjacent neighbors
			std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), coloring[p]);
			
			// Interior vertex
			if(properties[neighbor].interior())
			{
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\t\tinterior\n";
				#endif
				
				// Add neighbor to the outer face with incidence list starting at p
				before_p = properties[neighbor].initialize(neighbor, p, before_p, graph, embedding,
					face_location_sets);
				
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\t\t\t";
					std::cout << "before_p = " << before_p << "\n";
					std::cout << "\t\t\t\t";
					std::cout << "mark[" << neighbor << "] = " << properties[neighbor].get_face_location() << "\n";
				#endif
				
				// Remove p from incidence list
				properties[neighbor].remove_begin();
				
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\t\t\tINITIALIZE: ";
					print_incidence_range(neighbor, graph, embedding, properties);
				#endif
			}
			// Face vertex
			else
			{
				// Vertex prior to p on the outer face
				if(edge_iter == current_range.first)
				{
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\tbegin\n";
					#endif
					// Reassign x or y as needed
					if(x == p)
					{
						if(y == p)
						{
							before_y = properties[neighbor].set_face_location(before_y, face_location_sets);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "before_y = " << before_y << "\n";
							#endif
							new_y = neighbor;
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "mark[" << new_y << "] = "
									<< properties[new_y].get_face_location() << "\n";
							#endif
						}
						else
						{
							before_p = properties[neighbor].set_face_location(before_p, face_location_sets);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "before_p = " << before_p << "\n";
							#endif
							new_x = neighbor;
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "mark[" << new_x << "] = "
									<< properties[new_x].get_face_location() << "\n";
							#endif
						}
					}
					
					// Remove p from incidence list
					properties[neighbor].remove_end();
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\t\t";
						print_incidence_range(neighbor, graph, embedding, properties);
					#endif
				}
				else
				{
					int neighbor_face_location = properties[neighbor].get_face_location();
					
					before_p = properties[neighbor].set_face_location(before_p, face_location_sets);
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\t\t";
						std::cout << "before_p = " << before_p << "\n";
					#endif
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\t\t";
						std::cout << "mark[" << neighbor << "] = "
							<< properties[neighbor].get_face_location() << "\n";
					#endif
					
					// Reassign x or y as needed
					if(y == p)
					{
						if(x == p)
						{
							new_x = neighbor;
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "mark[" << new_x << "] = "
									<< properties[new_x].get_face_location() << "\n";
							#endif
						}
						else
						{
							before_y = properties[neighbor].set_face_location(before_y, face_location_sets);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "before_y = " << before_y << "\n";
							#endif
							new_y = neighbor;
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "mark[" << new_y << "] = "
									<< properties[new_y].get_face_location() << "\n";
							#endif
						}
					}
					
					if(edge_iter == pre_end)
					{
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\tend\n";
						#endif
					
						// Remove p from incidence list
						properties[neighbor].remove_begin();
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\t\t";
							print_incidence_range(neighbor, graph, embedding, properties);
						#endif
						
						vertex_descriptor new_p = neighbor;
						
						// Correction for weird case where x = y = p
						if(x == y)
						{
							before_x = before_p;
							before_p = -1;
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "before_x = before_p U before_y = " << before_x << "\n";
							#endif
						}
						// If we are removing the last path vertex, we must merge before_p and before_y
						else if(!properties[neighbor].colored())
						{
							before_y = face_location_sets.take_union(before_p, before_y);
							before_p = -1;
							new_p = new_x;
							
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "before_y = before_y U before_p = " << before_y << "\n";
							#endif
						}
						
						// Color graph with p removed
						path_list_color_recursive(graph, embedding, color_list, properties, face_location_sets,
							coloring, new_x, properties[new_x].get_current_range(), new_y,
							properties[new_y].get_current_range(), new_p, before_p, before_y, before_x);
					}
					else
					{
						// Split ranges along edge
						auto n_ranges = properties[neighbor].split_range_at(p);
						auto p_ranges = properties[p].split_range_at(neighbor);
					
						// Case 1: C[x,p) Cutvertex between x and p inclusive
						if(x == y)
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\tx == y\n";
							#endif
							
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
						
							properties[p].set_current_range(p_ranges.first);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, new_x,
								properties[new_x].get_current_range(), new_y,
								properties[new_y].get_current_range(), new_x, -1, before_y, before_p);
							
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							properties[neighbor].set_current_range(n_ranges.first);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
							
							// Color remaining "right" subgraph
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, x,
								properties[x].get_current_range(), y, properties[y].get_current_range(), p,
								-1, before_y, -1);
						}
						else if(compare_face_locations(neighbor_face_location, before_p, face_location_sets))
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\tbefore_p\n";
							#endif
							
							// Set ranges for "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							properties[neighbor].set_current_range(n_ranges.first);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
							
							// Color "right" subgraph
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, x,
								properties[x].get_current_range(), y, properties[y].get_current_range(), p,
								before_p, before_y, before_x);
							
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
						
							properties[p].set_current_range(p_ranges.first);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, neighbor,
								properties[neighbor].get_current_range(), neighbor,
								properties[neighbor].get_current_range(), neighbor, -1, before_p, -1);
						}
						// Case 2: C(p,y] Cutvertex between p and y inclusive
						else if(compare_face_locations(neighbor_face_location, before_y, face_location_sets))
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\tbefore_y\n";
							#endif
							
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
						
							properties[p].set_current_range(p_ranges.first);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Intersect face sections we are joining
							before_y = face_location_sets.take_union(before_p, before_y);
							before_p = -1;
							
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								std::cout << "union before_y and before_p = " << before_y << "\n";
							#endif
							
							// Color "left" subgraph with p removed
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, new_x, properties[new_x].get_current_range(),
								new_y, properties[new_y].get_current_range(), new_x, -1, before_y, before_x);
								
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							properties[neighbor].set_current_range(n_ranges.first);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
							
							// Color remaining "right" subgraph
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, neighbor,
								properties[neighbor].get_current_range(), neighbor,
								properties[neighbor].get_current_range(), neighbor, -1, before_y, -1);
						}
						//Case 3: C(y,x) Cutvertex between y and x exclusive
						else
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\tbefore_x\n";
							#endif
							
							// Reassign ranges for remaining "right" subgraph
							properties[p].set_current_range(p_ranges.second);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							properties[neighbor].set_current_range(n_ranges.first);
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
							
							// Color remaining "right" subgraph
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, neighbor,
								properties[neighbor].get_current_range(), y, properties[y].get_current_range(),
								p, -1, before_y, before_x);
							
							// Set ranges for "left" half to color
							properties[neighbor].set_current_range(n_ranges.second);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(neighbor, graph, embedding, properties);
							#endif
						
							properties[p].set_current_range(p_ranges.first);
						
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "\t\t\t\t";
								print_incidence_range(p, graph, embedding, properties);
							#endif
							
							// Remove p from incidence list
							properties[neighbor].remove_begin();
							properties[p].remove_end();
							
							// Color "left" subgraph with p removed
							path_list_color_recursive(graph, embedding, color_list, properties,
								face_location_sets, coloring, new_x, properties[new_x].get_current_range(),
								neighbor, properties[neighbor].get_current_range(), new_x, -1, before_p,
								before_x);
						}
					}
					
					break;
				}
			}
		}
	}
	
	template<
		typename index_graph, typename planar_embedding, typename list_color_property_map,
		typename edge_iterator = typename property_traits<planar_embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<index_graph>::vertex_descriptor,
		typename edge_iterator_pair = typename std::pair<edge_iterator, edge_iterator>
	>
	void print_incidence_range(vertex_descriptor vert, const index_graph & graph,
		const planar_embedding & embedding, const list_color_property_map & properties)
	{
		edge_iterator_pair current_range = properties[vert].get_current_range();
		vertex_descriptor begin_vert = get_incident_vertex(vert, *(current_range.first), graph);
		vertex_descriptor end_vert = get_incident_vertex(vert, *(--current_range.second), graph);
		
		edge_iterator_pair range = properties[vert].get_range();
		vertex_descriptor r_begin_vert = get_incident_vertex(vert, *(range.first), graph);
		vertex_descriptor r_end_vert = get_incident_vertex(vert, *(--range.second), graph);
		std::cout << "setting neighbor_range[" << vert << "] = (" << begin_vert << "," << end_vert
			<< ") (" << r_begin_vert << "," << r_end_vert << ")\n";
	}
	
	/*template<typename Graph, typename Embedding, typename IteratorMap,
		typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor>
	edge_iterator find_edge(const Graph & graph, const Embedding & embedding, IteratorMap & iterator_map,
		vertex_descriptor vertex, const std::pair<edge_iterator, edge_iterator> & range,
		vertex_descriptor goal_vertex)
	{
		// If we don't know the position, look it up and store it
		if(iterator_map.count(vertex) == 0 || iterator_map[vertex].count(goal_vertex) == 0)
		{
			auto ordering = embedding[vertex];
			for(auto edge_iter = range.first; edge_iter != range.second; ++edge_iter)
			{
				if(edge_iter == ordering.end()) edge_iter = ordering.begin();
				vertex_descriptor neighbor = get_incident_vertex(vertex, *edge_iter, graph);
			
				// Done if we found last vertex on path
				iterator_map[vertex][neighbor] = edge_iter;
			}
			
			if(iterator_map.count(vertex) == 0 || iterator_map[vertex].count(goal_vertex) == 0)
			{
				std::string error = "Find edge called for nonexistant neighbor " + std::to_string(goal_vertex);
				throw std::runtime_error(error);
			}
		}
		
		return iterator_map[vertex][goal_vertex];
	}*/
	
	/*template<typename Graph, typename Embedding, typename ColorList, typename Coloring, typename VertexIter>
	void path_list_color(const Graph & graph, const Embedding & embedding,  ColorList & color_list,
		Coloring & coloring, VertexIter begin, VertexIter end)
	{*/
		/*
		 * In this section we construct iterator property maps for the properties that must be tracked
		 * throghout the path 3-list-coloring algorithm. Property maps are used to allow essentially
		 * 'fast as possible' lookup alongside abstraction. For example if the graph has integer vertices
		 * then the vertex name may be used to directly lookup the property in a vector.
		 */
 /*
		// Define the storage type for integer property maps
		typedef iterator_property_map<
				std::vector<int>::iterator,
				typename property_map<Graph, vertex_index_t>::type
			> int_property_map;
		
		// Make state map
		std::vector<int> state_storage(num_vertices(graph));
		int_property_map state(state_storage.begin(), get(vertex_index, graph));
		
		// Make face location map
		std::vector<int> face_location_storage(num_vertices(graph));
		int_property_map face_location(face_location_storage.begin(), get(vertex_index, graph));
		
		// Define the storage type for edge iterator pair property map
		typedef typename property_traits<Embedding>::value_type::const_iterator edge_iterator;		
		
		typedef iterator_property_map<
				typename std::vector<std::pair<edge_iterator, edge_iterator> >::iterator,
				typename property_map<Graph, vertex_index_t>::type
			> edge_iterator_pair_map;
		
		// Make neighbor range map
		std::vector<std::pair<edge_iterator, edge_iterator> > neighbor_range_storage(num_vertices(graph));
		edge_iterator_pair_map neighbor_range(neighbor_range_storage.begin(), get(vertex_index, graph));

		// Define the storage type for edge iterator pair property map
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		std::unordered_map<
				vertex_descriptor,
				std::unordered_map<vertex_descriptor, edge_iterator>
			> iterator_map; */
		
		/*
		 * Here we set up the initial conditions for the graph based by iterating through the given
		 * outer face. All outer face vertices are set as between before_y, meaning between x and y with
		 * no colored path P in existance. We then set x and y to be the start and end of the path respectively.
		 * Neighbor ranges are initialized by finding neighboring face vertices in the incidence list. After
		 * these initializations, the recursive algorithm is applied to the graph to produce a path coloring.
		 */
/*	 
		 const int ON_FACE = 1;
		 const int BEFORE_P = 0, BEFORE_Y = 1, BEFORE_X = 2;
		 
		 for(auto vertex_iter = begin; vertex_iter != end; ++vertex_iter)
		 {
		 	auto previous = vertex_iter, next = vertex_iter;
		 	
		 	if(vertex_iter == begin)
		 	{
		 		previous = end;
		 		--previous;
		 	}
		 	else --previous;
		 	
		 	++next;
		 	if(next == end) next = begin;
		 	
		 	std::pair<edge_iterator, edge_iterator> iterator_range(embedding[*vertex_iter].begin(),
		 		embedding[*vertex_iter].end());
		 	
		 	edge_iterator first = find_edge(graph, embedding, iterator_map, *vertex_iter, iterator_range,
		 		*previous);
	 		edge_iterator second = find_edge(graph, embedding, iterator_map, *vertex_iter, iterator_range,
		 		*next);
	 		
	 		neighbor_range[*vertex_iter] = std::pair<edge_iterator, edge_iterator>(first, ++second);
	 		state[*vertex_iter] = ON_FACE;
	 		face_location[*vertex_iter] = BEFORE_Y;
		 }
		 
		 vertex_descriptor x = *begin, y = *(--end);
		 
		 path_list_color_recursive(graph, embedding, color_list, coloring, iterator_map, state, face_location,
		 	neighbor_range, x, neighbor_range[x], y, neighbor_range[y], x, BEFORE_P, BEFORE_Y, BEFORE_X);
	}
	
	template<typename Graph, typename Embedding, typename ColorList, typename Coloring, typename IteratorMap,
		typename StateMap, typename FaceLocationMap, typename IncidenceRangeMap,
		typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor,
		typename edge_iterator_pair = typename std::pair<edge_iterator, edge_iterator> >
	void path_list_color_recursive(const Graph & graph, const Embedding & embedding,  ColorList & color_list,
		Coloring & coloring, IteratorMap & iterator_map, StateMap & state, FaceLocationMap & face_location,
		IncidenceRangeMap & neighbor_range, vertex_descriptor x, edge_iterator_pair x_range,
		vertex_descriptor y, edge_iterator_pair y_range, vertex_descriptor p, int before_p, int before_y,
		int before_x)
	{
		//draw_list_color_state(graph, state);
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "x = " << x << ", y = " << y << ", p = " << p << "\n";
		#endif
		const int INTERIOR = 0, ON_FACE = 1, COLORED = 2;
		
		// Update neighbor_range for this subgraph
		neighbor_range[x] = x_range;
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\t";
			print_incidence_range(x, graph, embedding, neighbor_range);
		#endif
		neighbor_range[y] = y_range;
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\t";
			print_incidence_range(y, graph, embedding, neighbor_range);
		#endif
		
		// If colored path doesn't exist, we create one
		if(state[p] != COLORED)
		{
			auto path_color = color_list[x].front();
			vertex_descriptor path_end = x, next_vertex = x;
			coloring[x] = path_color;
			state[x] = COLORED;
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tcoloring path with path_color = " << path_color << "\n";
			#endif
			
			do
			{
				// Update the path end
				path_end = next_vertex;
				
				// Look counterclockwise through our current range of interior neighbors
				edge_iterator begin = neighbor_range[path_end].first, end = neighbor_range[path_end].second;
				auto ordering = embedding[path_end];
				for(auto edge_iter = begin; edge_iter != end; ++edge_iter)
				{
					if(edge_iter == ordering.end()) edge_iter = ordering.begin();
					vertex_descriptor neighbor = get_incident_vertex(path_end, *edge_iter, graph);
					
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\tchecking state of vertex " << neighbor << "\n";
					#endif
					
					// Check if vertex is on the outer face between x and y
					if(state[neighbor] == ON_FACE && face_location[neighbor] == before_y)
					{
						// Check if it may be colored path_color
						for(auto color : color_list[neighbor])
						{
							if(color == path_color)
							{
								// Color and append the vertex to the path
								next_vertex = neighbor;
								coloring[next_vertex] = path_color;
								state[next_vertex] = COLORED;
								
								#ifdef PLC_SHOW_ANNOTATIONS
									std::cout << "\t\t\tcolored " << next_vertex << "\n";
								#endif
								
								// Check if we took a chord
								edge_iterator prev_edge_iter = neighbor_range[next_vertex].first;
								vertex_descriptor prev_on_face =
									get_incident_vertex(next_vertex, *prev_edge_iter, graph);
								if(state[prev_on_face] != COLORED)
								{
									#ifdef PLC_SHOW_ANNOTATIONS
										std::cout << "\t\t\tlobe made by edge " << path_end
											<< "-" << next_vertex << "\n";
									#endif
									// Find chord edge in incidence list
									edge_iterator n_mid = find_edge(graph, embedding, iterator_map, next_vertex,
										neighbor_range[next_vertex], p);
									
									// Setup ranges for lobe x and y
									edge_iterator_pair lobe_x_range(neighbor_range[next_vertex].first, ++n_mid),
										lobe_y_range(edge_iter, neighbor_range[path_end].second);
									
									// Set end of range iterator one past chord edge
									edge_iterator new_pe_end = edge_iter;
									++new_pe_end;
									
									// Setup ranges for x and y in remaining subgraph
									edge_iterator_pair new_n_range(--n_mid, neighbor_range[next_vertex].second),
										new_pe_range(neighbor_range[path_end].first, new_pe_end);
									
									// Color lobe
									path_list_color_recursive(graph, embedding, color_list, coloring,
										iterator_map, state, face_location, neighbor_range, next_vertex,
										lobe_x_range, path_end, lobe_y_range, next_vertex,
										(before_y + 1) % 3, (before_y + 2) % 3, before_y);
									
									// Reassign ranges
									neighbor_range[next_vertex] = new_n_range;
									#ifdef PLC_SHOW_ANNOTATIONS
										std::cout << "\t\t\t";
										print_incidence_range(next_vertex, graph, embedding, neighbor_range);
									#endif
									neighbor_range[path_end] = new_pe_range;
									#ifdef PLC_SHOW_ANNOTATIONS
										std::cout << "\t\t\t";
										print_incidence_range(path_end, graph, embedding, neighbor_range);
									#endif
								}
								
								break;
							}
						}
					}
					
					// If we have found a new end vertex, stop looking at neighbors
					if(path_end != next_vertex) break;
				}
			}
			while(path_end != next_vertex);
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tdone coloring path\n";
			#endif
		}
		
		// Base Case 1: K_1
		if(x_range.first == x_range.second)
		{
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tbase case: K_1\n";
			#endif
			
			// Color x if not colored
			if(state[x] != COLORED)
			{
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\tcolored " << x << "\n";
				#endif
				state[x] = COLORED;
				coloring[x] = color_list[x].front();
			}
			
			return;
		}
		
		// Base Case 2: K_2
		if(x_range.first == --x_range.second)
		{
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\tbase case: K_2\n";
			#endif
			
			// Color x if not colored
			if(state[x] != COLORED)
			{
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\tcolored " << x << "\n";
				#endif
				state[x] = COLORED;
				coloring[x] = color_list[x].front();
			}
			
			// Color other vertex if not colored
			vertex_descriptor neighbor = get_incident_vertex(x, *x_range.first, graph);
			if(state[neighbor] != COLORED)
			{
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\tcolored " << neighbor << "\n";
				#endif
				state[neighbor] = COLORED;
				coloring[neighbor] = color_list[neighbor].front();
			}
			
			return;
		}
		
		// New parameters for coloring the graph after removing p
		vertex_descriptor new_x = x, new_y = y, new_p = x;
		
		// Grab important positions in p's incidence list
		edge_iterator begin = neighbor_range[p].first, end = neighbor_range[p].second;
		edge_iterator pre_end = end;
		--pre_end; 
		
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\tremoving p = " << p << "\n";
		#endif
		
		// Iterate counterclockwise through interior neighbors of p
		auto ordering = embedding[p];
		for(auto edge_iter = begin; edge_iter != end; ++edge_iter)
		{
			if(edge_iter == ordering.end()) edge_iter = ordering.begin();
			vertex_descriptor neighbor = get_incident_vertex(p, *edge_iter, graph);
			
			#ifdef PLC_SHOW_ANNOTATIONS
				std::cout << "\t\tlooking at neighbor = " << neighbor << "\n";
			#endif
			
			// Remove p's color from the list of all adjacent neighbors
			std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), coloring[p]);
			
			// Interior vertex
			if(state[neighbor] == INTERIOR)
			{
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\t\tinterior\n";
				#endif
				
				// Mark neighbor as on face before p since we are removing p
				state[neighbor] = ON_FACE;
				face_location[neighbor] = before_p;
				
				// Find edge to p in adjacency list
				auto n_ordering = embedding[neighbor];
				edge_iterator n_begin = find_edge(graph, embedding, iterator_map, neighbor,
					edge_iterator_pair(n_ordering.begin(), n_ordering.end()), p);
				edge_iterator n_end = n_begin;
				
				// Assign adjacency list to new face vertex
				if(++n_begin == n_ordering.end()) n_begin = n_ordering.begin();
				if(n_end == n_ordering.begin()) n_end = n_ordering.end();
				neighbor_range[neighbor] = edge_iterator_pair(n_begin, n_end);
				#ifdef PLC_SHOW_ANNOTATIONS
					std::cout << "\t\t\t\t";
					print_incidence_range(neighbor, graph, embedding, neighbor_range);
				#endif
			}
			else
			{
				// Vertex prior to p on the outer face
				if(edge_iter == begin)
				{
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\tbegin\n";
					#endif
					
					// Reassign x or y as needed
					if(x == p)
					{
						if(y == p) new_y = neighbor;
						else new_x = neighbor;
					}
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\t\t TEST ";
						print_incidence_range(neighbor, graph, embedding, neighbor_range);
					#endif
					// Reassign neighbor range
					if(--neighbor_range[neighbor].second == embedding[neighbor].begin())
					{
						neighbor_range[neighbor].second = embedding[neighbor].end();
					}
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\t\t";
						print_incidence_range(neighbor, graph, embedding, neighbor_range);
					#endif
				}
				else
				{
					// Reassign x or y as needed
					if(y == p)
					{
						if(x == p) new_x = neighbor;
						else new_y = neighbor;
					}
					
					if(edge_iter == pre_end)
					{
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\tend\n";
						#endif
						
						// Set new_p to the next vertex along the face
						new_p = neighbor;
					
						// Reassign neighbor range
						++neighbor_range[neighbor].first;
						if(neighbor_range[neighbor].first == embedding[neighbor].end())
						{
							neighbor_range[neighbor].first = embedding[neighbor].begin();
						}
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\t\t";
							print_incidence_range(neighbor, graph, embedding, neighbor_range);
						#endif
					}
					else
					{
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\tcutvertex ";
						#endif
						
						// Save old range
						edge_iterator_pair old_neighbor_range = neighbor_range[neighbor];
						
						// Find edge to p in adjacency list
						edge_iterator p_edge_iter =
							find_edge(graph, embedding, iterator_map, neighbor, neighbor_range[neighbor], p);
					
						// Find midpoint iterator
						edge_iterator new_y_begin = p_edge_iter;
						if(++new_y_begin == embedding[neighbor].end())
							new_y_begin == embedding[neighbor].begin();
						
						// Find range for new y in the block
						neighbor_range[neighbor] =
							edge_iterator_pair(new_y_begin, neighbor_range[neighbor].second);
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\t\t";
							print_incidence_range(neighbor, graph, embedding, neighbor_range);
						#endif
					
						// Case 1: C[x,p) Cutvertex between x and p inclusive
						if(face_location[neighbor] == before_p)
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "before_p\n";
							#endif
							
							path_list_color_recursive(graph, embedding, color_list, coloring, iterator_map,
								state, face_location, neighbor_range, neighbor, neighbor_range[neighbor],
								neighbor, neighbor_range[neighbor], neighbor, (before_p + 1) % 3,
								(before_p + 2) % 3, before_p);
						}
						// Case 2: C(p,y] Cutvertex between p and y inclusive
						else if(face_location[neighbor] == before_y)
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "before_y\n";
							#endif
							
							
							path_list_color_recursive(graph, embedding, color_list, coloring, iterator_map,
								state, face_location, neighbor_range, new_x, neighbor_range[new_x], new_y,
								neighbor_range[new_y], new_x, before_p, before_y, before_x);
								
							// Reassign y since it is no longer in the subgraph
							new_y = neighbor;
							// Reassign new_x as we have colored everything to the "left"
							new_x = neighbor;
						}
						//Case 3: C(y,x) Cutvertex between y and x exclusive
						else
						{
							#ifdef PLC_SHOW_ANNOTATIONS
								std::cout << "before_x\n";
							#endif
						
							path_list_color_recursive(graph, embedding, color_list, coloring, iterator_map,
								state, face_location, neighbor_range, new_x, neighbor_range[new_x], new_y,
								neighbor_range[new_y], new_x, before_p, before_y, before_x);
						
							// Reassign new_x as we have colored everything to the "left"
							new_x = neighbor;
						}
					
						// Reassign neighbor ranges for remaining graph to color
						neighbor_range[neighbor] = old_neighbor_range;
						neighbor_range[neighbor].second = ++p_edge_iter;
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\t";
							print_incidence_range(neighbor, graph, embedding, neighbor_range);
						#endif
						neighbor_range[p].first = edge_iter;
						#ifdef PLC_SHOW_ANNOTATIONS
							std::cout << "\t\t\t";
							print_incidence_range(p, graph, embedding, neighbor_range);
						#endif
					
						break;
					}
				}
			}
		}
		
		// Color what is left
		path_list_color_recursive(graph, embedding, color_list, coloring, iterator_map, state, face_location,
			neighbor_range, new_x, neighbor_range[new_x], new_y,  neighbor_range[new_y],
			new_p, before_p, before_y, before_x);
	}
	
	// A class to hold the coordinates of the straight line embedding
	struct coord_t
	{
	  std::size_t x;
	  std::size_t y;
	};
	
	template<typename Graph, typename State>
	void draw_list_color_state(const Graph & graph, const State & state)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
	
		// Define the storage type for the planar embedding
		typedef std::vector<
				std::vector<edge_descriptor>
			> embedding_storage_t;

		typedef iterator_property_map
			< typename embedding_storage_t::iterator, 
				typename property_map<Graph, vertex_index_t>::type
			> embedding_t;
		
		// Create the planar embedding
		embedding_storage_t embedding_storage(num_vertices(graph));
		embedding_t embedding(embedding_storage.begin(), get(vertex_index, graph));

		boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
		   boyer_myrvold_params::embedding = embedding);
		
		// Find a canonical ordering
		std::vector<vertex_descriptor> ordering;
		
		planar_canonical_ordering(graph, embedding, std::back_inserter(ordering));
		
		//Set up a property map to hold the mapping from vertices to coord_t's
		typedef std::vector< coord_t > straight_line_drawing_storage_t;
		typedef boost::iterator_property_map
			< straight_line_drawing_storage_t::iterator, 
				typename property_map<Graph, vertex_index_t>::type 
			> straight_line_drawing_t;

		straight_line_drawing_storage_t straight_line_drawing_storage
			(num_vertices(graph));
		straight_line_drawing_t straight_line_drawing
			(straight_line_drawing_storage.begin(), 
				get(vertex_index, graph)
			);



		// Compute the straight line drawing
		chrobak_payne_straight_line_drawing(graph, 
				embedding, 
				ordering.begin(),
				ordering.end(),
				straight_line_drawing
			);
	
		std::cout << draw_tikz_graph(graph, state, straight_line_drawing) << "\n";
	}
	
	template<typename Graph, typename Embedding, typename IncidenceRangeMap,
		typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor,
		typename edge_iterator_pair = typename std::pair<edge_iterator, edge_iterator> >
	void print_incidence_range(vertex_descriptor vert, const Graph & graph, const Embedding & embedding,
		const IncidenceRangeMap & neighbor_range)
	{
		edge_iterator begin = neighbor_range[vert].first;
		edge_iterator end = neighbor_range[vert].second;
		vertex_descriptor begin_vert = get_incident_vertex(vert, *begin, graph);
		vertex_descriptor end_vert = get_incident_vertex(vert, *(--end), graph);
		std::cout << "setting neighbor_range[" << vert << "] = (" << begin_vert << "," << end_vert << ")\n";
	}
*/
}

#endif