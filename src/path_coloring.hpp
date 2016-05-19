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

#define PLC_SHOW_ANNOTATIONS

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
			for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
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
			for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
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
			for(auto p_iter = p_begin; p_iter != p_end; p_iter++)
			{
				p_map[*p_iter] = p_iter;
			}
			
			// Check vertices in path q for neighbors in p
			for(auto q_iter = q_begin; q_iter != q_end; q_iter++)
			{
				// Look at each neighbor
				auto ordering = embedding[*q_iter];
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
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
			for(auto p_iter = p_begin; p_iter != p_end; p_iter++)
			{
				parent_map[*p_iter] = *p_iter;
			}
			
			// Mark vertices in path q
			for(auto q_iter = q_begin; q_iter != q_end; q_iter++)
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
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
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
	
	template<typename Graph, typename Embedding, typename IteratorMap,
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
			for(auto edge_iter = range.first; edge_iter != range.second; edge_iter++)
			{
				if(edge_iter == ordering.end()) edge_iter = ordering.begin();
				vertex_descriptor neighbor = get_incident_vertex(vertex, *edge_iter, graph);
			
				// Done if we found last vertex on path
				iterator_map[vertex][neighbor] = edge_iter;
			}
			
			if(iterator_map.count(vertex) == 0)
			{
				std::string error = "Find edge called for nonexistant neighbor " + std::to_string(goal_vertex);
				throw std::runtime_error(error);
			}
		}
		
		return iterator_map[vertex][goal_vertex];
	}
	
	template<typename Graph, typename Embedding, typename ColorList, typename Coloring, typename VertexIter>
	void path_list_color(const Graph & graph, const Embedding & embedding,  ColorList & color_list,
		Coloring & coloring, VertexIter begin, VertexIter end)
	{
		/*
		 * In this section we construct iterator property maps for the properties that must be tracked
		 * throghout the path 3-list-coloring algorithm. Property maps are used to allow essentially
		 * 'fast as possible' lookup alongside abstraction. For example if the graph has integer vertices
		 * then the vertex name may be used to directly lookup the property in a vector.
		 */
		 
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
			> iterator_map;
		
		/*
		 * Here we set up the initial conditions for the graph based by iterating through the given
		 * outer face. All outer face vertices are set as between before_y, meaning between x and y with
		 * no colored path P in existance. We then set x and y to be the start and end of the path respectively.
		 * Neighbor ranges are initialized by finding neighboring face vertices in the incidence list. After
		 * these initializations, the recursive algorithm is applied to the graph to produce a path coloring.
		 */
		 
		 const int ON_FACE = 1;
		 const int BEFORE_P = 0, BEFORE_Y = 1, BEFORE_X = 2;
		 
		 for(auto vertex_iter = begin; vertex_iter != end; vertex_iter++)
		 {
		 	auto previous = vertex_iter, next = vertex_iter;
		 	
		 	if(vertex_iter == begin)
		 	{
		 		previous = end;
		 		previous--;
		 	}
		 	else previous--;
		 	
		 	next++;
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
		neighbor_range[y] = y_range;
		
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
				for(auto edge_iter = begin; edge_iter != end; edge_iter++)
				{
					if(edge_iter == ordering.end()) edge_iter = ordering.begin();
					vertex_descriptor neighbor = get_incident_vertex(path_end, *edge_iter, graph);
					
					#ifdef PLC_SHOW_ANNOTATIONS
						std::cout << "\t\t\tchecking state of vertex " << neighbor << "\n";
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
									std::cout << "\t\tcolored " << next_vertex << "\n";
								#endif
								
								// Check if we took a chord
								edge_iterator prev_edge_iter = neighbor_range[next_vertex].first;
								vertex_descriptor prev_on_face =
									get_incident_vertex(next_vertex, *prev_edge_iter, graph);
								if(state[prev_on_face] != COLORED)
								{
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
									neighbor_range[path_end] = new_pe_range;
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
				state[x] = COLORED;
				coloring[x] = color_list[x].front();
			}
			
			// Color other vertex if not colored
			vertex_descriptor neighbor = get_incident_vertex(x, *x_range.first, graph);
			if(state[neighbor] != COLORED)
			{
				state[neighbor] = COLORED;
				coloring[neighbor] = color_list[neighbor].front();
			}
			
			return;
		}
		
		// New parameters for coloring the graph after removing p
		vertex_descriptor new_x = x, new_y = y, new_p = x;
		
		// Grab important positions in p's incidence list
		edge_iterator begin = neighbor_range[p].first, end = neighbor_range[p].second, pre_end = end;
		pre_end--; 
		
		#ifdef PLC_SHOW_ANNOTATIONS
			std::cout << "\tremoving p = " << p << "\n";
		#endif
		
		// Iterate counterclockwise through interior neighbors of p
		auto ordering = embedding[p];
		for(auto edge_iter = begin; edge_iter != end; edge_iter++)
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
				edge_iterator n_begin = find_edge(graph, embedding, iterator_map, neighbor,
					edge_iterator_pair(embedding[neighbor].begin(), embedding[neighbor].end()), p);
				edge_iterator n_end = n_begin;
				
				// Assign adjacency list to new face vertex
				auto n_ordering = embedding[neighbor];
				if(++n_begin == ordering.end()) n_begin = ordering.begin();
				if(n_end == ordering.begin()) n_end = ordering.end();
				neighbor_range[neighbor] = edge_iterator_pair(n_begin, n_end);
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
					
					// Reassign neighbor range
					neighbor_range[neighbor].second--;
					if(neighbor_range[neighbor].second == embedding[neighbor].begin())
					{
						neighbor_range[neighbor].second = embedding[neighbor].end();
					}
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
						neighbor_range[neighbor].first++;
						if(neighbor_range[neighbor].first == embedding[neighbor].end())
						{
							neighbor_range[neighbor].first = embedding[neighbor].begin();
						}
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
						if(++new_y_begin == embedding[neighbor].end()) new_y_begin == embedding[neighbor].begin();
						
						// Find range for new y in the block
						neighbor_range[neighbor] =
							edge_iterator_pair(new_y_begin, neighbor_range[neighbor].second);
					
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
						neighbor_range[p].first = edge_iter;
					
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
}

#endif