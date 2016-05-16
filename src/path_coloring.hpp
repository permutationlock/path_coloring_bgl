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
	
	template<typename Graph, typename Embedding, typename IteratorMap
		typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor>
	edge_iterator find_edge(const Graph & graph, const Embedding & embedding, IteratorMap & iterator_map,
		vertex_descriptor vertex, const std::pair<edge_iterator, edge_iterator> & range,
		vertex_descriptor goal_vertex)
	{
		// Dynamically save all lookups for future lookups
		static std::unordered_map<vertex_descriptor,
			std::unordered_map<vertex_descriptor, edge_iterator> > iterator_map;
		
		if(iterator_map.count(vertex) && iterator_map[vertex].count(goal_vertex))
		{
			auto ordering = embedding[vertex];
			for(auto edge_iter = range.first; edge_iter != range.end; edge_iter++)
			{
				if(edge_iter == ordering.end()) edge_iter = ordering.begin();
				vertex_descriptor neighbor = get_incident_vertex(vertex, *edge_iter, graph);
			
				// Done if we found last vertex on path
				iterator_map[vertex][neighbor] = edge_iter;
			}
		}
		
		return iterator_map[vertex][goal_vertex];
	}
	
	template<typename Graph, typename Embedding, typename ColorList, typename Coloring, typename IteratorMap
			typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
			typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor,
			typename color_type = typename property_traits<Coloring>::value_type,
			typename state_map = typename std::unordered_map<vertex_descriptor, bool>,
			typename face_location_map = typename std::unordered_map<vertex_descriptor, bool>,
			typename edge_iterator_pair = typename std::pair<edge_iterator, edge_iterator>,
			typename incidence_range_map = typename std::unordered_map<vertex_descriptor, edge_iterator_pair>
		>
	void path_list_color_recursive(const Graph & graph, const Embedding & embedding, ColorList & color_list,
		Coloring & coloring, state_map & state, face_location_map & face_location,
		incidence_range_map & neighbor_range, vertex_descriptor x, edge_iterator_pair x_range,
		vertex_descriptor y, edge_iterator_pair y_range, vertex_descriptor p, int before_p, int before_y,
		int before_x)
	{
		// Update neighbor_range for this subgraph
		neighbor_range[x] = x_range;
		neighbor_range[y] = y_range;
		
		// Base Case 1: K_1
		if(x_begin == x_end)
		{
			// Color x if not colored
			if(!state[x]) color[x] = color_list[x].first();
			
			return;
		}
		
		// Base Case 2: K_2
		if(++x_begin == x_end)
		{
			// Color x and/or y if not colored
			if(!state[x]) color[x] = color_list[x].first();
			if(!state[y]) color[y] = color_list[y].first();
			
			return;
		}
		
		// If colored path doesn't exist, we create one
		if(!state[p])
		{
			vertex_descriptor path_end = x, next_vertex = x;
			
			do
			{
				// Update the path end
				path_end = next_vertex;
				
				// Color and mark the new path vertex
				coloring[path_end] = path_color;
				state[path_end] = true;
				
				// Look counterclockwise through our current range of interior neighbors
				edge_iterator begin = neighbor_range[path_end].first, end = neighbor_range[path_end].second;
				auto ordering = embedding[path_end];
				for(auto edge_iter = begin; edge_iter != end; edge_iter++)
				{
					if(edge_iter == ordering.end()) edge_iter = ordering.begin();
					vertex_descriptor neighbor = get_incident_vertex(path_end, *edge_iter, graph);
					
					// Check if vertex is on the outer face between x and y
					if(state.count(neighbor) == 1 && face_location[neighbor] == before_y)
					{
						// Should never find a colored vertex while making path
						if(state[neighbor]) throw std::runtime_error("Colored vertex found on face.");
						
						// Check if it may be colored path_color
						for(auto color : color_list[neighbor])
						{
							if(color == path_color)
							{
								// Assign the next vertex in the path
								next_vertex = neighbor;
								
								// Check if we took a chord
								prev_edge_iter = neighbor_range[next_vertex].first;
								vertex_descriptor prev_on_cycle =
									get_incident_vertex(next_vertex, *(--prev_edge_iter), graph);
								if(state[prev_on_cycle] != COLORED)
								{
									// Find cut point in incidence list
									edge_iterator n_mid = find_edge(graph, embedding, next_vertex,
										neighbor_range[next_vertex], p);
									
									// Setup ranges for lobe x and y
									edge_iterator_pair new_x_range(n_begin, ++n_mid)
										new_y_range(edge_iter, neibhorhood[path_end].second);
									
									// Color lobe
									path_list_color_recursive(graph, embedding, color_list, coloring, state,
										face_location, neighbor_range, next_vertex, new_x_range,
										path_end, new_y_range, next_vertex,
										(before_y + 1) % 3, (before_y + 2) % 3, before_y);
									
									// Reassign range for end of chord
									neighborhood[next_vertex].first = --n_mid;
									
									// Correction to make sure list wraps correctly
									if(edge_iter == ordering.begin()) edge_iter = ordering.end();
									
									// Reassign range for  beginning of chord
									neighborhood[path_end].second = edge_iter;
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
		}
		
		
		// New parameters for coloring the graph after removing p
		vertex_descriptor new_x = x, new_y = y, new_p = p;
		
		// Grab important positions in p's incidence list
		edge_iterator begin = neighbor_range[p].first, end = neighbor_range[p].second, pre_end = end;
		pre_end--; 
		
		// Iterate counterclockwise through interior neighbors of p
		auto ordering = embedding[p];
		for(auto edge_iter = begin; edge_iter != end; edge_iter++)
		{
			if(edge_iter == ordering.end()) edge_iter = ordering.begin();
			
			vertex_descriptor neighbor = get_incident_vertex(p, *edge_iter, graph);
			
			// Remove p's color from the list of all adjacent neighbors
			std::remove(color_list[neighbor].begin(), color_list[neighbor].end(), color[p]);
			
			// Interior vertex
			if(state.count(neighbor) == 0)
			{
				// Mark neighbor as on face before p since we are removing p
				state[neighbor] = false;
				face_location[neighbor] = before_p;
				
				// Find edge to p in adjacency list
				edge_iterator n_begin = find_edge(graph, embedding, neighbor, neighbor_range[neighbor], p),
					n_end = n_begin;
				
				// Assign adjacency list to new face vertex
				auto n_ordering = embedding[neighbor];
				if(++n_begin == ordering.end()) n_begin = ordering.begin();
				if(n_end = ordering.begin()) n_end = ordering.end();
				neighbor_range[neighbor] = edge_iterator_pair(n_begin, --n_end);
			}
			else
			{
				// Vertex prior to p on the outer face
				if(edge_iter == begin)
				{
					// Reassign x if it is being removed
					if(new_x == p)
					{
						new_x = neighbor;
					}
					
					// Reassign neighbor range
					if(neighbor_range[neighbor].second != embedding[neighbor].end())
					{
						neighbor_range[neighbor].second--;
					}
					else
					{
						neighbor_range[neighbor].second = embedding[neighbor].begin()
					}
				}
				// Vertex following p on the outer face
				else if(edge_iter == pre_end)
				{
					// Reassign y if it is being removed
					if(y == p)
					{
						new_y = neighbor;
					}
					
					// Set new_p to the next vertex along the face
					new_p = neighbor;
					
					// Reassign neighbor range
					neighbor_range[neighbor].first++;
					if(neighbor_range[neighbor].first == embedding[neighbor].end())
					{
						neighbor_range[neighbor].first = embedding[neighbor].begin();
					}
				}
				// Cutvertex
				else
				{
					// Find edge to p in adjacency list
					edge_iterator p_edge_iter =
						find_edge(graph, embedding, neighbor, neighbor_range[neighbor], p);
					
					// Find midpoint iterator
					edge_iterator new_y_begin = p_edge_iter;
					if(++new_y_begin == embedding[neighbor].end()) new_y_begin == embedding[neighbor].begin();
					
					// Find range for new y in the block
					edge_iterator_pair new_y_range(new_y_begin, neighbor_range[neighbor].second);
					
					// Case 1: C[x,p) Cutvertex between x and p inclusive
					if(face_location[neighbor] == before_p)
					{
						path_list_color_recursive(graph, embedding, color_list, coloring, state, face_location,
							neighbor_range, neighbor, new_y_range, neighbor, new_y_range, neighbor,
							(before_p + 1) % 3, (before_p + 2) % 3, before_p);
					}
					// Case 2: C(p,y] Cutvertex between p and y inclusive
					else if(face_location[neighbor] == before_y)
					{
						path_list_color_recursive(graph, embedding, color_list, coloring, state, face_location,
							neighbor_range, new_x, neighbor_range[new_x], new_y, new_y_range, new_x,
							before_p, before_y, before_x);
					}
					//Case 3: C(y,x) Cutvertex between y and x exclusive
					else
					{
						// Must explicitly assing range to cutvertex since it won't be x or y in this call
						neighbor_range[neighbor] = new_y_range;
						
						path_list_color_recursive(graph, embedding, color_list, coloring, state, face_location,
							neighbor_range, new_x, neighbor_range[new_x], new_y,  neighbor_range[new_y],
							new_x, before_p, before_y, before_x);
						
						// Reassign y since it is no longer in the subgraph
						new_y = neighbor;
					}
					
					// Reassign neighbor ranges for remaining graph to color
					neighbor_range[neighbor].second = ++p_edge_iter;
					neighbor_range[p].first = edge_iter;
					
					// Reassign new_x as we have colored everything to the "left"
					new_x = neighbor;
				}
			}
		}
		
		// Color what is left
		path_list_color_recursive(graph, embedding, color_list, coloring, state, face_location,
			neighbor_range, new_x, neighbor_range[new_x], new_y,  neighbor_range[new_y],
			new_p, before_p, before_y, before_x);
	}
}

#endif