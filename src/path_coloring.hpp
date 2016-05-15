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
	/*template<typename Graph, typename Embedding, typename ColorList, typename Coloring,
		typename VertexList>
	void hartman_path_list_color_block(const Graph & graph, const Embedding & embedding,
		VertexList & top, VertexList & bottom,
		const ColorList & color_list, Coloring & coloring)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		typedef typename property_traits<Coloring>::value_type color_type;
		
		// Select the first color in x's list
		color_type path_color = *(color_list[*p_begin].begin());
		
		// Find and color path of outer face vertices between x and y
		std::list<VertexIter> path;
		{
			std::unordered_map<vertex_descriptor, VertexIter> p_map;
		
			for(auto p_iter = p_begin; p_iter != p_end; p_iter++)
			{
				p_map[*p_iter] = p_iter;
			}
			
			path.push_back(p_begin);
			coloring[*p_begin] = path_color;
			
			bool done = false;
			while(!done)
			{
				done = true;
				
				// Iterate through incidentedges of the end of the path
				auto ordering = embedding[*(path.back())];
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
				{
					// Find neighbor vertex
					vertex_descriptor neighbor = get_incident_vertex(*(path.back()), *edge_iter, graph);
				
					// If neighbor is on the outer face between p_0 and q_0
					if(p_map.count(neighbor) != 0)
					{
						path.push_back(p_map.at(neighbor));
						coloring[neighbor] = path_color;
						done = false;
					}
				}
			}
		}
		
		// Determine lobes and center
		{
			for(std::size_t i = 0; i < path.size(); i++)
			{
				// Check to see if we take a chord
				auto p_iter = path[i];
				if(path[i+1] != p_iter++)
				{
					// Deal with the lobe that was chopped off by the chord
					std::list<vertex_descriptor> new_path = { *path[i], *path[i+1] };
					hartman_remove_path_vertex(graph, embedding, p_iter, path[i+1], p_iter, p_iter,
						new_path.begin(), new_path.end(), color_list, coloring, path_color);
				}
			}
		}
		
		// Remove path and deal with the interior
		{
			// Find path to remove
			std::list<vertex_descriptor> new_path;
			for(VertexIter p_iter : path)
			{
				new_path.push_back(*p_iter);
			}
			hartman_remove_path_vertex(graph, embedding, path.back(), p_end(), q_begin(), q_end(),
				new_path.begin(), new_path.end(), color_list, coloring, path_color);
		}
	}*/
	
	/*template<typename Graph, typename Embedding, typename VertexList, typename VertexMap,
		typename ColorList, typename Coloring>
	void path_list_color(const Graph & graph, const Embedding & embedding,
		VertexList & outer_face, VertexList & colored_path, VertexList::iterator p,
		VertexList::iterator x, VertexList::iterator y, VertexMap & marked
		ColorList & color_list, Coloring & coloring)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		typedef typename property_traits<Coloring>::value_type color_type;
		
		// If path is null, create path
		if(colored_path.size() == 0)
		{
			// Retrieve the first color in L(x)
			color_type path_color = *(color_list[*x].begin());
			
			bool done = false;
			while(!done)
			{
				done = true;
				
				// Iterate through incidentedges of the end of the path
				auto ordering = embedding[*(path.back())];
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
				{
					// Find neighbor vertex
					vertex_descriptor neighbor = get_incident_vertex(*(path.back()), *edge_iter, graph);
				
					// If neighbor is on the outer face between p_0 and q_0
					if(p_map.count(neighbor) != 0)
					{
						path.push_back(p_map.at(neighbor));
						coloring[neighbor] = path_color;
						done = false;
					}
				}
			}
		}
		
		//
	}*/
	
	/*
	template<typename VertexDescriptor, typename EdgeIter>
	struct SubgraphTraits
	{
		VertexDescriptor x, p, y;
		std::size_t x_begin, x_end
		EdgeIter y_begin, y_end;
		unsigned int bp, by, bx;
		
		SubgraphTraits(VertexDescriptor x, EdgeIter x_begin, EdgeIter x_end, VertexDescriptor y,
			EdgeIter y_begin, EdgeIter y_ends, VertexDescriptor p, unsigned int bp, unsigned int by,
			unsigned int bx) : x(x), p(p), y(y), x_begin(x_begin), x_end(x_end), y_begin(y_begin),
			y_end(y_end), bp(bp), by(by), bx(bx) {}
	};
	
	template<typename Graph, typename Embedding, typename VertexList, typename ColorList, typename Coloring>
	void path_list_color(const Graph & graph, const Embedding & embedding,
		VertexList & outer_face, ColorList & color_list, Coloring & coloring)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
		typedef typename property_traits<Coloring>::value_type color_type;
		typedef typename boost::property_traits<Embedding>::value_type::const_iterator edge_iterator;
		
		enum State { interior = 0, on_face = 1, colored = 2 };
		
		// Tracks the range in which uncolored neighbors lie for the given vertex
		std::unordered_map<vertex_descriptor, std::pair<edge_iterator, edge_iterator> > neighbor_range;
		
		// Tracks the current state of each vertex in the graph
		std::unordered_map<vertex_descriptor, State> state;
		std::unordered_map<vertex_descriptor, unsigned int> face_location;
		
		// Setup initial conditions for the outer face of the graph
		for(auto face_iter = outer_face.begin(); face_iter != outer_face.end(); face_iter++)
		{
			vertex_iter vertex = *face_iter;
			
			// Mark vertex as between x and y
			state[*face_iter] = Sate::on_face;
			face_location[*face_iter] = 0;
			
			// Computer next and previous vertices on outer face
			vertex_iter previous, next;
			if(face_iter + 1 == outer_face.end())
			{
				next = outer_face.front();
			}
			else
			{
				next = *(face_iter + 1);
			}
			if(face_iter == outer_face.begin())
			{
				previous = outer_face.back();
			}
			else
			{
				previous = *(face_iter - 1);
			}
			
			// Note the index of previous and next vertices on the outer face in incidence list
			edge_iterator begin, end;
			auto ordering = embedding[vertex];
			for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); edge_iter++)
			{
				vertex_descriptor neighbor = get_incident_vertex(*q_iter, *edge_iter, graph);
				if(neighbor == previous)
				{
					begin = edge_iter; 
				}
				else if(neighbor == next)
				{
					end = edge_iter;
					end++;
				}
			}
			neighbor_range[vertex] = std::pair<edge_iterator, edge_iterator>(begin, end);
		}
		
		// Tracks all remaining graphs to be colored
		std::stack<SubgraphTraits<vertex_descriptor> > subgraphs;
		
		// Start with no path, outer face of graph, and set x and y to the start and end of the face
		subgraphs.push(SubgraphTraits(outer_face.front(), neighbor_range[outer_face.front()].first,
			neighbor_range[outer_face.front()].second, outer_face.back(),
			neighbor_range[outer_face.back()].first, neighbor_range[outer_face.back()].second,
			outer_face.front(), 2, 0, 1));
		
		// Loop until all subgraphs have been colored
		while(!subgraphs.empty())
		{
			// Get current graph info
			SubgraphTraits subgraph = subgraphs.front();
			subgraphs.pop();
			
			// Update neighbor_range for this iteration
			neighbor_range[subgraph.x] =
				std::pair<edge_iterator, edge_iterator>(subgraph.x_begin, subgraph.x_end);
			neighbor_range[subgraph.y] =
				std::pair<edge_iterator, edge_iterator>(subgraph.y_begin, subgraph.y_end);
			
			// If a colored path doesn't exist, we must create one
			if(state.at(subgraph.p) != State::colored)
			{
				// Start path at x and use the first (and possibly only) color in x's list
				subgraph.p = x;
				vertex_descriptor path_end = subgraph.x, next_vertex = subgraph.x;
				color_type path_color = color_list[subgraph.x].first();
				
				// Color an induced path of vertices along the outer face between x and y with path color
				do
				{
					// Update the path end
					path_end = next_vertex;
					
					// Color and mark next vertex
					coloring[next_vertex] = path_color;
					state[next_vertex] = State::colored;
					
					// Look counterclockwise through our current range of interior neighbors
					edge_iterator begin = neighbor_range[path_end].first, end = neighbor_range[path_end].second;
					auto ordering = embedding[path_end];
					for(auto edge_iter = begin; edge_iter != end; edge_iter++)
					{
						if(edge_iter == ordering.end()) edge_iter = ordering.begin();
						vertex_descriptor neighbor = get_incident_vertex(path_end, *edge_iter, graph);
						
						// Check if we have hit previous vertex in the path
						if(state.at(neighbor) == State::colored)
						{
							break;
						}
						// Check if vertex is on the outer face between x and y
						else if(state.at(neighbor) == State::on_face &&
							face_location.at(neighbor) == FaceLoc::before_y)
						{
							// Check if it may be colored path_color
							for(auto color : color_list[neighbor])
							{
								if(color == path_color)
								{
									// Assign the next vertex in the path
									next_vertex = neighbor;
									
									// Check if we took a chord
									prev_edge_iter = neighbor_range.at(next_vertex).end;
									vertex_descriptor prev_on_cycle =
										get_incident_vertex(next_vertex, *(--prev_edge_iter), graph);
									if(state.at(prev_on_cycle) != State::colored)
									{
										// Find cut point in incidence list
										edge_iter n_begin = neighbor_range.at(next_vertex).first,
											n_end = neighbor_range.at(next_vertex).second, n_mid;
										for(n_mid = ++n_begin; n_mid != prev_on_cycle; n_mid++)
										{
											if(n_mid == ordering.end()) n_mid = ordering.begin();
											vertex_descriptor n_neighbor =
												get_incident_vertex(next_vertex, *n_mid, graph);
											
											// Find chord edge
											if(prev_on_cycle == n_neighbor)	break;
										}
										
										// Reassign range for end of chord
										neighborhood[next_vertex].first = mid_edge;
										
										// Set up the lobe to be colored
										subgraphs.push(Subgraph(next_vertex, n_begin, n_mid, path_end,
											edge_iter, neighborhood.at(path_end).second, next_vertex,
											subgraph.by, (subgraph.by + 1) % 3, (subgraph.by + 2) % 3));
										
										// Reassign range for  beginning of chord
										neighborhood[path_end].second = edge_iter;
									}
									
									break;
								}
							}
						}
					}
				}
				while(path_end != next_vertex);
			}
			
			// Deal with removing p
			{
				edge_iterator begin = neighbor_range.at(subgraph.p).first,
					end = neighbor_range.at(subgraph.p).second;
				
				// Complete graph on one vertex
				if(begin == end)
				{
					break;
				}
				// Complete graph on two vertices
				else if(begin == --end)
				{
					vertex_descriptor neighbor = get_incident_vertex(subgraph.p, *begin, graph);
					
					// Color the second vertex if uncolored
					if(state.at(neighbor) != State::colored)
					{
						coloring[neighbor] = color_list[neighbor].first();
					}
				}
				// Graph on 3 or more vertices
				{
					auto ordering = embedding[subgraph.p];
					for(auto edge_iter = --end; edge_iter != ++begin; edge_iter--)
					{
						vertex_descriptor neighbor = get_incident_vertex(subgraph.p, *edge_iter, graph);
				
						if(state.at(neighbor) == State::interior)
						{
							state[neighbor] = State::before_p;
						}
						else if(state.at(neighbor) == State::on_face)
						{
					
						}
				
						if(edge_iter == ordering.begin()) edge_iter = ordering.end();
					}
				}
			}
		}
	}*/
	
	template<typename Graph, typename Embedding,
		typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor>
	edge_iterator find_edge(const Graph & graph, const Embedding & embedding, vertex_descriptor vertex,
		const std::pair<edge_iterator, edge_iterator> & range,  vertex_descriptor goal_vertex)
	{
		auto ordering = embedding[vertex];
		for(auto edge_iter = range.first; edge_iter != range.end; edge_iter++)
		{
			if(edge_iter == ordering.end()) edge_iter = ordering.begin();
			vertex_descriptor neighbor = get_incident_vertex(vertex, *edge_iter, graph);
			
			// Done if we found last vertex on path
			if(neighbor == goal_vertex) return edge_iter;
		}
		
		throw std::runtime_error("find_edge: vertex not found in incidence list.");
		return begin;
	}
	
	template<typename Graph, typename Embedding, typename ColorList, typename Coloring,
		typename edge_iterator = typename property_traits<Embedding>::value_type::const_iterator,
		typename vertex_descriptor = typename graph_traits<Graph>::vertex_descriptor,
		typename color_type = typename property_traits<Coloring>::value_type,
		typename state_map = typename std::unordered_map<vertex_descriptor, bool>,
		typename face_location_map = typename std::unordered_map<vertex_descriptor, bool>,
		typename incidence_range_map =
			typename std::unordered_map<vertex_descriptor, std::pair<edge_iterator, edge_iterator> > >
	void path_list_color_recursive(const Graph & graph, const Embedding & embedding, ColorList & color_list,
		Coloring & coloring, state_map & state, face_location_map & face_location,
		incidence_range_map & neighbor_range, vertex_descriptor x, edge_iterator x_begin, edge_iterator x_end,
		vertex_descriptor y, edge_iterator y_begin, edge_iterator y_end, vertex_descriptor p, int before_p,
		int before_x, int before_y)
	{
		typedef std::pair<edge_iterator, edge_iterator> edge_iterator_pair 
		// Update neighbor_range for this subgraph
		neighbor_range[x] = edge_iterator_pair(x_begin, x_end);
		neighbor_range[y] = edge_iterator_pair(y_begin, y_end);
		
		// Base Case 1: K_1
		if(x_begin == x_end)
		{
			// Color x if not colored
			if(!state.at(x)) color[x] = color_list.at(x).first();
			
			return;
		}
		
		// Base Case 2: K_2
		if(++x_begin == x_end)
		{
			// Color x and/or y if not colored
			if(!state.at(x)) color[x] = color_list.at(x).first();
			if(!state.at(y)) color[y] = color_list.at(y).first();
			
			return;
		}
		
		// If colored path doesn't exist, we create one
		if(!state.at(p))
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
					if(state.count(neighbor) == 1 && face_location.at(neighbor) == before_y)
					{
						// Should never find a colored vertex while making path
						if(state.at(neighbor)) throw std::runtime_error("Colored vertex found on face.");
						
						// Check if it may be colored path_color
						for(auto color : color_list[neighbor])
						{
							if(color == path_color)
							{
								// Assign the next vertex in the path
								next_vertex = neighbor;
								
								// Check if we took a chord
								prev_edge_iter = neighbor_range.at(next_vertex).first;
								vertex_descriptor prev_on_cycle =
									get_incident_vertex(next_vertex, *(--prev_edge_iter), graph);
								if(state.at(prev_on_cycle) != COLORED)
								{
									// Find cut point in incidence list
									edge_iterator n_mid = find_edge(graph, embedding, next_vertex,
										neighbor_range.at(next_vertex), p);
									
									// Color lobe
									path_list_color_recursive(graph, embedding, color_list, coloring, state,
										face_location, neighbor_range, next_vertex, n_begin, n_mid, path_end,
										edge_iter, neibhorhood.at(path_end).second, next_vertex,
										(before_y + 1) % 3, (before_y + 2) % 3, before_y);
									
									// Reassign range for end of chord
									neighborhood[next_vertex].first = mid_edge;
									
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
			
			// Interior vertex
			if(state.count(neighbor) == 0)
			{
				// Mark neighbor as on face before p since we are removing p
				state[neighbor] = false;
				face_location[neighbor] = before_p;
				
				// Find edge to p in adjacency list
				edge_iterator n_begin = find_edge(graph, embedding, neighbor,
					neighbor_range.at(neighbor), p), n_end = n_begin;
				
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
					// Set new_p to the next vertex along the face
					new_p = neighbor;
					
					// Reassign neighbor range
					neighbor_range[neighbor].first++;
					if(neighbor_range[neighbor].first == embedding[neighbor].end())
					{
						neighbor_range[neighbor].first = embedding[neighbor].begin();
					}
				}
				
				// Case 1: C[x,p) Cutvertex between x and p inclusive
				else if(face_location[neighbor] == before_p)
				{
					
				}
				// Case 2: C(p,y] Cutvertex between p and y inclusive
				else if(face_location[neighbor] == before_y)
				{
				
				}
				//Case 3: C(y,x) Cutvertex between y and x exclusive
				else if(face_location[neighbor] == before_x)
				{
				
				}
			}
		}
	}
}

#endif