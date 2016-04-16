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
	template<typename Graph, typename Embedding, typename Coloring, typename VertexIter>
	void poh_path_color(const Graph & graph, const Embedding & embedding,
		VertexIter p_begin, VertexIter p_end, VertexIter q_begin, VertexIter q_end, Coloring & coloring,
		typename property_traits<Coloring>::value_type new_color)
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
	template<typename Graph, typename Embedding, typename ColorList, typename Coloring,
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
		std::vector<VertexIter> path;
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
				    hartman_remove_path_vertex(graph, embedding, p_iter, path[i+1], p_iter, p_iter, new_path.begin(),
				        new_path.end(), color_list, coloring, path_color);
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
		    hartman_remove_path_vertex(graph, embedding, path.back(), p_end(), q_begin(), q_end(), new_path.begin(),
		        new_path.end(), color_list, coloring, path_color);
		}
	}
	
	template<typename Graph, typename Embedding, typename ColorList, typename Coloring,
		typename VertexIter>
	void hartman_list_color_remove_path(const Graph & graph, const Embedding & embedding,
		VertexList top, VertexList bottom, VertexList::iterator p_begin, VertexList::iterator p_end,
		const ColorList & color_list, Coloring & coloring)
	{
		typedef typename graph_traits<Graph>::vertex_descriptor vertex_descriptor;
		typedef typename property_traits<Coloring>::value_type color_type;
    }
}

#endif