/*
 * path_coloring.h
 * Author: Aven Bross
 * 
 * Implementation of several algorithms for path coloring graphs.
 */

#ifndef __PATH_COLORING_H
#define __PATH_COLORING_H

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

	template<typename Graph>
	typename graph_traits<Graph>::vertex_descriptor
		getIncidentVertex(typename graph_traits<Graph>::vertex_descriptor vertex,
		typename graph_traits<Graph>::edge_descriptor edge, const Graph & graph)
	{
		auto adj = source(edge, graph);
		if(adj == vertex)
			adj = target(edge, graph);
		return adj;
	}

	/*
	 * Preconditions: graph is a weakly triangulated planar graph the neighbors of each
	 * vertex arranged in clockwise order for some embedding. p and q are paths with
	 * p above q (in clockwise orientation scheme) such that joining the beginnings and
	 * ends of the paths forms a chordless cycle, and with p and q each colored a
	 * different color. new_color is a third color different from the color of p and q.
	 * 
	 * Postconditions: returns a copy of g with the face enclosed by the cycle pq
	 * 3-colored such that each color class induces a disjoint union of paths.
	 */
	template<typename Graph, typename Embedding, typename Coloring, typename VertexIter>
	void poh_path_color(const Graph & graph, const Embedding & embedding,
		VertexIter p_0, VertexIter p_1, VertexIter q_0, VertexIter q_1, Coloring & coloring,
		typename property_traits<Coloring>::value_type new_color)
	{
		typedef graph_traits<Graph> GraphTraits;
		typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
		typedef typename GraphTraits::edge_descriptor edge_descriptor;
		
		// Base case: K_2
		if(p_1 - p_0 == 1 && q_1 - q_0 == 1)
		{
			return;
		}
	
		// Find inner triangulation vertex for beginning of paths
		vertex_descriptor t_0 = *p_0;
		
		// Iterate incident edges in clockwise embedded order
		for(auto iter = embedding[*p_0].begin(); iter != embedding[*p_0].end(); iter++)
		{
			// Locate vertex on the end of edge
			vertex_descriptor neighbor = getIncidentVertex(*p_0, *iter, graph);
			
			if(neighbor == *q_0)	// Find begining vertex of second path
			{
				// Triangulation vertex is one back counterclockwise
				if(iter == embedding[*p_0].begin())
					t_0 = getIncidentVertex(*p_0, *(embedding[*p_0].end() - 1), graph);	
				else
					t_0 = getIncidentVertex(*p_0, *(iter - 1), graph);
				break;
			}
		}
	
		if(t_0 == *p_0) throw std::runtime_error("no triangulation vertex found at start of paths.");
	
		// Case 1: First triangle hits p or q
		if(p_1 - p_0 > 1 && t_0 == *(p_0 + 1))
		{
			std::cout << "  --First triangle hit p--\n";
			// Removing triangle removes first vertex of p
			poh_path_color(graph, embedding, p_0 + 1, p_1, q_0, q_1, coloring, new_color);
			
			return;
		}
		else if(q_1 - q_0 > 1  && t_0 == *(q_0 + 1))
		{
			std::cout << "  --First triangle hit q--\n";
			// Removing triangle removes first vertex of q
			poh_path_color(graph, embedding, p_0, p_1, q_0 + 1, q_1, coloring, new_color);
		
			return;
		}
	
		// Find inner triangulation vertex for end of paths
		vertex_descriptor t_1 = *q_1;
	
		// Iterate incident edges in clockwise embedded order
		for(auto iter = embedding[*(q_1 - 1)].begin(); iter != embedding[*(q_1 - 1)].end(); iter++)
		{
			// Locate vertex on the end of edge
			vertex_descriptor neighbor = getIncidentVertex(*(q_1 - 1), *iter, graph);
			
			if(neighbor == *(p_1 - 1))	// Find begining vertex of second path
			{
				// Triangulation vertex is one back counterclockwise
				if(iter == embedding[*(q_1 - 1)].begin())
					t_1 = getIncidentVertex(*(q_1 - 1), *(embedding[*(q_1 - 1)].end() - 1), graph);	
				else
					t_1 = getIncidentVertex(*(q_1 - 1), *(iter - 1), graph);
				break;
			}
		}
	
		if(t_1 == *q_1) throw std::runtime_error("no triangulation vertex found at end of paths.");
	
		// Case 2: Second triangle hits p or q
		if(p_1 - p_0 > 1 && t_1 == *(p_1 - 2))
		{	
			std::cout << "  --Second triangle hit p--\n";
			// Removing triangle removes last vertex of p
			poh_path_color(graph, embedding, p_0, p_1 - 1, q_0, q_1, coloring, new_color);
		
			return;
		}
		else if(q_1 - q_0 > 1  && t_1 == *(q_1 - 2))
		{
			std::cout << "  --Second triangle hit q--\n";
			// Removing triangle removes last vertex of q
			poh_path_color(graph, embedding, p_0, p_1, q_0, q_1 - 1, coloring, new_color);
		
			return;
		}
		
		// Case 3: Search for path-bridging edge
		{
			std::cout << "  --Finding bridging edge--\n";
			std::unordered_map<vertex_descriptor,VertexIter> p_map;
		
			// Mark all edges from vertices on p
			for(auto p_iter = p_0; p_iter != p_1; p_iter++)
			{
				p_map[*p_iter] = p_iter;
			}
		
			for (auto q_iter = q_0; q_iter != q_1; ++q_iter)
			{
				for (edge_descriptor inc_edge : embedding[*q_iter])
				{
					vertex_descriptor q_nbr = getIncidentVertex(*q_iter, inc_edge, graph);
					
					if(p_map.count(q_nbr) != 0)
					{
						// Check if our edge is one of the outside edges
						if(q_iter == q_0 && p_map[q_nbr] == p_0)
							continue;
						if(q_iter == q_1 - 1 && p_map[q_nbr] == p_1 - 1)
							continue;
					
						// Recurse on "left" and "right" halves
						poh_path_color(graph, embedding, p_0, p_map[q_nbr] + 1, q_0,q_iter + 1,
							coloring, new_color);
						poh_path_color(graph, embedding, p_map[q_nbr], p_1, q_iter, q_1,
							coloring, new_color);
					
						return;
					}
				}
			}
		}
		
		// Case 4: Find splitting path
		{
			std::cout << "  --Finding splitting path--\n";
			std::unordered_map<vertex_descriptor,vertex_descriptor> parent_map;
			std::queue<vertex_descriptor> bfs_queue;
		
			// Perform BFS from t_1 inside the cycle and locate t_1-t_0 path
			bfs_queue.push(t_1);
			parent_map[t_1] = t_1;
	
			for(auto p_iter = p_0; p_iter != p_1; p_iter++)
			{
				parent_map[*p_iter] = *p_iter;
			}
		
			for(auto q_iter = q_0; q_iter != q_1; q_iter++)
			{
				parent_map[*q_iter] = *q_iter;
			}
	
			vertex_descriptor curr_vertex;
		
			while(!bfs_queue.empty() && curr_vertex != t_0)
			{
				curr_vertex = bfs_queue.front();
				bfs_queue.pop();
		
				for(edge_descriptor inc_edge : embedding[curr_vertex])
				{
					vertex_descriptor neighbor = getIncidentVertex(curr_vertex, inc_edge, graph);
					if(parent_map.count(neighbor) == 0)
					{
						bfs_queue.push(neighbor);
						parent_map[neighbor] = curr_vertex;
					}
				}
			}
		
			// A t_1 - t_0 path should be guaranteed, throw if not found
			if(curr_vertex != t_0)
			{
				throw std::runtime_error("No splitting path or bridging edge.");
			}
		
			// Backtrack to construct t_0 - t_1 path and color using new_color
			std::vector<vertex_descriptor> splitting_path;
			splitting_path.push_back(curr_vertex);
			coloring[curr_vertex] = new_color;
		
			while(parent_map[curr_vertex] != curr_vertex)
			{
				curr_vertex = parent_map[curr_vertex];
				splitting_path.push_back(curr_vertex);
				coloring[curr_vertex] = new_color;
			}
		
			// Recurse on top and bottom halves, now coloring with correct colors
			poh_path_color(graph, embedding, p_0, p_1, splitting_path.begin(),
				splitting_path.end(), coloring, coloring[*q_0]);
			poh_path_color(graph, embedding, splitting_path.begin(), splitting_path.end(),
				q_0, q_1, coloring, coloring[*p_0]);
		}
	}
}

#endif