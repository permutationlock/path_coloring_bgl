#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>
#include <string>
#include <unordered_map>
#include <set>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_connected.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>

#include "path_coloring.hpp"
#include "draw_tikz_graph.hpp"

using namespace boost;

//Comment line below to hide passing tests.
#define SHOW_PASSES

bool failed=false;

// A class to hold the coordinates of the straight line embedding
struct coord_t
{
  std::size_t x;
  std::size_t y;
};


template<typename AdjacencyGraph, typename Coloring>
void test_path_coloring(const AdjacencyGraph & graph, const Coloring & coloring)
{
	typedef graph_traits<AdjacencyGraph> GraphTraits;
	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::vertex_iterator vertex_iterator;
	typedef typename GraphTraits::adjacency_iterator adjacency_iterator;
	typedef typename boost::property_traits<Coloring>::value_type color_type;
	
	std::set<vertex_descriptor> visited;
	
	// Iterate over each vertex
	vertex_iterator v_iter, v_end;
    for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++)
    {
		vertex_descriptor curr_vertex = *v_iter;
		color_type curr_color = coloring[curr_vertex];
		
		// If we have not yet visited this vertex, start bfs on same color vertices
		if(visited.count(curr_vertex) == 0)
		{
			std::queue<vertex_descriptor> bfs_queue;
			visited.insert(curr_vertex);
			bfs_queue.push(curr_vertex);
			
			std::size_t extra = 1;
			
			while(!bfs_queue.empty())
			{
				vertex_descriptor v = bfs_queue.front();
				bfs_queue.pop();
				
				if(visited.count(curr_vertex) == 0)
				{
					std::size_t count = 0;
					
					adjacency_iterator n_iter, n_end;
					for(tie(n_iter,n_end) = adjacent_vertices(curr_vertex, graph);
						n_iter != n_end; n_iter++)
					{
						vertex_descriptor n = *n_iter;
						
						// If this vertex is unvisited
						if(visited.count(n) == 0 && coloring[n] == curr_color)
						{
							count++;
							visited.insert(n);
							bfs_queue.push(n);
						}
						
						
					}
					
					// Find more than one unvisited neighbor, not a path
					// First vertex may have two neighbors, so we use extra
					if(count > 1 + extra)
					{
						std::string error = "vertex " + std::to_string(v) + " forms a color "
							+ std::to_string(curr_color) + " branch or cycle.";
						throw std::runtime_error(error);
					}
					extra = 0;
				}
			}
		}
	}
}

void test_poh_color()
{
	std::cout<<"Poh coloring"<<std::endl;
	
	{
		// Define graph properties
		typedef adjacency_list
			<	vecS,
				vecS,
				undirectedS,
				property<vertex_index_t, int>,
				property<edge_index_t, int>
			> graph;



		// Define the storage type for the planar embedding
		typedef std::vector< std::vector< graph_traits<graph>::edge_descriptor > > 
			embedding_storage_t;

		typedef iterator_property_map
			< embedding_storage_t::iterator, 
				property_map<graph, vertex_index_t>::type
			> embedding_t;
		
		try
		{
			// Construct triangulated plane graph on 9 vertices
			graph g(9);
			
			add_edge(0,1,g);
			add_edge(0,8,g);
			add_edge(8,1,g);
			add_edge(2,1,g);
			add_edge(3,1,g);
			add_edge(2,8,g);
			add_edge(3,0,g);
			add_edge(2,5,g);
			add_edge(2,4,g);
			add_edge(2,3,g);
			add_edge(3,6,g);
			add_edge(3,4,g);
			add_edge(4,5,g);
			add_edge(4,6,g);
			add_edge(5,8,g);
			add_edge(5,7,g);
			add_edge(5,6,g);
			add_edge(6,0,g);
			add_edge(6,7,g);
			add_edge(7,8,g);
			add_edge(7,0,g);
			
			// Create the planar embedding
			embedding_storage_t embedding_storage(num_vertices(g));
			embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));

			boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
				boyer_myrvold_params::embedding = embedding);
				
			/*// Print embedding ordered adjacency lists
			for(int i=0; i<9; i++)
			{                      
				std::cout << "Adjacency list for " << i << ".\n";       
				for(auto iter = embedding[i].begin(); iter != embedding[i].end(); iter++)
				{
					std::cout << "{" << source(*iter, g) << ", " << target(*iter, g)
						<< "} adj = " << get_incident_vertex(i,*iter,g) << "\n"; 
				}
				std::cout << "\n";
			}*/
			
			// Find a canonical ordering
			std::vector<graph_traits<graph>::vertex_descriptor> ordering;
			planar_canonical_ordering(g, embedding, std::back_inserter(ordering));
			
			// Create property map to hold the coloring
			std::map<graph::vertex_descriptor,int> color_map;
			associative_property_map< std::map<graph::vertex_descriptor, int> >
				color_property_map(color_map);
			
			// Create vectors to hold chordless paths
			std::vector<graph::vertex_descriptor> p;
			std::vector<graph::vertex_descriptor> q;
			
			// Initialize path p with top vertex of outer face and color 1
			p.push_back(ordering.back());
			color_property_map[ordering.back()] = 1;

			// Initialize path q with bottom vertices of outer face and color 2
			q.push_back(ordering[0]);
			q.push_back(ordering[1]);
			color_property_map[ordering[0]] = 2;
			color_property_map[ordering[1]] = 2;
			
			// Remaining color is 3
			int color = 0;
			
			// Call Poh algorithm
			poh_path_color(g, embedding, p.begin(), p.end(), q.begin(), q.end(),
				color_property_map, color);
			
			// Test correctness of path coloring
			test_path_coloring(g, color_property_map);
			
			/*//Set up a property map to hold the mapping from vertices to coord_t's
			typedef std::vector< coord_t > straight_line_drawing_storage_t;
			typedef boost::iterator_property_map
				<	straight_line_drawing_storage_t::iterator, 
					property_map<graph, vertex_index_t>::type 
				> straight_line_drawing_t;

			straight_line_drawing_storage_t straight_line_drawing_storage
				(num_vertices(g));
			straight_line_drawing_t straight_line_drawing
				(straight_line_drawing_storage.begin(), 
				get(vertex_index,g));



			// Compute the straight line drawing
			chrobak_payne_straight_line_drawing(g, 
				embedding, 
				ordering.begin(),
				ordering.end(),
				straight_line_drawing
				);
			
			std::cout << draw_tikz_graph(g, color_property_map, straight_line_drawing) << "\n";*/
			
			#ifdef SHOW_PASSES
				std::cout<<"    PASS 9 vertex Poh path color."<<std::endl;
			#endif
		}
		catch(std::exception& error)
		{
			std::cout<<"    FAIL 9 vertex Poh path color ("<<error.what()<<")."<<std::endl;
			failed=true;
		}
		catch(...)
		{
			std::cout<<"    FAIL 9 vertex Poh path color (unknown error)."<<std::endl;
			failed=true;
		}
	}
}

int main()
{
	test_poh_color();

	if(failed)
		std::cout<<"THERE ARE FAILING TESTS"<<std::endl;
	else
		std::cout<<"ALL TESTS PASSED"<<std::endl;


	return 0;
}