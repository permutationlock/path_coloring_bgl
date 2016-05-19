#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <unordered_map>
#include <set>

// Basic graph headers
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Random graph headers
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>

// Planar graph headers
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/make_connected.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>

// Local project headers
#include "path_coloring.hpp"
#include "draw_tikz_graph.hpp"

using namespace boost;

// Comment line below to hide passing tests.
#define SHOW_PASSES

// Comment line below to hide tikz drawing printouts
#define SHOW_VISUALIZATION

bool failed=false;

// A class to hold the coordinates of the straight line embedding
//struct coord_t
//{
//  std::size_t x;
//  std::size_t y;
//};

template<typename Graph>
void make_triangulated(Graph & graph)
{
	typedef typename graph_traits<Graph>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<Graph, vertex_index_t>::type
		> embedding_t;

	make_connected(graph);

	//Initialize the interior edge index
	typename property_map<Graph, edge_index_t>::type e_index = get(edge_index, graph);
	typename graph_traits<Graph>::edges_size_type edge_count = 0;
	typename graph_traits<Graph>::edge_iterator ei, ei_end;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);

	// Create the planar embedding
	embedding_storage_t embedding_storage(num_vertices(graph));
	embedding_t embedding(embedding_storage.begin(), get(vertex_index, graph));

	if(!boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
						   boyer_myrvold_params::embedding = embedding))
	{
		std::string error = "Non-planar graph.";
		throw std::logic_error(error);
	}            

	// Re-initialize the edge index, since we just added a few edges
	edge_count = 0;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);

	make_biconnected_planar(graph, embedding);

	boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
						   boyer_myrvold_params::embedding = embedding);

	// Re-initialize the edge index, since we just added a few edges
	edge_count = 0;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);

	make_maximal_planar(graph, embedding);

	// Re-initialize the edge index, since we just added a few edges
	edge_count = 0;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);
}

template<typename Graph>
void draw_graph_no_color(Graph & graph)
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
	
	std::map<vertex_descriptor,int> color_map;
  	boost::associative_property_map< std::map<vertex_descriptor, int> >
		color_property_map(color_map);
		
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
	
	std::cout << draw_tikz_graph(graph, color_property_map, straight_line_drawing) << "\n";
}

template<typename Graph, typename Coloring>
void draw_graph_color(const Graph & graph, const Coloring & coloring)
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
	
	std::cout << draw_tikz_graph(graph, coloring, straight_line_drawing) << "\n";
}

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
			bfs_queue.push(curr_vertex);
			
			std::size_t extra = 1;
			
			while(!bfs_queue.empty())
			{
				vertex_descriptor v = bfs_queue.front();
				bfs_queue.pop();
				
				if(visited.count(curr_vertex) == 0)
				{
					std::size_t new_neighbor_count = 0, old_neighbor_count = 0;
					
					adjacency_iterator n_iter, n_end;
					for(tie(n_iter,n_end) = adjacent_vertices(curr_vertex, graph);
						n_iter != n_end; n_iter++)
					{
						vertex_descriptor n = *n_iter;
						
						// If this vertex is unvisited
						if(visited.count(n) == 0 && coloring[n] == curr_color)
						{
							new_neighbor_count++;
							visited.insert(n);
							bfs_queue.push(n);
						}
						else if(coloring[n] == curr_color)
						{
							old_neighbor_count++;
						}
					}
					
					// Find more than one unvisited neighbor, not a path
					// First vertex may have two neighbors, so we use extra
					if(new_neighbor_count > 1 + extra)
					{
						std::string error = "vertex " + std::to_string(v) + " forms a color "
							+ std::to_string(curr_color) + " branch.";
						throw std::runtime_error(error);
					}
					if(old_neighbor_count > 1)
					{
						std::string error = "vertex " + std::to_string(v) + " forms a color "
							+ std::to_string(curr_color) + " cycle.";
						throw std::runtime_error(error);
					}
					extra = 0;
				}
			}
		}
	}
}

// Apply Poh algorithm to given graph and verify it works
template<typename AdjacencyGraph>
void poh_color_test(const AdjacencyGraph & graph)
{
	typedef typename graph_traits<AdjacencyGraph>::vertex_descriptor vertex_descriptor;
	typedef typename graph_traits<AdjacencyGraph>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<AdjacencyGraph, vertex_index_t>::type
		> embedding_t;
	
	// Create the planar embedding
	embedding_storage_t embedding_storage(num_vertices(graph));
	embedding_t embedding(embedding_storage.begin(), get(vertex_index, graph));

	boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
		boyer_myrvold_params::embedding = embedding);
	
	// Find a canonical ordering
	std::vector<vertex_descriptor> ordering;
	planar_canonical_ordering(graph, embedding, std::back_inserter(ordering));
	
	// Create property map to hold the coloring
	std::map<vertex_descriptor, int> color_map;
	associative_property_map< std::map<vertex_descriptor, int> >
		color_property_map(color_map);
	
	// Create vectors to hold chordless paths
	std::vector<vertex_descriptor> p;
	std::vector<vertex_descriptor> q;
	
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
	poh_path_color(graph, embedding, p.begin(), p.end(), q.begin(), q.end(),
		color_property_map, color);
	
	// Test correctness of path coloring
	test_path_coloring(graph, color_property_map);
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, color_property_map);
	#endif
}

// Apply Poh algorithm to given graph and verify it works
template<typename AdjacencyGraph>
void path_list_color_test(const AdjacencyGraph & graph)
{
	typedef typename graph_traits<AdjacencyGraph>::vertex_descriptor vertex_descriptor;
	typedef typename graph_traits<AdjacencyGraph>::vertex_iterator vertex_iterator;
	typedef typename graph_traits<AdjacencyGraph>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<AdjacencyGraph, vertex_index_t>::type
		> embedding_t;
	
	// Create the planar embedding
	embedding_storage_t embedding_storage(num_vertices(graph));
	embedding_t embedding(embedding_storage.begin(), get(vertex_index, graph));

	boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
		boyer_myrvold_params::embedding = embedding);
	
	// Find a canonical ordering
	std::vector<vertex_descriptor> ordering;
	planar_canonical_ordering(graph, embedding, std::back_inserter(ordering));
	
	// Create property map to hold the coloring
	std::map<vertex_descriptor, int> color_map;
	associative_property_map< std::map<vertex_descriptor, int> >
		coloring(color_map);
	
	// Create vectors to hold chordless paths
	std::vector<vertex_descriptor> outer_face;
	
	outer_face.push_back(ordering[1]);
	outer_face.push_back(ordering[0]);
	outer_face.push_back(ordering.back());
	
	std::unordered_map<vertex_descriptor, std::list<int> > color_list;
	
	// Iterate over each vertex
	vertex_iterator v_iter, v_end;
	for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++)
	{
		color_list[*v_iter] = {0, 1, 2};
	}
	
	// Call Poh algorithm
	path_list_color(graph, embedding, color_list, coloring, outer_face.begin(), outer_face.end());
	
	// Test correctness of path coloring
	test_path_coloring(graph, coloring);
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, coloring);
	#endif
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
			> Graph;
		
		typedef erdos_renyi_iterator<minstd_rand, Graph> ERGen;
		
		boost::minstd_rand gen;
		
		for(std::size_t order = 10; order < 25; order++)
		{
			bool found_planar = false;
			while(!found_planar)
			{
				try
				{
					// Construct a random trriangulated graph
					//std::cout << "Generating graph.\n";
					Graph graph(ERGen(gen, order, 2 * order - 4), ERGen(), order);
					
					//std::cout << "Triangulating graph.\n";
					make_triangulated(graph);
					
					found_planar = true;
					
					//draw_graph_no_color(graph);
					
					//std::cout << "Testing planarity.\n";
					poh_color_test(graph);
					
					#ifdef SHOW_PASSES
						std::cout<<"    PASS " << order << " vertex Poh path color."<<std::endl;
					#endif
				}
				catch(std::logic_error error)
				{
					// Generated a non-planar graph, ignore this case
				}
				catch(std::exception& error)
				{
					std::cout<<"    FAIL " << order << " vertex Poh path color ("<<error.what()<<")."<<std::endl;
					failed=true;
				}
				catch(...)
				{
					std::cout<<"    FAIL " << order << " vertex Poh path color (unknown error)."<<std::endl;
					failed=true;
				}
			}
		}
	}
}

void test_list_path_color()
{
	std::cout<<"List path coloring"<<std::endl;
	
	{
		// Define graph properties
		typedef adjacency_list
			<	vecS,
				vecS,
				undirectedS,
				property<vertex_index_t, int>,
				property<edge_index_t, int>
			> Graph;
	
		try
		{
			// Construct a random trriangulated graph
			//std::cout << "Generating graph.\n";
			Graph graph(10);
	
			//std::cout << "Triangulating graph.\n";
			make_triangulated(graph);
	
			draw_graph_no_color(graph);
	
			//std::cout << "Testing planarity.\n";
			path_list_color_test(graph);
	
			#ifdef SHOW_PASSES
				std::cout<<"    PASS " << 3 << " vertex list path color."<<std::endl;
			#endif
		}
		catch(std::logic_error error)
		{
			// Generated a non-planar graph, ignore this case
		}
		catch(std::exception& error)
		{
			std::cout<<"    FAIL " << 3 << " vertex list path color ("<<error.what()<<")."<<std::endl;
			failed=true;
		}
		catch(...)
		{
			std::cout<<"    FAIL " << 3 << " vertex list path color (unknown error)."<<std::endl;
			failed=true;
		}
	}
}

int main()
{
	test_list_path_color();

	if(failed)
		std::cout<<"THERE ARE FAILING TESTS"<<std::endl;
	else
		std::cout<<"ALL TESTS PASSED"<<std::endl;


	return 0;
}