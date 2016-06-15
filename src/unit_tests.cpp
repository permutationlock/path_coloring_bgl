#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <unordered_map>
#include <set>
#include <chrono>
#include <random>

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
#include "poh_color.hpp"
#include "hartman_skrekovski_color.hpp"
#include "draw_tikz_graph.hpp"

using namespace boost;

// Comment line below to hide passing tests.
//#define SHOW_PASSES

// Comment line below to hide tikz drawing printouts
//#define SHOW_VISUALIZATION

// Comment line below to hide test timings.
//#define SHOW_TIMINGS

// Comment line below to hide color list assignment printouts
//#define SHOW_COLOR_LISTS

bool failed=false;

// A class to hold the coordinates of the straight line embedding
struct coord_t {
  std::size_t x;
  std::size_t y;
};

typedef std::chrono::high_resolution_clock Timer;

template<typename Graph>
void make_triangulated(Graph & graph) {
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
void draw_graph_no_color(Graph & graph) {
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
void draw_graph_color(const Graph & graph, const Coloring & coloring) {
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

template<typename index_graph, typename color_map>
void test_path_coloring(const index_graph & graph, const color_map & coloring) {
	typedef typename graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	typedef typename graph_traits<index_graph>::vertex_iterator vertex_iterator;
	typedef typename graph_traits<index_graph>::adjacency_iterator adjacency_iterator;
	typedef typename property_traits<color_map>::value_type color_type;
	
	std::set<vertex_descriptor> visited;
	
	// Iterate over each vertex
	vertex_iterator v_iter, v_end;
	for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
		vertex_descriptor curr_vertex = *v_iter;
		color_type curr_color = coloring[curr_vertex];
		
		// If we have not yet visited this vertex, start bfs on same color vertices
		if(visited.count(curr_vertex) == 0) {
			std::queue<vertex_descriptor> bfs_queue;
			bfs_queue.push(curr_vertex);
			
			std::size_t extra = 1;
			
			while(!bfs_queue.empty()) {
				vertex_descriptor v = bfs_queue.front();
				bfs_queue.pop();
				
				if(visited.count(v) == 0) {
					std::size_t new_neighbor_count = 0, old_neighbor_count = 0;
					
					adjacency_iterator n_iter, n_end;
					for(tie(n_iter, n_end) = adjacent_vertices(v, graph); n_iter != n_end; n_iter++) {
						vertex_descriptor n = *n_iter;
						
						// If this vertex is unvisited
						if(visited.count(n) == 0 && coloring[n] == curr_color) {
							new_neighbor_count++;
							bfs_queue.push(n);
						}
						else if(coloring[n] == curr_color) {
							old_neighbor_count++;
						}
					}
					
					// Find more than one unvisited neighbor, not a path
					// First vertex may have two neighbors, so we use extra
					if(new_neighbor_count > 1 + extra) {
						std::string error = "vertex " + std::to_string(v) + " forms a color "
							+ std::to_string(curr_color) + " branch.";
						throw std::runtime_error(error);
					}
					if(old_neighbor_count > 1) {
						std::string error = "vertex " + std::to_string(v) + " forms a color "
							+ std::to_string(curr_color) + " cycle.";
						throw std::runtime_error(error);
					}
					extra = 0;
					
					visited.insert(v);
				}
			}
		}
	}
}

// Apply Poh algorithm to given graph and verify it works
template<typename index_graph>
void poh_color_test(const index_graph & graph) {
	typedef typename graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	typedef typename graph_traits<index_graph>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<index_graph, vertex_index_t>::type
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
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, color_property_map);
	#endif
	
	// Test correctness of path coloring
	test_path_coloring(graph, color_property_map);
}

// Apply Poh algorithm to given graph and verify it works
template<typename index_graph>
void path_list_color_test(const index_graph & graph, std::size_t num_colors) {
	typedef typename graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	typedef typename graph_traits<index_graph>::vertex_iterator vertex_iterator;
	typedef typename graph_traits<index_graph>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<index_graph, vertex_index_t>::type
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
	typedef iterator_property_map<
			typename std::vector<int>::iterator,
			typename property_map<index_graph, vertex_index_t>::type
		> color_property_map;
	std::vector<int> color_storage(num_vertices(graph));
	color_property_map coloring(color_storage.begin(), get(vertex_index, graph));
	
	// Create vectors to hold ordered outer face
	std::vector<vertex_descriptor> outer_face;
	
	outer_face.push_back(ordering[1]);
	outer_face.push_back(ordering[0]);
	outer_face.push_back(ordering.back());
	
	// Create property map to hold the coloring
	typedef iterator_property_map<
			typename std::vector<std::list<int> >::iterator,
			typename property_map<index_graph, vertex_index_t>::type
		> color_list_property_map;
	std::vector<std::list<int> > color_list_storage(num_vertices(graph));
	color_list_property_map color_list(color_list_storage.begin(), get(vertex_index, graph));
	
	std::mt19937 generator;
	std::uniform_int_distribution<int> distribution(0, num_colors - 1);
	
	// Iterate over each vertex
	vertex_iterator v_iter, v_end;
	for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
		std::vector<int> random_colors(3);
		while(random_colors[0] == random_colors[1] || random_colors[0] == random_colors[2] ||
			random_colors[1] == random_colors[2])
		{
			for(std::size_t i = 0; i < 3; ++i) {
				random_colors[i] = distribution(generator);
			}
		}
		
		#ifdef SHOW_COLOR_LISTS
			std::cout << "color_list[" << *v_iter << "] = { ";
			for(std::size_t i = 0; i < 3; ++i) {
				std::cout << random_colors[i] << ((i != 2) ? ", " : "");
			}
			std::cout << "\n";
		#endif
		
		std::copy(random_colors.begin(), random_colors.end(), std::back_inserter(color_list[*v_iter]));
	}
	
	#ifdef SHOW_TIMINGS
		auto start = Timer::now();
	#endif
	
	// Call path 3-list-color algorithm
	hartman_path_list_color(graph, embedding, color_list, coloring, outer_face.begin(), outer_face.end());
	
	#ifdef SHOW_TIMINGS
		auto end = Timer::now();
		std::cout << "Time = " << (end - start).count() << "\n";
	#endif
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, coloring);
	#endif
	
	// Test correctness of path coloring
	test_path_coloring(graph, coloring);
}

void test_poh_color() {
	std::cout<<"Path 3-coloring"<<std::endl;
	
	{
		// Define graph properties
		typedef adjacency_list
			<	setS,
				vecS,
				undirectedS,
				property<vertex_index_t, std::size_t>,
				property<edge_index_t, std::size_t>
			> index_graph;
		
		typedef erdos_renyi_iterator<minstd_rand, index_graph> ERGen;
		
		boost::minstd_rand gen;
		gen.seed(8573);
		
		for(std::size_t order = 7; order < 100; order++) {
			for(std::size_t seed = 0; seed < 5; seed++) {
				bool found_planar = false;
				std::size_t count = 4;
				
				while(!found_planar) {
					try {
						// Construct a random trriangulated graph
						//std::cout << "Generating graph.\n";
						index_graph graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
				
						++count;
						
						//std::cout << "Triangulating graph.\n";
						make_triangulated(graph);
				
						found_planar = true;
				
						//draw_graph_no_color(graph);
				
						//std::cout << "Testing planarity.\n";
						poh_color_test(graph);
				
						#ifdef SHOW_PASSES
							std::cout<<"    PASS " << order << " vertex path 3-color."<<std::endl;
						#endif
					}
					catch(std::logic_error error) {
						// Generated a non-planar graph, ignore this case
					}
					catch(std::exception& error) {
						std::cout<<"    FAIL " << order << " vertex path 3-color ("<<error.what()<<")."<<std::endl;
						failed=true;
					}
					catch(...) {
						std::cout<<"    FAIL " << order << " vertex path 3-color (unknown error)."<<std::endl;
						failed=true;
					}
				}
			}
		}
	}
}

void test_list_path_color()
{
	std::cout<<"Path 3-list-coloring"<<std::endl;
	
	{
		// Define graph properties
		typedef adjacency_list
			<	setS,
				vecS,
				undirectedS,
				property<vertex_index_t, std::size_t>,
				property<edge_index_t, std::size_t>
			> index_graph;
	
		typedef erdos_renyi_iterator<minstd_rand, index_graph> ERGen;
		
		boost::minstd_rand gen;
		gen.seed(8573);
		
		for(std::size_t order = 7; order < 100; order++) {
			for(std::size_t colors = 3; colors < 9; ++colors) {
				for(std::size_t seed = 0; seed < 5; seed++) {
					bool found_planar = false;
					std::size_t count = 4;
					
					while(!found_planar) {
						try {
							// Construct a random trriangulated graph
							//std::cout << "Generating graph.\n";
							index_graph graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
							
							++count;
							
							//std::cout << "Triangulating graph.\n";
							make_triangulated(graph);
					
							found_planar = true;
					
							//draw_graph_no_color(graph);
					
							//std::cout << "Testing planarity.\n";
							path_list_color_test(graph, colors);
					
							#ifdef SHOW_PASSES
								std::cout<<"    PASS " << order << " vertex, " << colors << " colors, test "
									<< seed << " path 3-list-color."<<std::endl;
							#endif
						}
						catch(std::logic_error error) {
							// Generated a non-planar graph, ignore this case
						}
						catch(std::exception& error) {
							std::cout<<"    FAIL " << order << " vertex, " << colors << " colors, test "
									<< seed << " path 3-list-color ("<<error.what()<<")."<<std::endl;
							failed=true;
						}
						catch(...) {
							std::cout<<"    FAIL " << order << " vertex, " << colors << " colors, test "
									<< seed << " path 3-list-color (unknown error)."<<std::endl;
							failed=true;
						}
					}
				}
			}
		}
	}
}

int main() {
	//test_poh_color();
	test_list_path_color();

	if(failed)
		std::cout<<"THERE ARE FAILING TESTS"<<std::endl;
	else
		std::cout<<"ALL TESTS PASSED"<<std::endl;
	
	return 0;
}