/*
 * unit_tests.cpp
 * Author: Aven Bross
 *
 * Unit testing suite for path coloring.
 */

// STL headers
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
#include "../path_coloring/poh_color.hpp"
#include "../path_coloring/poh_color_bfs.hpp"
#include "../path_coloring/augmented_embedding.hpp"
#include "../path_coloring/hartman_skrekovski_choose.hpp"
#include "../visualization/draw_tikz_graph.hpp"

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

typedef std::chrono::high_resolution_clock nanosecond_timer;

template<typename index_graph>
void make_triangulated(index_graph & graph) {
	typedef typename graph_traits<index_graph>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<index_graph, vertex_index_t>::type
		> embedding_t;

	make_connected(graph);
	
	//Initialize the interior edge index
	typename property_map<index_graph, edge_index_t>::type e_index = get(edge_index, graph);
	typename graph_traits<index_graph>::edges_size_type edge_count = 0;
	typename graph_traits<index_graph>::edge_iterator ei, ei_end;
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

	make_biconnected_planar(graph, embedding);
	
	// Re-initialize the edge index, since we just added a few edges
	edge_count = 0;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);
	
	boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
						   boyer_myrvold_params::embedding = embedding);

	make_maximal_planar(graph, embedding);
	
	// Re-initialize the edge index, since we just added a few edges
	edge_count = 0;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);

}

template<typename index_graph>
void draw_graph_no_color(index_graph & graph) {
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
			typename property_map<index_graph, vertex_index_t>::type 
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

template<typename index_graph, typename color_map>
void draw_graph_color(const index_graph & graph, const color_map & coloring) {
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
		
	//Set up a property map to hold the mapping from vertices to coord_t's
	typedef std::vector< coord_t > straight_line_drawing_storage_t;
	typedef boost::iterator_property_map
		< straight_line_drawing_storage_t::iterator, 
			typename property_map<index_graph, vertex_index_t>::type 
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

template<
		typename graph_t, typename planar_embedding_t,
		typename augmented_embedding_t
	>
void test_augmented_embedding(
		const graph_t & graph, const planar_embedding_t & planar_embedding,
		const augmented_embedding_t & augmented_embedding
	)
{
	// Type definitions
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename graph_traits<graph_t>::vertex_iterator vertex_iterator_t;
	
	// Iterate over each vertex
	vertex_iterator_t v_iter, v_end;
	for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
		vertex_t v = *v_iter;
		
		auto edge_iter = planar_embedding[v].begin();
		auto adjacency_node_iter = augmented_embedding[v].begin();
		
		while(edge_iter != planar_embedding[v].end() &&
			adjacency_node_iter != augmented_embedding[v].end())
		{
			vertex_t u_0 = get_incident_vertex(v, *edge_iter, graph);
			vertex_t u_1 = adjacency_node_iter -> vertex;
			
			if(u_0 != u_1) {
				std::string error = "vertex " + std::to_string(v) + " has neighbor "
					+ std::to_string(u_0) + " in embedding, but neighbor "
					+ std::to_string(u_1) + " in augmented_embedding.";
				throw std::runtime_error(error);
			}
			
			vertex_t v_0 = adjacency_node_iter -> iterator -> vertex;
			
			if(v_0 != v) {
				std::string error = "edge " + std::to_string(v) + " to "
					+ std::to_string(u_0) + " back iterator points to "
					+ std::to_string(v_0) + " in augmented_embedding.";
				throw std::runtime_error(error);
			}
			
			++edge_iter;
			++adjacency_node_iter;
		}
		
		if(edge_iter != planar_embedding[v].end()) {
			std::string error = "vertex " + std::to_string(v)
				+ " augmented adjacency list too long.";
			throw std::runtime_error(error);
		}
		else if(adjacency_node_iter != augmented_embedding[v].end()) {
			std::string error = "vertex " + std::to_string(v)
				+ " augmented adjacency list too short.";
			throw std::runtime_error(error);
		}
	}
}

template<typename graph_t, typename planar_embedding_t>
void test_augmented_embedding_construction(
		const graph_t & graph, const planar_embedding_t & planar_embedding
	)
{
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	
	struct adjacency_node_t {
		vertex_t vertex;
		typename std::vector<adjacency_node_t>::iterator iterator;
	};
	
	typedef boost::iterator_property_map<
			typename std::vector<std::vector<adjacency_node_t>>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> augmented_embedding_t;
	
	std::vector<std::vector<adjacency_node_t>>
		adjacency_node_storage(boost::num_vertices(graph));
	augmented_embedding_t augmented_embedding(
			adjacency_node_storage.begin(),
			boost::get(boost::vertex_index, graph)
		);
	
	augment_embedding(graph, planar_embedding, augmented_embedding);
	
	test_augmented_embedding(graph, planar_embedding, augmented_embedding);
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
void poh_color_bfs_test(const index_graph & graph) {
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
	
	// Initialize path p with top vertex of outer face
	p.push_back(ordering[0]);

	// Initialize path q with bottom vertices of outer face
	q.push_back(ordering[1]);
	q.push_back(ordering.back());
	
	// Call Poh algorithm
	#ifdef SHOW_TIMINGS
		auto start = nanosecond_timer::now();
		
		for(std::size_t i = 0; i < 1000; ++i) {
			poh_color_bfs(graph, embedding, color_property_map, p.begin(), p.end(), q.begin(), q.end(), 1, 2, 3);
		}
	
		auto end = nanosecond_timer::now();
		std::cout << "Time = " << (end - start).count() / 1000 << "ns\n";
	#else
		poh_color_bfs(graph, embedding, color_property_map, p.begin(), p.end(), q.begin(), q.end(), 1, 2, 3);
	#endif
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, color_property_map);
	#endif
	
	// Test correctness of path coloring
	test_path_coloring(graph, color_property_map);
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
	
	// Initialize path p with top vertex of outer face
	p.push_back(ordering[0]);

	// Initialize path q with bottom vertices of outer face
	q.push_back(ordering[1]);
	q.push_back(ordering.back());
	
	// Call Poh algorithm
	#ifdef SHOW_TIMINGS
		auto start = nanosecond_timer::now();
		
		for(std::size_t i = 0; i < 1000; ++i) {
			poh_color_bfs(graph, embedding, color_property_map, p.begin(), p.end(), q.begin(), q.end(), 1, 2, 3);
		}
	
		auto end = nanosecond_timer::now();
		std::cout << "Time = " << (end - start).count() / 1000 << "ns\n";
	#else
		poh_color(graph, embedding, color_property_map, p.begin(), p.end(), q.begin(), q.end(), 1, 2, 3);
	#endif
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, color_property_map);
	#endif
	
	// Test correctness of path coloring
	test_path_coloring(graph, color_property_map);
}


// Apply Hartman-Skrekovski algorithm to given graph and verify it works
template<typename graph_t>
void path_choose_test(const graph_t & graph, std::size_t num_colors) {
	typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename graph_traits<graph_t>::vertex_iterator vertex_iterator_t;
	typedef typename graph_traits<graph_t>::edge_descriptor edge_descriptor;
	
	// Define the storage type for the planar embedding
	typedef std::vector<
			std::vector<edge_descriptor>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<graph_t, vertex_index_t>::type
		> embedding_t;
	
	// Define the storage type for nodes in the augmented embedding
	struct adjacency_node_t {
		vertex_t vertex;
		typename std::vector<adjacency_node_t>::iterator iterator;
	};
	
	typedef boost::iterator_property_map<
			typename std::vector<std::vector<adjacency_node_t>>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> augmented_embedding_t;
	
	// Create the planar embedding
	embedding_storage_t embedding_storage(num_vertices(graph));
	embedding_t embedding(embedding_storage.begin(), get(vertex_index, graph));

	boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
		boyer_myrvold_params::embedding = embedding);
	
	// Find a canonical ordering
	std::vector<vertex_t> ordering;
	planar_canonical_ordering(graph, embedding, std::back_inserter(ordering));
	
	std::vector<std::vector<adjacency_node_t>>
		adjacency_node_storage(boost::num_vertices(graph));
	augmented_embedding_t augmented_embedding(
			adjacency_node_storage.begin(),
			boost::get(boost::vertex_index, graph)
		);
	
	vertex_iterator_t v_iter, v_end;
	
	// Reserve deg(v) space for each v in G
	for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end;
			v_iter++
		)
	{
		vertex_t v = *v_iter;
		augmented_embedding[v].reserve(out_degree(v, graph));
	}
	
	augment_embedding(graph, embedding, augmented_embedding);
	
	// Create property map to hold the coloring
	typedef iterator_property_map<
			typename std::vector<int>::iterator,
			typename property_map<graph_t, vertex_index_t>::type
		> color_property_map;
	std::vector<int> color_storage(num_vertices(graph));
	color_property_map coloring(color_storage.begin(), get(vertex_index, graph));
	
	// Set up clockwise outer face using properties of a canonical ordering
	std::vector<vertex_t> outer_face = { ordering[1], ordering[0], ordering.back() };
	
	// Create property map to hold the coloring
	typedef iterator_property_map<
			typename std::vector<std::list<int> >::iterator,
			typename property_map<graph_t, vertex_index_t>::type
		> color_list_property_map;
	std::vector<std::list<int> > color_list_storage(num_vertices(graph));
	color_list_property_map color_list(color_list_storage.begin(), get(vertex_index, graph));
	
	std::mt19937 generator;
	std::uniform_int_distribution<int> distribution(0, num_colors - 1);
	
	// Iterate over each vertex
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
	
	// Call path 3-choose algorithm
	#ifdef SHOW_TIMINGS
		auto start = nanosecond_timer::now();
		
		for(std::size_t i = 0; i < 1000; ++i) {
			hartman_skrekovski_choose(graph, augmented_embedding, color_list,
				coloring, outer_face.begin(), outer_face.end());
		}
	
		auto end = nanosecond_timer::now();
		std::cout << "Time = " << (end - start).count() / 1000 << "ns\n";
	#else
		hartman_skrekovski_choose(graph, augmented_embedding, color_list,
				coloring, outer_face.begin(), outer_face.end());
	#endif
	
	#ifdef SHOW_VISUALIZATION
		draw_graph_color(graph, coloring);
	#endif
	
	// Test correctness of path coloring
	test_path_coloring(graph, coloring);
}


void test_augmenting_embeddings() {
	std::cout << "Augmenting Planar Embeddings" << std::endl;
	
	// Define graph properties
	typedef adjacency_list
		<	setS,
			vecS,
			undirectedS,
			property<vertex_index_t, std::size_t>,
			property<edge_index_t, std::size_t>
		> index_graph;
	
	typedef erdos_renyi_iterator<minstd_rand, index_graph> ERGen;
	
	typedef typename graph_traits<index_graph>::edge_descriptor edge_t;
	
	typedef std::vector<
			std::vector<edge_t>
		> embedding_storage_t;

	typedef iterator_property_map
		< typename embedding_storage_t::iterator, 
			typename property_map<index_graph, vertex_index_t>::type
		> embedding_t;
	
	boost::minstd_rand gen;
	gen.seed(8573);
	
	for(std::size_t order = 4; order <= 100; ++order) {
		bool found_planar = false;
		std::size_t count = 4;
		
		while(!found_planar) {
			try {
				// Construct a random trriangulated graph
				index_graph graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
		
				++count;
				
				make_triangulated(graph);
		
				found_planar = true;
				
				// Create the planar embedding
				embedding_storage_t embedding_storage(num_vertices(graph));
				embedding_t embedding(embedding_storage.begin(), get(vertex_index, graph));
	
				boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
					boyer_myrvold_params::embedding = embedding);
		
				//draw_graph_no_color(graph);
		
				test_augmented_embedding_construction(graph, embedding);
		
				#ifdef SHOW_PASSES
					std::cout<<"    PASS " << order << " vertex augmented embedding construction."<<std::endl;
				#endif
			}
			catch(std::logic_error error) {
				// Generated a non-planar graph, ignore this case
			}
			catch(std::exception& error) {
				std::cout<<"    FAIL " << order << "vertex augmented embedding construction ("<<error.what()<<")."<<std::endl;
				failed=true;
			}
			catch(...) {
				std::cout<<"    FAIL " << order << " vertex augmented embedding construction (unknown error)."<<std::endl;
				failed=true;
			}
		}
	}
}

void test_poh_color_bfs() {
	std::cout<<"Path 3-coloring (Poh w/BFS)"<<std::endl;
	
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
	
	for(std::size_t order = 4; order <= 100; ++order) {
		bool found_planar = false;
		std::size_t count = 4;
		
		while(!found_planar) {
			try {
				// Construct a random trriangulated graph
				index_graph graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
		
				++count;
				
				make_triangulated(graph);
		
				found_planar = true;
		
				//draw_graph_no_color(graph);
		
				poh_color_bfs_test(graph);
		
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

void test_poh_color() {
	std::cout<<"Path 3-coloring (Poh w/outer face tracing)"<<std::endl;
	
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
	
	for(std::size_t order = 4; order <= 100; ++order) {
		bool found_planar = false;
		std::size_t count = 4;
		
		while(!found_planar) {
			try {
				// Construct a random trriangulated graph
				index_graph graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
				
				++count;
				
				make_triangulated(graph);
				
				found_planar = true;
				
				//draw_graph_no_color(graph);
				
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


void test_path_choose()
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
	
	for(std::size_t colors = 3; colors < 9; ++colors) {
		std::cout << "Path 3-choosing (Hartman-Skrekovski) with " << colors << " colors" << std::endl;
	
		for(std::size_t order = 4; order <= 100; ++order) {
			bool found_planar = false;
			std::size_t count = 4;
			
			while(!found_planar) {
				try {
					// Construct a random trriangulated graph
					index_graph graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
					
					++count;
					
					make_triangulated(graph);
			
					found_planar = true;
			
					//draw_graph_no_color(graph);
			
					path_choose_test(graph, colors);
			
					#ifdef SHOW_PASSES
						std::cout<<"    PASS " << order << " vertex, " << colors
							<< " colors path 3-choose."<<std::endl;
					#endif
				}
				catch(std::logic_error error) {
					// Generated a non-planar graph, ignore this case
				}
				catch(std::exception& error) {
					std::cout<<"    FAIL " << order << " vertex, " << colors
							<< " colors path 3-choose ("<<error.what()<<")."<<std::endl;
					failed=true;
				}
				catch(...) {
					std::cout<<"    FAIL " << order << " vertex, " << colors
							<< " colors path 3-choose (unknown error)."<<std::endl;
					failed=true;
				}
			}
		}
	}
}

int main() {
	test_augmenting_embeddings();
	test_poh_color_bfs();
	test_poh_color();
	test_path_choose();

	if(failed)
		std::cout<<"THERE ARE FAILING TESTS"<<std::endl;
	else
		std::cout<<"ALL TESTS PASSED"<<std::endl;
	
	return 0;
}
