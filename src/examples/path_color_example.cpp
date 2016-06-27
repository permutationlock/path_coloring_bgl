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
#include "../hartman_skrekovski_color/hartman_skrekovski_color.hpp"
#include "../visualization/draw_tikz_graph.hpp"

// Define graph properties
typedef boost::adjacency_list<
		boost::setS,
		boost::vecS,
		boost::undirectedS,
		boost::property<boost::vertex_index_t, std::size_t>,
		boost::property<boost::edge_index_t, std::size_t>
	> index_graph;

// Vertex and edge types
typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
typedef typename boost::graph_traits<index_graph>::edge_descriptor edge_descriptor;

// Define types for the planar embedding
typedef std::vector<
		std::vector<edge_descriptor>
	> embedding_storage_t;
typedef boost::iterator_property_map<
		typename embedding_storage_t::iterator, 
		typename boost::property_map<index_graph, boost::vertex_index_t>::type
	> embedding_t;
	
// Define types for the coloring
typedef boost::iterator_property_map<
		std::vector<int>::iterator,
		typename boost::property_map<index_graph, boost::vertex_index_t>::type
	> color_map_t;
	
// Define types for the color list map
typedef boost::iterator_property_map<
		std::vector<std::list<int> >::iterator,
		typename boost::property_map<index_graph, boost::vertex_index_t>::type
	> color_list_map_t;

// A class to hold the coordinates of the straight line embedding
struct coord_t {
  std::size_t x;
  std::size_t y;
};

int main() {
	// Construct a planar graph somehow (here we manually add edges)
	index_graph graph(5);
	boost::add_edge(0, 1, graph);
	boost::add_edge(1, 2, graph);
	boost::add_edge(2, 0, graph);
	boost::add_edge(1, 3, graph);
	boost::add_edge(0, 3, graph);
	boost::add_edge(2, 3, graph);
	boost::add_edge(0, 4, graph);
	boost::add_edge(2, 4, graph);
	boost::add_edge(3, 4, graph);
	
	// Create the planar embedding
	embedding_storage_t embedding_storage(num_vertices(graph));
	embedding_t embedding(embedding_storage.begin(), boost::get(boost::vertex_index, graph));
	boost::boyer_myrvold_planarity_test(
			boost::boyer_myrvold_params::graph = graph,
			boost::boyer_myrvold_params::embedding = embedding
		);
	
	// Construct a canonical ordering to find the outer triangle
	std::vector<vertex_descriptor> ordering;
	boost::planar_canonical_ordering(graph, embedding, std::back_inserter(ordering));
	
	// Set up outer face using properties of a canonical ordering
	std::vector<vertex_descriptor> outer_face = { ordering[1], ordering[0], ordering.back() };
	
	// Define a map to write colors to
	std::vector<int> color_storage(boost::num_vertices(graph));
	color_map_t coloring(color_storage.begin(), boost::get(boost::vertex_index, graph));
	
	// Define our set of color lists
	std::vector<std::list<int> > color_list_storage(boost::num_vertices(graph));
	color_list_map_t color_list(color_list_storage.begin(), boost::get(boost::vertex_index, graph));
	
	// Assign a list of 3 colors to each vertexsssss
	typename boost::graph_traits<index_graph>::vertex_iterator v_iter, v_end;
	for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; ++v_iter) {
		for(int color = 0; color < 3; ++color) {
			color_list[*v_iter].push_back(color);
		}
	}
	
	// Manually set up edge indices (necessary for the way I compute the augmented embedding)
	auto e_index = boost::get(boost::edge_index, graph);
	typename boost::graph_traits<index_graph>::edges_size_type edge_count = 0;
	typename boost::graph_traits<index_graph>::edge_iterator ei, ei_end;
	for(boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei)
		put(e_index, *ei, edge_count++);
	
	// Path 3-list-color the graph with the Hartman-Skrekovski algorithm
	hartman_path_list_color(
			graph,
			embedding,
			color_list,
			coloring,
			outer_face.begin(),
			outer_face.end()
		);
	
	//Set up a property map to hold the mapping from vertices to coord_t's
	typedef std::vector< coord_t > straight_line_drawing_storage_t;
	typedef boost::iterator_property_map
		< straight_line_drawing_storage_t::iterator, 
			typename boost::property_map<index_graph, boost::vertex_index_t>::type 
		> straight_line_drawing_t;
	
	straight_line_drawing_storage_t straight_line_drawing_storage
		(boost::num_vertices(graph));
	straight_line_drawing_t straight_line_drawing
		(straight_line_drawing_storage.begin(), 
			boost::get(boost::vertex_index, graph)
		);
	
	// Draw the plane graph to get integer plane coordinates
	boost::chrobak_payne_straight_line_drawing(graph, 
			embedding, 
			ordering.begin(),
			ordering.end(),
			straight_line_drawing
		);
	
	// Pint the tikz drawing of the graph
	std::cout << draw_tikz_graph(graph, coloring, straight_line_drawing) << "\n";
	
	return 0;
}