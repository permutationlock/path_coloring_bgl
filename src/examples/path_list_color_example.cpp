/*
 * path_list_color_example.cpp
 * Author: Aven Bross
 *
 * Example of path list coloring a graph.
 */

// STL headers
#include <iostream>
#include <vector>

// Basic graph headers
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Planar graph headers
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

// Local project headers
#include "../hartman_skrekovski_color/hartman_skrekovski_color.hpp"
#include "../visualization/draw_tikz_graph.hpp"

// Define graph properties
typedef boost::adjacency_list<
		boost::setS,
		boost::vecS,
		boost::undirectedS,
		boost::property<boost::vertex_index_t, std::size_t>,
		boost::no_property
	> index_graph;

// Vertex and edge types
typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
typedef typename boost::graph_traits<index_graph>::edge_descriptor edge_descriptor;

// Define types for the planar embedding
typedef std::vector<std::vector<edge_descriptor> > embedding_storage_t;
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

// Define types for the straight line drawing
typedef std::vector< coord_t > straight_line_drawing_storage_t;
typedef boost::iterator_property_map<
		straight_line_drawing_storage_t::iterator, 
		typename boost::property_map<index_graph, boost::vertex_index_t>::type 
	> straight_line_drawing_t;

int main() {
	/*
	 * First we must construct a weakly plane graph. In boost this means creating two structures:
	 * a graph and an embedding (a property map ordering the edges around each vertex) described
	 * in the boost PlanarEmbedding concept. We must also construct a list of the vertices on the
	 * outer face in clockwise order.
	 */
	
	// Construct a weakly triangulated planar graph somehow (here we manually add edges)
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
	
	// Create the plane embedding
	embedding_storage_t embedding_storage(num_vertices(graph));
	embedding_t embedding(embedding_storage.begin(), boost::get(boost::vertex_index, graph));
	boost::boyer_myrvold_planarity_test(
			boost::boyer_myrvold_params::graph = graph,
			boost::boyer_myrvold_params::embedding = embedding
		);
	
	// We construct a canonical ordering to find the outer face (triangle)
	std::vector<vertex_descriptor> ordering;
	boost::planar_canonical_ordering(graph, embedding, std::back_inserter(ordering));
	
	// Set up clockwise outer face using properties of a canonical ordering
	std::vector<vertex_descriptor> outer_face = { ordering[1], ordering[0], ordering.back() };
	
	
	/*
	 * To call the algoithm we first set up the structures to store the coloring and color lists.
	 * Here each vertex recieves the same 3 color list, but any assignment providing a list of 2 or
	 * more colors to vertices on the outer face, and 3 or more to interior vertices, will suffice.
	 */
	
	// Define a map to write colors to
	std::vector<int> color_storage(boost::num_vertices(graph));
	color_map_t coloring(color_storage.begin(), boost::get(boost::vertex_index, graph));
	
	// Define our set of color lists
	std::vector<std::list<int> > color_list_storage(boost::num_vertices(graph));
	color_list_map_t color_list(color_list_storage.begin(), boost::get(boost::vertex_index, graph));
	
	// Assign a list of 3 colors to each vertex
	typename boost::graph_traits<index_graph>::vertex_iterator v_iter, v_end;
	for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; ++v_iter) {
		for(int color = 0; color < 3; ++color) {
			color_list[*v_iter].push_back(color);
		}
	}
	
	// Path 3-list-color the graph with the Hartman-Skrekovski algorithm
	hartman_skrekovski_color(
			graph,
			embedding,
			color_list,
			coloring,
			outer_face.begin(),
			outer_face.end()
		);
	
	
	/*
	 * Extra steps to produce a tikz drawing of the plane graph and coloring.
	 */
	
	//Set up a property map to hold the mapping from vertices to coord_t's
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
	
	// Print the tikz drawing of the graph
	std::cout << draw_tikz_graph(graph, coloring, straight_line_drawing) << "\n";
	
	return 0;
}