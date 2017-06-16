/*
 * hartman_skrekovski_example.cpp
 * Author: Aven Bross
 *
 * An example using Hartman-Skrekovski to path list-color a plane graph.
 */

// STL headers
#include <algorithm>
#include <iostream>
#include <list>
#include <vector>
#include <utility>

// Basic graph headers
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Planar graph headers
#include <boost/graph/boyer_myrvold_planar_test.hpp>

// Local project headers
#include "../../path_coloring/hartman_skrekovski_color.hpp"
#include "../../path_coloring/augmented_embedding.hpp"
#include "../../path_coloring/incidence_list_helpers.hpp"

using namespace boost;

/* 
 * -----------------------------------------------------------------------------
 *                         Template type definitions
 * -----------------------------------------------------------------------------
 */

// Define the graph type for all test graphs
typedef adjacency_list<
        setS,
        vecS,
        undirectedS,
        property<vertex_index_t, std::size_t>,
        property<edge_index_t, std::size_t>
    > graph_t;

// Vertex and edge types
typedef typename graph_traits<graph_t>::vertex_descriptor vertex_t;
typedef typename graph_traits<graph_t>::edge_descriptor edge_t;

// Vertex and edge index map types
typedef typename property_map<graph_t, vertex_index_t>::const_type
    vertex_index_map_t;

// Vertex property map type for an integer property
typedef std::vector<int> integer_property_storage_t;
typedef iterator_property_map<
        typename integer_property_storage_t::iterator, 
        vertex_index_map_t
    > integer_property_map_t;

// Vertex property map type for a planar embedding
typedef std::vector<std::vector<edge_t>> planar_embedding_storage_t;
typedef iterator_property_map<
        typename planar_embedding_storage_t::iterator, 
        vertex_index_map_t
    > planar_embedding_t;

// Vertex property map to store color lists
typedef std::vector<std::list<int>> color_list_storage_t;
typedef iterator_property_map<
        typename color_list_storage_t::iterator,
        vertex_index_map_t
    > color_list_map_t;

// Vertex property for an augmented embedding
struct adjacency_node_t {
        vertex_t vertex;
        typename std::vector<adjacency_node_t>::iterator iterator;
    };
typedef typename std::vector<std::vector<adjacency_node_t>>
    augmented_embedding_storage_t;
typedef iterator_property_map<
        typename augmented_embedding_storage_t::iterator,
        vertex_index_map_t
    > augmented_embedding_t;

// Vertex property map type for the neighbor ranges of planar_embedding_t
typedef typename property_traits<augmented_embedding_t>::value_type
        ::const_iterator embedding_iterator_t;
typedef typename std::vector<
        std::pair<embedding_iterator_t, embedding_iterator_t>
    > neighbor_range_storage_t;
typedef iterator_property_map<
        typename neighbor_range_storage_t::iterator,
        vertex_index_map_t
    > neighbor_range_map_t;


/* 
 * -----------------------------------------------------------------------------
 *               Main: path 3-color a plane graph
 * -----------------------------------------------------------------------------
 */

int main() {
    graph_t graph(5);
    
    // First we will construct a planar graph.
    boost::add_edge(0, 1, graph);
	boost::add_edge(1, 2, graph);
	boost::add_edge(2, 0, graph);
	boost::add_edge(1, 3, graph);
	boost::add_edge(0, 3, graph);
	boost::add_edge(2, 3, graph);
	boost::add_edge(0, 4, graph);
	boost::add_edge(2, 4, graph);
	boost::add_edge(3, 4, graph);
	
	// We choose our outer face to be the triangle 012
	std::vector<vertex_t> cycle = { 0, 1, 2 };
    
	// Create the planar embedding
	planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    for(std::size_t v = 0; v < num_vertices(graph); ++v) {
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );
    
    // Create the augmented embedding
    augmented_embedding_storage_t augmented_embedding_storage(
            num_vertices(graph)
        );
    augmented_embedding_t augmented_embedding(
            augmented_embedding_storage.begin(), get(vertex_index, graph)
        );
    for(std::size_t v = 0; v < num_vertices(graph); ++v) {
        augmented_embedding[v].reserve(out_degree(v, graph));
    }
    augment_embedding(graph, planar_embedding, augmented_embedding);
    
    // Print embedding ordered adjacency list
    std::cout << "Embedding ordered adajacency lists:\n";
    for(std::size_t v = 0; v < num_vertices(graph); ++v) {
        std::cout << "    Adj[" << v << "] = ";
        for(auto node_iter = augmented_embedding[v].begin();
            node_iter != augmented_embedding[v].end(); ++node_iter)
        {
            if(node_iter != augmented_embedding[v].begin())
                std::cout << " -> ";
            std::cout << node_iter -> vertex;
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Create a vertex property map to store color lists
	color_list_storage_t color_list_storage(num_vertices(graph));
    color_list_map_t color_list_map(
            color_list_storage.begin(), get(vertex_index, graph)
        );
    
    // Give each vertex a list of the appropriate number of colors
    color_list_map[0] = { 1, 2 };
    color_list_map[1] = { 2, 3 };
    color_list_map[2] = { 1, 4 };
    color_list_map[3] = { 1, 3, 4 };
    color_list_map[4] = { 1, 2, 4 };
    
    // Print color lists
    std::cout << "Color list assignment:\n";
    for(std::size_t v = 0; v < num_vertices(graph); ++v) {
        std::cout << "    L[" << v << "] = ";
        for(auto color_iter = color_list_map[v].begin();
            color_iter != color_list_map[v].end(); ++color_iter)
        {
            if(color_iter != color_list_map[v].begin())
                std::cout << " -> ";
            std::cout << *color_iter;
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Construct a vertex property for neighbor ranges
    neighbor_range_storage_t neighbor_range_storage(num_vertices(graph));
    neighbor_range_map_t neighbor_range_map(
            neighbor_range_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property map to store face location marks
    integer_property_storage_t face_location_storage(num_vertices(graph));
    integer_property_map_t face_location_map(
            face_location_storage.begin(), get(vertex_index, graph)
        );
    
    
    // Call Poh with the given cycle and color set { 1, 2, 3 }
    hartman_skrekovski_color(
            graph, augmented_embedding, color_list_map,
            neighbor_range_map, face_location_map,
            cycle.begin(), cycle.end()
        );
    
    // Print the coloring
    std::cout << "The path L-list-coloring:\n";
    for(std::size_t v = 0; v < num_vertices(graph); ++v) {
        std::cout << "    color[" << v << "] = "
            << color_list_map[v].front() << "\n";
    }
    
    return 0;
}

