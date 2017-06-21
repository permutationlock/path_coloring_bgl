/*
 * poh_example.cpp
 * Author: Aven Bross
 *
 * An example using Poh with face walking to path 3-color a plane graph.
 */

// STL headers
#include <algorithm>
#include <iostream>
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
#include "../../path_coloring/poh_color.hpp"
#include "../../path_coloring/incidence_list_helpers.hpp"

using namespace boost;

/* 
 * -----------------------------------------------------------------------------
 *                         Template type definitions
 * -----------------------------------------------------------------------------
 */

// Define the graph type for all test graphs
typedef adjacency_list<
        vecS,
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

// Vertex property map type for the neighbor ranges of planar_embedding_t
typedef typename property_traits<planar_embedding_t>::value_type
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
    // First we will construct a planar graph.
    graph_t graph(5);
    
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
    std::vector<vertex_t> path_1 = { 0 };
    std::vector<vertex_t> path_2 = { 1, 2 };
    
    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    for(std::size_t v=0; v < num_vertices(graph); ++v) {
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );
    
    // Print embedding ordered adjacency list
    std::cout << "Embedding ordered adajacency lists:\n";
    for(std::size_t v=0; v < num_vertices(graph); ++v) {
        std::cout << "    Adj[" << v << "] = ";
        for(auto edge_iter = planar_embedding[v].begin();
            edge_iter != planar_embedding[v].end(); ++edge_iter)
        {
            if(edge_iter != planar_embedding[v].begin())
                std::cout << " -> ";
            std::cout << get_incident_vertex(v, *edge_iter, graph);
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    
    // Create a vertex property map for the coloring
    integer_property_storage_t color_map_storage(num_vertices(graph));
    integer_property_map_t color_map(
            color_map_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property map to store vertex marks
    integer_property_storage_t mark_storage(num_vertices(graph));
    integer_property_map_t mark_map(
            mark_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property map to store neighbor ranges
    neighbor_range_storage_t neighbor_range_storage(num_vertices(graph));
    neighbor_range_map_t neighbor_range_map(
            neighbor_range_storage.begin(), get(vertex_index, graph)
        );
    
    // Call Poh with the given paths and structurs and color set { 1, 2, 3 }
    poh_color(
            graph, planar_embedding, color_map, neighbor_range_map,
            mark_map, path_1.begin(), path_1.end(),
            path_2.begin(), path_2.end(),
            1, 2, 3
        );
    
    // Print the coloring
    std::cout << "The path 3-coloring:\n";
    for(std::size_t v=0; v < num_vertices(graph); ++v) {
        std::cout << "    color[" << v << "] = " << color_map[v] << "\n";
    }
    
    return 0;
}

