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
#include <set>
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
#include "../path_coloring/hartman_skrekovski_color.hpp"
#include "../visualization/draw_tikz_graph.hpp"

using namespace boost;


/* 
 * -----------------------------------------------------------------------------
 *                             Test option macros
 * -----------------------------------------------------------------------------
 */

// Comment line below to hide passing tests.
//#define SHOW_PASSES

// Comment line below to hide tikz drawing printouts
//#define SHOW_VISUALIZATION

// Comment line below to hide color list assignment printouts
//#define SHOW_COLOR_LISTS


bool failed=false;

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

// Vertex iterator type
typedef typename graph_traits<graph_t>::vertex_iterator vertex_iterator_t;
typedef typename graph_traits<graph_t>::edge_iterator edge_iterator_t;

// Vertex and edge index map types
typedef typename property_map<graph_t, vertex_index_t>::const_type
    vertex_index_map_t;
typedef typename property_map<graph_t, edge_index_t>::const_type
    edge_index_map_t;

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

// Vertex property map type for the neighbor ranges of augmented_embedding_t
typedef typename property_traits<augmented_embedding_t>::value_type
        ::const_iterator augmented_embedding_iterator_t;
typedef typename std::vector<
        std::pair<
                augmented_embedding_iterator_t, augmented_embedding_iterator_t
            >
    > augmented_neighbor_range_storage_t;
typedef iterator_property_map<
            typename augmented_neighbor_range_storage_t::iterator,
            vertex_index_map_t
        > augmented_neighbor_range_map_t;

// Vertex property for a graph drawing; gives each vertex cartesians coordinates
struct coord_t {
        std::size_t x;
        std::size_t y;
    };
typedef std::vector<coord_t> straight_line_drawing_storage_t;
typedef iterator_property_map<
        straight_line_drawing_storage_t::iterator, 
        vertex_index_map_t
    > straight_line_drawing_t;

// Random graph generator
typedef erdos_renyi_iterator<minstd_rand, graph_t> ERGen;

/* 
 * -----------------------------------------------------------------------------
 *                         Planar graph triangulation
 * -----------------------------------------------------------------------------
 */

void make_triangulated(graph_t & graph) {
    make_connected(graph);
    
    // Initialize the interior edge index
    auto edge_index_map = get(edge_index, graph);
    typename graph_traits<graph_t>::edges_size_type edge_count = 0;
    edge_iterator_t edge_iter, edge_end;
    for(tie(edge_iter, edge_end) = edges(graph); edge_iter != edge_end;
        ++edge_iter)
    {
        put(edge_index_map, *edge_iter, edge_count++);
    }

    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    if(!boyer_myrvold_planarity_test(boyer_myrvold_params::graph = graph,
        boyer_myrvold_params::embedding = planar_embedding))
    {
        std::string error = "Non-planar graph.";
        throw std::logic_error(error);
    }

    make_biconnected_planar(graph, planar_embedding);
    
    // Re-initialize the edge index, since we just added a few edges
    edge_count = 0;
    for(tie(edge_iter, edge_end) = edges(graph); edge_iter != edge_end;
        ++edge_iter)
    {
        put(edge_index_map, *edge_iter, edge_count++);
    }
    
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );

    make_maximal_planar(graph, planar_embedding);
    
    // Re-initialize the edge index, since we just added a few edges
    edge_count = 0;
    for(tie(edge_iter, edge_end) = edges(graph); edge_iter != edge_end;
        ++edge_iter)
    {
        put(edge_index_map, *edge_iter, edge_count++);
    }
}


/* 
 * -----------------------------------------------------------------------------
 *                            Planar graph drawing
 * -----------------------------------------------------------------------------
 */

void draw_graph_no_color(graph_t & graph) {
    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );
    
    // Create a vertex property map for the coloring
    integer_property_storage_t color_map_storage(num_vertices(graph));
      integer_property_map_t color_map(
              color_map_storage.begin(), get(vertex_index, graph)
        );
        
    // Find a canonical ordering
    std::vector<vertex_t> ordering;
    planar_canonical_ordering(
            graph, planar_embedding, std::back_inserter(ordering)
        );
        
    // Compute the straight line drawing
    straight_line_drawing_storage_t straight_line_drawing_storage(
            num_vertices(graph)
        );
    straight_line_drawing_t straight_line_drawing(
            straight_line_drawing_storage.begin(), 
            get(vertex_index, graph)
        );
    chrobak_payne_straight_line_drawing(graph, 
            planar_embedding, 
            ordering.begin(),
            ordering.end(),
            straight_line_drawing
        );
    
    std::cout << draw_tikz_graph(
            graph, color_map, straight_line_drawing
        ) << "\n";
}

void draw_graph_color(
        const graph_t & graph, const integer_property_map_t & coloring
    )
{
    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );
        
    // Find a canonical ordering
    std::vector<vertex_t> ordering;
    planar_canonical_ordering(
            graph, planar_embedding, std::back_inserter(ordering)
        );
        
    // Compute the straight line drawing
    straight_line_drawing_storage_t straight_line_drawing_storage(
            num_vertices(graph)
        );
    straight_line_drawing_t straight_line_drawing(
            straight_line_drawing_storage.begin(), 
            get(vertex_index, graph)
        );
    chrobak_payne_straight_line_drawing(graph, 
            planar_embedding, 
            ordering.begin(),
            ordering.end(),
            straight_line_drawing
        );
    
    std::cout << draw_tikz_graph(graph, coloring, straight_line_drawing)
        << "\n";
}


/* 
 * -----------------------------------------------------------------------------
 *                        Augmented embedding tests
 * -----------------------------------------------------------------------------
 */

void test_augmented_embedding(
        const graph_t & graph, const planar_embedding_t & planar_embedding,
        const augmented_embedding_t & augmented_embedding
    )
{
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
                std::string error = "vertex " + std::to_string(v)
                    + " has neighbor " + std::to_string(u_0)
                    + " in embedding, but neighbor " + std::to_string(u_1)
                    + " in augmented_embedding.";
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

void test_augmented_embedding_construction(
        const graph_t & graph, const planar_embedding_t & planar_embedding
    )
{
    // Compute an augmented embedding
    std::vector<std::vector<adjacency_node_t>>
        adjacency_node_storage(num_vertices(graph));
    augmented_embedding_t augmented_embedding(
            adjacency_node_storage.begin(),
            get(vertex_index, graph)
        );
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        augmented_embedding[v].reserve(out_degree(v, graph));
    }
    augment_embedding(graph, planar_embedding, augmented_embedding);
    
    // Run the test
    test_augmented_embedding(graph, planar_embedding, augmented_embedding);
}


/* 
 * -----------------------------------------------------------------------------
 *                           Path coloring test
 * -----------------------------------------------------------------------------
 */

void test_path_coloring(
        const graph_t & graph, const integer_property_map_t & coloring
    )
{
    // Type definitions
    typedef typename property_traits<integer_property_map_t>::value_type
        color_t;
    typedef typename graph_traits<graph_t>::adjacency_iterator
        adjacency_iterator_t;
    
    std::set<vertex_t> visited;
    
    // Iterate over each vertex
    vertex_iterator_t v_iter, v_end;
    for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t curr_vertex = *v_iter;
        color_t curr_color = coloring[curr_vertex];
        
        // If we have not yet visited this vertex, start bfs on same color vertices
        if(visited.count(curr_vertex) == 0) {
            std::queue<vertex_t> bfs_queue;
            bfs_queue.push(curr_vertex);
            
            std::size_t extra = 1;
            
            while(!bfs_queue.empty()) {
                vertex_t v = bfs_queue.front();
                bfs_queue.pop();
                
                if(visited.count(v) == 0) {
                    std::size_t new_neighbor_count = 0, old_neighbor_count = 0;
                    
                    adjacency_iterator_t n_iter, n_end;
                    for(tie(n_iter, n_end) = adjacent_vertices(v, graph); n_iter != n_end; n_iter++) {
                        vertex_t n = *n_iter;
                        
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

void test_list_coloring(
        const graph_t & graph, const color_list_map_t & original_color_list_map,
        const color_list_map_t & color_list_map
    )
{
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        if(color_list_map[v].size() < 1) {
            std::string error = "vertex "
                + std::to_string(v) + " has an empty color list.";
            throw std::runtime_error(error);
        }
        
        if(color_list_map[v].size() > 1) {
            std::string error = "vertex "
                + std::to_string(v) + " has more than one color in its list.";
            throw std::runtime_error(error);
        }
        
        auto v_color = color_list_map[v].front();
        auto color_iter = std::find(
                original_color_list_map[v].begin(),
                original_color_list_map[v].end(), v_color
            );
        
        if(color_iter == original_color_list_map[v].end()) {
            std::string error = "vertex "
                + std::to_string(v) + " has the color "
                + std::to_string(v_color) + " not in the original list.";
            throw std::runtime_error(error);
        }
    }
}


/* 
 * -----------------------------------------------------------------------------
 *                     Algorithm tests for single graph
 * -----------------------------------------------------------------------------
 */

void poh_color_bfs_test(const graph_t & graph) {
    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );
    
    // Create a vertex property map for the coloring
    integer_property_storage_t color_map_storage(num_vertices(graph));
      integer_property_map_t color_map(
              color_map_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property for marking vertices
    integer_property_storage_t mark_storage(num_vertices(graph));
    integer_property_map_t mark_map(
                mark_storage.begin(), get(vertex_index, graph)
            );
    
    // Construct a vertex property for storing parents to track the BFS tree
    std::vector<vertex_t> parent_storage(num_vertices(graph));
    iterator_property_map<
            typename std::vector<vertex_t>::iterator,
            vertex_index_map_t
        > parent_map(
                parent_storage.begin(), get(vertex_index, graph)
            );
        
    // Find a canonical ordering
    std::vector<vertex_t> ordering;
    planar_canonical_ordering(
            graph, planar_embedding, std::back_inserter(ordering)
        );
    
    // Create vectors to hold chordless paths
    std::vector<vertex_t> p;
    std::vector<vertex_t> q;
    
    // Initialize path p with top vertex of outer face
    p.push_back(ordering[0]);

    // Initialize path q with bottom vertices of outer face
    q.push_back(ordering[1]);
    q.push_back(ordering.back());
    
    // Call Poh algorithm
    poh_color_bfs(
            graph, planar_embedding, color_map, mark_map,
            parent_map, p.begin(), p.end(), q.begin(), q.end(), 1, 2, 3
        );
    
    #ifdef SHOW_VISUALIZATION
        draw_graph_color(graph, color_map);
    #endif
    
    // Test correctness of path coloring
    test_path_coloring(graph, color_map);
}

void poh_color_test(const graph_t & graph) {
    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        planar_embedding[v].reserve(out_degree(v, graph));
    }
    boyer_myrvold_planarity_test(
            boyer_myrvold_params::graph = graph,
            boyer_myrvold_params::embedding = planar_embedding
        );
    
    // Create a vertex property map for the coloring
    integer_property_storage_t color_map_storage(num_vertices(graph));
    integer_property_map_t color_map(
            color_map_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property for marking vertices
    integer_property_storage_t mark_storage(num_vertices(graph));
    integer_property_map_t mark_map(
            mark_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property for neighbor ranges
    neighbor_range_storage_t neighbor_range_storage(num_vertices(graph));
    neighbor_range_map_t neighbor_range_map(
            neighbor_range_storage.begin(), get(vertex_index, graph)
        );
    
    // Find a canonical ordering
    std::vector<vertex_t> ordering;
    planar_canonical_ordering(
            graph, planar_embedding, std::back_inserter(ordering)
        );
    
    // Create vectors to hold chordless paths
    std::vector<vertex_t> p;
    std::vector<vertex_t> q;
    
    // Initialize path p with top vertex of outer face
    p.push_back(ordering[0]);

    // Initialize path q with bottom vertices of outer face
    q.push_back(ordering[1]);
    q.push_back(ordering.back());
    
    // Call Poh algorithm
    poh_color(
            graph, planar_embedding, color_map, neighbor_range_map,
            mark_map, p.begin(), p.end(), q.begin(), q.end(), 1, 2, 3
        );
    
    #ifdef SHOW_VISUALIZATION
        draw_graph_color(graph, color_map);
    #endif
    
    // Test correctness of path coloring
    test_path_coloring(graph, color_map);
}


void hartman_skrekovski_test(const graph_t & graph, std::size_t num_colors) {
    // Create the planar embedding
    planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
    planar_embedding_t planar_embedding(
            planar_embedding_storage.begin(), get(vertex_index, graph)
        );
    vertex_iterator_t v_iter, v_end;
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
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
    for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        augmented_embedding[v].reserve(out_degree(v, graph));
    }
    augment_embedding(graph, planar_embedding, augmented_embedding);
    
    // Construct a vertex property for marking vertices
    integer_property_storage_t face_location_storage(num_vertices(graph));
    integer_property_map_t face_location_map(
            face_location_storage.begin(), get(vertex_index, graph)
        );
    
    // Construct a vertex property for neighbor ranges
    augmented_neighbor_range_storage_t neighbor_range_storage(
            num_vertices(graph)
        );
    augmented_neighbor_range_map_t neighbor_range_map(
            neighbor_range_storage.begin(), get(vertex_index, graph)
        );
    
    // Find a canonical ordering
    std::vector<vertex_t> ordering;
    planar_canonical_ordering(
            graph, planar_embedding, std::back_inserter(ordering)
        );
    
    // Grab a triangle to use as the outer face from the ordering
    std::vector<vertex_t> outer_face = {
            ordering[1], ordering[0], ordering.back()
        };
    
    // Create property map to hold the list assignment
    color_list_storage_t color_list_storage(num_vertices(graph));
    color_list_map_t color_list_map(
            color_list_storage.begin(), get(vertex_index, graph)
        );
    color_list_storage_t original_color_list_storage(num_vertices(graph));
    color_list_map_t original_color_list_map(
            original_color_list_storage.begin(), get(vertex_index, graph)
        );
    
    // Uniform distribution to generate random colors
    std::mt19937 generator;
    std::uniform_int_distribution<int> distribution(1, num_colors);
    
    // Iterate over each vertex
    for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        std::vector<int> random_colors(3);
        while(
                random_colors[0] == random_colors[1] ||
                random_colors[0] == random_colors[2] ||
                random_colors[1] == random_colors[2]
            )
        {
            for(std::size_t i = 0; i < 3; ++i) {
                random_colors[i] = distribution(generator);
            }
        }
        
        #ifdef SHOW_COLOR_LISTS
            std::cout << "color_list_map[" << *v_iter << "] = { ";
            for(std::size_t i = 0; i < 3; ++i) {
                std::cout << random_colors[i] << ((i != 2) ? ", " : "");
            }
            std::cout << "}\n";
        #endif
        
        std::copy(
                random_colors.begin(), random_colors.end(),
                std::back_inserter(color_list_map[*v_iter])
            );
        std::copy(
                random_colors.begin(), random_colors.end(),
                std::back_inserter(original_color_list_map[*v_iter])
            );
    }
    
    // Call path list-coloring algorithm
    hartman_skrekovski_color(
            graph, augmented_embedding, color_list_map,
            neighbor_range_map, face_location_map,
            outer_face.begin(), outer_face.end()
        );
    
    // Test correctness of list coloring
    test_list_coloring(graph, original_color_list_map, color_list_map);
    
    // Create a vertex property map for the coloring
    integer_property_storage_t color_map_storage(num_vertices(graph));
    integer_property_map_t color_map(
            color_map_storage.begin(), get(vertex_index, graph)
        );
    
    // Grab the color from each vertex's list
    for (tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
        vertex_t v = *v_iter;
        color_map[v] = color_list_map[v].front();
    }
    
    #ifdef SHOW_VISUALIZATION
        draw_graph_color(graph, color_map);
    #endif
    
    // Test correctness of path coloring
    test_path_coloring(graph, color_map);
}


/* 
 * -----------------------------------------------------------------------------
 *             Run each algorithm test on a set of generated graphs
 * -----------------------------------------------------------------------------
 */

void test_augmenting_embeddings() {
    std::cout << "Augmenting Planar Embeddings" << std::endl;
    
    minstd_rand gen;
    gen.seed(8573);
    
    for(std::size_t order = 4; order <= 200; order+=10) {
        bool found_planar = false;
        std::size_t count = 4;
        
        while(!found_planar) {
            try {
                // Construct a random trriangulated graph
                graph_t graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
        
                ++count;
                
                make_triangulated(graph);
        
                found_planar = true;
                
                // Create the planar embedding
                planar_embedding_storage_t planar_embedding_storage(num_vertices(graph));
                planar_embedding_t planar_embedding(
                        planar_embedding_storage.begin(), get(vertex_index, graph)
                    );
                vertex_iterator_t v_iter, v_end;
                for(tie(v_iter, v_end) = vertices(graph); v_iter != v_end; v_iter++) {
                    vertex_t v = *v_iter;
                    planar_embedding[v].reserve(out_degree(v, graph));
                }
                boyer_myrvold_planarity_test(
                        boyer_myrvold_params::graph = graph,
                        boyer_myrvold_params::embedding = planar_embedding
                    );
        
                //draw_graph_no_color(graph);
        
                test_augmented_embedding_construction(graph, planar_embedding);
        
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
    
    minstd_rand gen;
    gen.seed(8573);
    
    for(std::size_t order = 4; order <= 200; order+=10) {
        bool found_planar = false;
        std::size_t count = 4;
        
        while(!found_planar) {
            try {
                // Construct a random trriangulated graph
                graph_t graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
        
                ++count;
                
                make_triangulated(graph);
        
                found_planar = true;
        
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
    
    minstd_rand gen;
    gen.seed(8573);
    
    for(std::size_t order = 4; order <= 200; order+=10) {
        bool found_planar = false;
        std::size_t count = 4;
        
        while(!found_planar) {
            try {
                // Construct a random trriangulated graph
                graph_t graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
                
                ++count;
                
                make_triangulated(graph);
                
                found_planar = true;
                
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


void test_hartman_skrekovski_color()
{
    minstd_rand gen;
    gen.seed(8573);
    
    for(std::size_t colors = 3; colors < 14; colors += 2) {
        std::cout << "Path 3-choosing (Hartman-Skrekovski) with " << colors << " colors" << std::endl;
    
        for(std::size_t order = 4; order <= 200; order+=10) {
            bool found_planar = false;
            std::size_t count = 4;
            
            while(!found_planar) {
                try {
                    // Construct a random trriangulated graph
                    graph_t graph(ERGen(gen, order, 2 * order - count), ERGen(), order);
                    
                    ++count;
                    
                    make_triangulated(graph);
            
                    found_planar = true;
            
                    hartman_skrekovski_test(graph, colors);
            
                    #ifdef SHOW_PASSES
                        std::cout<<"    PASS " << order << " vertex, " << colors
                            << " colors path 3-list-color."<<std::endl;
                    #endif
                }
                catch(std::logic_error error) {
                    // Generated a non-planar graph, ignore this case
                }
                catch(std::exception& error) {
                    std::cout<<"    FAIL " << order << " vertex, " << colors
                            << " colors path 3-list-color ("<<error.what()<<")."<<std::endl;
                    failed=true;
                }
                catch(...) {
                    std::cout<<"    FAIL " << order << " vertex, " << colors
                            << " colors path 3-list-color (unknown error)."<<std::endl;
                    failed=true;
                }
            }
        }
    }
}


/* 
 * -----------------------------------------------------------------------------
 *                             Main: run tests
 * -----------------------------------------------------------------------------
 */

int main() {
    test_augmenting_embeddings();
    test_poh_color_bfs();
    test_poh_color();
    test_hartman_skrekovski_color();

    if(failed)
        std::cout<<"THERE ARE FAILING TESTS"<<std::endl;
    else
        std::cout<<"ALL TESTS PASSED"<<std::endl;
    
    return 0;
}

