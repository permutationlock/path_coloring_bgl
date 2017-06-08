/*
 * hartman_skrekovski_choose.hpp
 * Author: Aven Bross
 * 
 * Implementation of Hartman-Skrekovski path 3-choosing algorithm.
 */

#ifndef __HARTMAN_SKREKOVSKI_CHOOSE_HPP
#define __HARTMAN_SKREKOVSKI_CHOOSE_HPP

// STL headers
#include <vector>
#include <utility>
#include <algorithm>

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Local project headers
#include "disjoint_set.hpp"
#include "incidence_list_helpers.hpp"

namespace {
	// The three states a vertex may have
	const int INTERIOR = 0, ON_FACE = 1, COLORED = 2;

	/* 
	 * set_face_location
	 *
	 * inputs: a vertex v, the new face location for v, a read-write-able vertex
	 *     property map storing the face location of each vertex, a disjoint set
	 *     structure to store face location sets.
	 * 
	 * outputs: sets the face location of v to the given face location, if the
	 * given location doesn't exist, that is, the integer new_face_location is
	 * negative, a new face location is created in the disjoint set structure.
	 */
	template<
			typename vertex_t, typename face_location_map_t
		>
	inline int set_face_location(
			vertex_t v, int new_face_location,
			face_location_map_t & face_location_map,
			disjoint_set_t & face_location_sets
		)
	{
		if(!face_location_sets.exists(new_face_location)) {
			face_location_map[v] = face_location_sets.make_next();
		}
		else {
			face_location_map[v] = face_location_sets.find(new_face_location);
		}
		return face_location_map[v];
	}
	
	/*
	 * hartman_skrekovski_choose
	 *
	 * assumptions: There is a cycle C=v_0v_1...v_k in the weakly triangulated
	 *     being colored. Suppose x, y, and p are in C such that p is between x
	 *     and y clockwise around C. Suppose all vertices other than x, y, and p
	 *     are uncolored. Suppose vertices x, y, and p have lists of size 1 or
	 *     more, other vertices in C have lists of size 2 or more, and vertices
	 *     in the subgraph bounded by C have lists of size 3 or more.
	 *     Additionally, if p is colored some color c, assume vertices between x
	 *     and p clockwise, including x, don't have the color c in their lists.
	 * 
	 * inputs: A weakly triangulated planar graph with vertex indices
	 *     (predefined boost property), a valid augmented planar embedding of
	 *     the graph (modeling the boost AugmentedEmbedding concept defined in
	 *     documentation), a read-write-able vertex property map assigning each
	 *     vertex an integer to track face location, a disjoint set structure to
	 *     store face location sets, a read-write-able vertex property map
	 *     assigning each vertex an integer to track state, a read-able vertex
	 *     property map assigning a pair of iterators from its augmented
	 *     embedding list to track the range of valid neighbors, a read-write-
	 *     able vertex property map to which the coloring will be assigned, the
	 *     vertices x, y, and p on the outer cycle of the current subgraph, and
	 *     the integers for the current face location marks.
	 * 
	 * outputs: The coloring vertex property will contain a valid path list-
	 *     coloring of the subgraph bounded by C such that x, y, and p each
	 *     receive at most one same color neighbor. Moreover, if x=y=p then
	 *     x,y,p receive no same color neighbors.
	 */

	template<
			typename graph_t, typename augmented_embedding_t,
			typename face_location_map_t, typename state_map_t,
			typename neighbor_range_map_t, typename color_list_map_t,
			typename color_map_t, typename vertex_t
		>
	void hartman_skrekovski_choose_recursive(
			const graph_t & graph,
			const augmented_embedding_t & augmented_embedding,
			face_location_map_t & face_location_map,
			disjoint_set_t & face_location_sets, state_map_t & state_map,
			neighbor_range_map_t & neighbor_range_map,
			color_list_map_t & color_list_map, color_map_t & color_map,
			vertex_t x, vertex_t y, vertex_t p,
			int before_p, int before_y, int before_x
		)
	{
		// If p isn't colored yet we are starting a new path
		if(state_map[p] != COLORED) {
			// Color p the first color in its list
			state_map[p] = COLORED;
			color_map[p] = color_list_map[p].front();
		}
		
		// Track the vertices that will become x and y in the first cycle
		vertex_t new_x = x, new_y = y;
		auto neighbor_iter = neighbor_range_map[p].first;
		
		// Iterate through p's adjacency list
		do {
			// Wrap adjacency list
			if(neighbor_iter == augmented_embedding[p].end())
				neighbor_iter = augmented_embedding[p].begin();
			
			// Grab neighbor n, and the iterator to p in n's adjacency list
			vertex_t n = neighbor_iter -> vertex;
			auto back_iter = neighbor_iter -> iterator;
			int n_location = face_location_map[n];
			
			// The case n is not in C
			if(state_map[n] == INTERIOR) {
				state_map[n] = ON_FACE;
				
				// Note that n is between x and p on the first cycle
				before_p = set_face_location(
						n, before_p, face_location_map, face_location_sets
					);
				
				// Initialize n's neighbor range to start at the edge pn
				initialize_neighbor_range(
						n, back_iter, neighbor_range_map, augmented_embedding
					);
				
				// Remove the edge pn from n's neighbor range
				remove_first_neighbor(
						n, neighbor_range_map, augmented_embedding
					);
				
				color_list_map[n].remove(color_map[p]);
			}
			// The case n is immediately counterclockwise to p in C
			else if(neighbor_iter == neighbor_range_map[p].first) {
				// If p has a single neighbor the current subgraph is K_2
				if(neighbor_iter == neighbor_range_map[p].second) {
					if(n != x && n != y) {
						color_list_map[n].remove(color_map[p]);
					}
					
					if(state_map[n] != COLORED) {
						state_map[n] = COLORED;
						color_map[n] = color_list_map[n].front();
					}
					
					break;
				}
				// If n is y, swap p and y so the path runs clockwise
				else if(n == y) {
					if(state_map[y] != COLORED) {
						state_map[y] = COLORED;
						color_map[y] = color_map[p];
					}
					
					if(x == p) new_x = y;
					
					hartman_skrekovski_choose_recursive(
							graph, augmented_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							new_x, p, y,
							-1, -1, before_y
						);
					
					break;
				}
				else {
					// If we are removing x, setup new_x as n
					if(x == p) {
						new_x = n;
						
						// Note that n is between x and p on the first cycle
						before_p = set_face_location(
								n, before_p, face_location_map,
								face_location_sets
							);
					}
					
					color_list_map[n].remove(color_map[p]);
					
					// Remove the edge pn from n's adjacency list
					remove_last_neighbor(
							n, neighbor_range_map, augmented_embedding
						);
				}
			}
			else {
				// Divide the neighbor range of n at the edge pn
				auto n_ranges = split_neighbor_range(
						n, back_iter, neighbor_range_map, augmented_embedding
					);
				
				// p has been removed from the first cycle so we may remove
				//     those neighbors from p's neighbor range
				neighbor_range_map[p].first = neighbor_iter;
				
				// The case we are removing y
				if(p == y) {
					color_list_map[n].remove(color_map[p]);
					
					// If p=x=y then we are removing a single vertex path
					if(p == x) {
						// Color the subgraph bounded by the first cycle
						neighbor_range_map[n] = n_ranges.second;
						hartman_skrekovski_choose_recursive(
								graph, augmented_embedding,
								face_location_map, face_location_sets,
								state_map, neighbor_range_map, 
								color_list_map, color_map,
								new_x, n, new_x,
								-1, before_p, before_y
							);
						
						// If the edge pn is a chord, we must color the rest
						if(back_iter != n_ranges.first.first) {
							// Color the subgraph bounded by the second cycle
							neighbor_range_map[n] = n_ranges.first;
							hartman_skrekovski_choose_recursive(
									graph, augmented_embedding,
									face_location_map, face_location_sets,
									state_map, neighbor_range_map, 
									color_list_map, color_map,
									n, p, p,
									-1, -1, before_y
								);
						}
						
						break;
					}
					// The case the edge pn is in the cycle, i.e. not a chord
					else if(neighbor_iter == neighbor_range_map[p].second) {
						// Color the subgraph bounded by the first cycle
						neighbor_range_map[n] = n_ranges.second;
						hartman_skrekovski_choose_recursive(
								graph, augmented_embedding,
								face_location_map, face_location_sets,
								state_map, neighbor_range_map, 
								color_list_map, color_map,
								new_x, n, new_x,
								-1, before_p, before_x
							);
						
						break;
					}
					// Otherwise we condition on the location of n as usual
					
					// Remember n will be our new y vertex for the first cycle
					new_y = n;
				}
				
				// The case n is in C[p,y] 
				if(face_location_sets.compare(n_location, before_y) || n == y) {
					vertex_t new_p = n;
					
					// The case n may be colored and added to the path
					if(
						color_list_map[n].end() != std::find(
									color_list_map[n].begin(),
									color_list_map[n].begin(), color_map[p]
								)
							&& state_map[n] != COLORED
						)
					{
						color_map[n] = color_map[p];
						state_map[n] = COLORED;
					}
					// The case n will not be added to the path
					else if(
							state_map[n] != COLORED ||
							color_map[n] != color_map[p]
						)
					{
						// Start new path at the new x on the first cycle
						new_p = new_x;
						
						// Combine the before_p and before_y segments
						before_y = face_location_sets.take_union(
								before_p, before_y
							);
						
						// Note before_p is gone
						before_p = -1;
					}
					
					// Color the subgraph bounded by the first cycle
					neighbor_range_map[n] = n_ranges.second;
					hartman_skrekovski_choose_recursive(
							graph, augmented_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							new_x, new_y, new_p,
							before_p, before_y, before_x
						);
					
					// If the edge pn is a chord, handle the rest
					if(back_iter != n_ranges.first.first) {
						// Color the subgraph bounded by the second cycle
						neighbor_range_map[n] = n_ranges.first;
						hartman_skrekovski_choose_recursive(
								graph, augmented_embedding,
								face_location_map, face_location_sets,
								state_map, neighbor_range_map, 
								color_list_map, color_map,
								p, n, p,
								-1, before_y, -1
							);
					}
				}
				// The case n is in C[y,x]
				else if(face_location_sets.compare(n_location, before_x)) {
					color_list_map[n].remove(color_map[p]);
					
					// Color the subgraph bounded by the first cycle
					neighbor_range_map[n] = n_ranges.second;
					hartman_skrekovski_choose_recursive(
							graph, augmented_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							new_x, n, new_x,
							-1, before_p, before_x
						);
					
					// Color the subgraph bounded by the second cycle
					neighbor_range_map[n] = n_ranges.first;
					hartman_skrekovski_choose_recursive(
							graph, augmented_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							n, y, p,
							-1, before_y, before_x
						);
				}
				// The case n is in C[x,p]
				else {
					color_list_map[n].remove(color_map[p]);
					
					// Color the subgraph bounded b the second cycle
					neighbor_range_map[n] = n_ranges.first;
					hartman_skrekovski_choose_recursive(
							graph, augmented_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							x, y, p,
							before_p, before_y, before_x
						);
					
					// Color the subgraph bounded by the first cycle
					neighbor_range_map[n] = n_ranges.second;
					hartman_skrekovski_choose_recursive(
							graph, augmented_embedding,
							face_location_map, face_location_sets,
							state_map, neighbor_range_map, 
							color_list_map, color_map,
							n, n, n,
							-1, before_p, -1
						);
				}
				
				break;
			}
		} while(neighbor_iter++ != neighbor_range_map[p].second);
	}
}


/*
 * hartman_skrekovski_choose
 * 
 * inputs: A weakly triangulated planar graph with vertex indices (predefined
 *     boost property), a valid augmented planar embedding of the graph
 *     (modeling the boost AugmentedEmbedding concept defined in documentation),
 *     a read-able vertex property map assigning a range of colors to each
 *     vertex (each vertex must recieve a range of at least 3 colors if
 *     interior, and at least two colors if on the outer face), a read-write-
 *     able vertex property map to which the coloring will be assigned, and
 *     finally a pair of iterators for the list of vertices on the outer face of
 *     the graph in clockwise order.
 * 
 * outputs: The coloring will be a valid assignment of colors from the input
 *     color lists such that each color class induces a disjoint union of paths.
 */
 
template<
		typename graph_t, typename augmented_embedding_t,
		typename color_list_map_t, typename color_map_t,
		typename face_iterator_t
	>
void hartman_skrekovski_choose(
		const graph_t & graph,
		const augmented_embedding_t & augmented_embedding,
		const color_list_map_t & color_list_map,
		color_map_t & color_map,
		face_iterator_t face_begin,
		face_iterator_t face_end
	)
{
	// Type definitions
	typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
	typedef typename color_map_t::value_type color_t;
	typedef typename boost::graph_traits<graph_t>::vertex_iterator
		vertex_iterator_t;
	typedef typename boost::property_traits<augmented_embedding_t>::value_type
		::const_iterator neighbor_iterator_t;
	
	// Construct a vertex property for start/stop positions in adjacency lists
	std::vector<std::pair<neighbor_iterator_t, neighbor_iterator_t>>
		neighbor_range_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			typename std::vector<std::pair<neighbor_iterator_t,
				neighbor_iterator_t>>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> neighbor_range_map(
				neighbor_range_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
	
	// Construct a vertex property for the location of vertices on outer face
	std::vector<int> face_location_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			std::vector<int>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> face_location_map(
				face_location_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
	
	// Constuct a vertex property for the state of each vertex
	std::vector<int> state_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			std::vector<int>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> state_map(
				state_storage.begin(), boost::get(boost::vertex_index, graph)
			);
	
	// Construct a vertex property to store color lists for each vertex
	std::vector<std::list<color_t>> color_list_storage(
			boost::num_vertices(graph)
		);
	boost::iterator_property_map<
			typename std::vector<std::list<color_t>>::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> color_list_map_copy(
				color_list_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
	
	// Copy color lists
	vertex_iterator_t vertex_iter, vertex_end;
	for(
			boost::tie(vertex_iter, vertex_end) = vertices(graph);
			vertex_iter != vertex_end; vertex_iter++
		)
	{
		vertex_t v = *vertex_iter;
		
		std::copy(
			color_list_map[v].begin(), color_list_map[v].end(),
			std::back_inserter(color_list_map_copy[v])
		);
	}
	
	// Setup face location sets (we will have only the region before_y)
	disjoint_set_t face_location_sets;
	int before_y = -1;
	
	// Initialize vertices on outer face
	for(auto face_iter = face_begin; face_iter != face_end; ++face_iter)
	{
		auto next = face_iter;
		if(++next == face_end) next = face_begin;
		
		vertex_t l = *face_iter;
		vertex_t v = *next;
		
		state_map[v] = ON_FACE;
		
		before_y = set_face_location(
				v, before_y, face_location_map, face_location_sets
			);
		
		auto back_iter = find_neighbor_iterator(
				v, l, augmented_embedding, graph
			);
		initialize_neighbor_range(
				v, back_iter, neighbor_range_map, augmented_embedding
			);
	}
	
	vertex_t x = *face_begin, y = *(--face_end);
	
	// Construct the path choosing from given color lists
	hartman_skrekovski_choose_recursive(
			graph, augmented_embedding,
			face_location_map, face_location_sets,
			state_map, neighbor_range_map, 
			color_list_map_copy, color_map,
			x, y, x,
			-1, before_y, -1
		);
}

#endif
