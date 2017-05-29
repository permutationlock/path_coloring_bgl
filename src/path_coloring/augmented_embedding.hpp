/*
 * disjoint_set.hpp
 * Author: Aven Bross
 *
 * Produces an augmented embedding structure from a planar embedding
 */

#ifndef __AUGMENTED_EMBEDDING_HPP
#define __AUGMENTED_EMBEDDING_HPP

// STL headers
#include <vector>

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

// Local project headers
#include "incidence_list_helpers.hpp"


template<
		typename graph_t, typename planar_embedding_t,
		typename augmented_embedding_t
	>
void augment_embedding(
		const graph_t & graph, const planar_embedding_t & planar_embedding,
		augmented_embedding_t & augmented_embedding
	)
{
	// Type definitions
	typedef typename boost::property_traits<augmented_embedding_t>::value_type
		::value_type adjacency_node_t;
	typedef typename boost::property_traits<augmented_embedding_t>::value_type
		::iterator adjacency_node_iterator_t;
	typedef typename boost::graph_traits<graph_t>::vertex_iterator
		vertex_iterator_t;
	typedef typename boost::graph_traits<graph_t>::vertex_descriptor
		vertex_t;
	
	vertex_iterator_t v_begin, v_iter, v_end;
	
	// Create the embedding ordered adjacency list with uninitialized iterators
	for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end;
			v_iter++
		)
	{
		vertex_t v = *v_iter;
		
		for(auto edge_iter = planar_embedding[v].begin();
				edge_iter != planar_embedding[v].end(); ++edge_iter
			)
		{
			vertex_t u = get_incident_vertex(v, *edge_iter, graph);
			adjacency_node_t node;
			node.vertex = u;
			augmented_embedding[v].push_back(node);
		}
	}
	
	// Construct worklist data structure
	std::vector<std::vector<adjacency_node_t>>
		work_list_storage(boost::num_vertices(graph));
	boost::iterator_property_map<
			typename std::vector<std::vector<adjacency_node_t>>
				::iterator,
			typename boost::property_map<graph_t, boost::vertex_index_t>
				::const_type
		> work_list_map(
				work_list_storage.begin(),
				boost::get(boost::vertex_index, graph)
			);
	
	// Reserve deg(v) space for each v in G
	for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end;
			v_iter++
		)
	{
		vertex_t v = *v_iter;
		work_list_map[v].reserve(out_degree(v, graph));
	}
	
	// Construct the worklist for each vertex
	for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end;
			v_iter++
		)
	{
		vertex_t v = *v_iter;
		
		for(auto adjacency_node_iter = augmented_embedding[v].begin();
				adjacency_node_iter != augmented_embedding[v].end();
				++adjacency_node_iter
			)
		{
			vertex_t u = adjacency_node_iter -> vertex;
			
			adjacency_node_t node;
			node.vertex = v;
			node.iterator = adjacency_node_iter;
			work_list_map[u].push_back(node);
		}
	}
	
	// Find iterators using worklists
	boost::tie(v_begin, v_iter) = boost::vertices(graph);
	while(v_iter != v_begin) {
		--v_iter;
		
		vertex_t v = *v_iter;
		
		for(auto work_list_iter = work_list_map[v].begin();
				work_list_iter != work_list_map[v].end();
				++work_list_iter
			)
		{
			// Grab the vertex and iterator for this worklist node
			vertex_t u = work_list_iter -> vertex;
			adjacency_node_iterator_t u_v_iterator
				= work_list_iter -> iterator;
			
			// The node for v must be the last node in u's worklist, so grab it
			adjacency_node_iterator_t v_u_iterator
				= work_list_map[u].back().iterator;
			
			// Assign u's iterator to v's node and vice versa
			u_v_iterator -> iterator = v_u_iterator;
			v_u_iterator -> iterator = u_v_iterator;
			
			// Remove the node for v from u's worklist
			work_list_map[u].pop_back();
		}
	}
}


#endif

