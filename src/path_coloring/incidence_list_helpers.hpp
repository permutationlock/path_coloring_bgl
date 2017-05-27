/*
 * incidence_list_helpers.h
 * Author: Aven Bross
 * 
 * Simple helper functions for working with incidence lists
 */

#ifndef __INCIDENCE_LIST_HELPERS_HPP
#define __INCIDENCE_LIST_HELPERS_HPP

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>

template<
		typename graph_t, typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename edge_t = typename boost::graph_traits<graph_t>::edge_descriptor
	>
inline vertex_t get_incident_vertex(
		vertex_t vertex, edge_t edge, const graph_t & graph
	)
{
	vertex_t neighbor = boost::source(edge, graph);
	if(neighbor == vertex)
		neighbor = boost::target(edge, graph);
	return neighbor;
}

template<
		typename graph_t, typename planar_embedding_t, typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename edge_iterator_t = typename boost
			::property_traits<planar_embedding_t>::value_type::const_iterator
	>
inline edge_iterator_t find_edge_iterator(
		vertex_t vertex, vertex_t target,
		const planar_embedding_t & planar_embedding, const graph_t & graph
	)
{
	for(edge_iterator_t edge_iter = planar_embedding[vertex].begin();
			edge_iter != planar_embedding[vertex].end(); ++edge_iter
		)
	{
		vertex_t neighbor = get_incident_vertex(vertex, *edge_iter, graph);
		if(neighbor == target) {
			return edge_iter;
		}
	}
	return planar_embedding[vertex].end();
}

template<
		typename graph_t, typename planar_embedding_t, typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename edge_iterator_t = typename boost
			::property_traits<planar_embedding_t>::value_type::const_iterator
	>
inline edge_iterator_t find_edge_iterator_restricted(
		vertex_t vertex, vertex_t target, edge_iterator_t begin,
		edge_iterator_t end, const planar_embedding_t & planar_embedding,
		const graph_t & graph
	)
{
	for(edge_iterator_t edge_iter = begin; edge_iter != end; ++edge_iter) {
		if(edge_iter == planar_embedding[vertex].end())
			edge_iter = planar_embedding[vertex].begin();
		
		vertex_t neighbor = get_incident_vertex(vertex, *edge_iter, graph);
		
		if(neighbor == target) {
			return edge_iter;
		}
	}
	return planar_embedding[vertex].end();
}

template<
		typename graph_t, typename augmented_embedding_t, typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename neighbor_iterator_t = typename boost
			::property_traits<augmented_embedding_t>::value_type::const_iterator
	>
inline neighbor_iterator_t find_neighbor_iterator(
		vertex_t v, vertex_t u,
		const augmented_embedding_t & augmented_embedding, const graph_t & graph
	)
{
	for(auto neighbor_iter = augmented_embedding[v].begin();
			neighbor_iter != augmented_embedding[v].end(); ++neighbor_iter
		)
	{
		vertex_t n = neighbor_iter -> vertex;
		if(n == u) {
			return neighbor_iter;
		}
	}
	return augmented_embedding[v].end();
}

template<
		typename graph_t, typename augmented_embedding_t, typename vertex_t
			= typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename neighbor_iterator_t = typename boost
			::property_traits<augmented_embedding_t>::value_type::const_iterator
	>
inline neighbor_iterator_t find_neighbor_iterator_restricted(
		vertex_t v, vertex_t u, neighbor_iterator_t begin,
		neighbor_iterator_t end,
		const augmented_embedding_t & augmented_embedding, const graph_t & graph
	)
{
	for(neighbor_iterator_t neighbor_iter = begin; neighbor_iter != end;
			++neighbor_iter
		)
	{
		if(neighbor_iter == augmented_embedding[v].end())
			neighbor_iter = augmented_embedding[v].begin();
		
		vertex_t n = neighbor_iter -> vertex;
		if(n == u) {
			return neighbor_iter;
		}
	}
	return augmented_embedding[v].end();
}

template<
		typename vertex_t, typename neighbor_range_map_t,
		typename embedding_t
	>
inline void remove_first_neighbor(
		vertex_t v, neighbor_range_map_t & neighbor_range_map,
		const embedding_t & embedding
	)
{
	if(++neighbor_range_map[v].first == embedding[v].end()) {
		neighbor_range_map[v].first = embedding[v].begin();
	}
}

template<
		typename vertex_t, typename neighbor_range_map_t,
		typename embedding_t
	>
inline void remove_last_neighbor(
		vertex_t v, neighbor_range_map_t & neighbor_range_map,
		const embedding_t & embedding
	)
{
	if(neighbor_range_map[v].second == embedding[v].begin()) {
		neighbor_range_map[v].second = embedding[v].end();
	}
	--neighbor_range_map[v].second;
}

template<
		typename vertex_t, typename embedding_t,
		typename neighbor_iterator_t
			= typename boost::property_traits<embedding_t>::value_type
		    	::const_iterator,
		typename neighbor_range_map_t, typename neighbor_range_t
		    = typename boost::property_traits<neighbor_range_map_t>::value_type
	>
inline std::pair<neighbor_range_t, neighbor_range_t> split_neighbor_range(
		vertex_t v, neighbor_iterator_t mid_iter,
		neighbor_range_map_t & neighbor_range_map,
		embedding_t & embedding
	)
{
	neighbor_range_t range_0, range_1;
	range_0.first = neighbor_range_map[v].first;
	range_0.second = mid_iter;
	
	if(++mid_iter == embedding[v].end())
		mid_iter = embedding[v].begin();
		
	range_1.first = mid_iter;
	range_1.second = neighbor_range_map[v].second;
	
	return std::pair<neighbor_range_t, neighbor_range_t>(range_0, range_1);
}

template<
		typename vertex_t, typename neighbor_range_map_t,
		typename embedding_t, typename neighbor_iterator_t
		    = typename boost::property_traits<embedding_t>::value_type
		    	::const_iterator
	>
inline void initialize_neighbor_range(
		vertex_t v, neighbor_iterator_t start_iter,
		neighbor_range_map_t & neighbor_range_map,
		const embedding_t & embedding
	)
{
	neighbor_range_map[v].first = start_iter;
	neighbor_range_map[v].second = start_iter;
	remove_last_neighbor(v, neighbor_range_map, embedding);
}

#endif
