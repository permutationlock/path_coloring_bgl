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
inline edge_iterator_t find_neighbor_iterator(
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
inline edge_iterator_t find_neighbor_iterator_restricted(
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
		typename vertex_t, typename neighbor_range_map_t,
		typename planar_embedding_t
	>
inline void remove_first_edge(
		vertex_t v, neighbor_range_map_t & neighbor_range_map,
		const planar_embedding_t & planar_embedding
	)
{
	if(++neighbor_range_map[v].first == planar_embedding[v].end()) {
		neighbor_range_map[v].first = planar_embedding[v].begin();
	}
}

template<
		typename vertex_t, typename neighbor_range_map_t,
		typename planar_embedding_t
	>
inline void remove_last_edge(
		vertex_t v, neighbor_range_map_t & neighbor_range_map,
		const planar_embedding_t & planar_embedding
	)
{
	if(neighbor_range_map[v].second == planar_embedding[v].begin()) {
		neighbor_range_map[v].second = planar_embedding[v].end();
	}
	--neighbor_range_map[v].second;
}

template<
		typename vertex_t, typename edge_iterator_t,
		typename neighbor_range_map_t, typename planar_embedding_t,
		typename neighbor_range_t
			= typename std::pair<edge_iterator_t, edge_iterator_t>
	>
inline std::pair<neighbor_range_t, neighbor_range_t> split_range(
		vertex_t v, edge_iterator_t mid_iter,
		neighbor_range_map_t & neighbor_range_map,
		planar_embedding_t & planar_embedding
	)
{
	neighbor_range_t first(neighbor_range_map[v].first, mid_iter);
	
	if(++mid_iter == planar_embedding[v].end())
		mid_iter = planar_embedding[v].begin();
	neighbor_range_t second(mid_iter, neighbor_range_map[v].second);
	
	return std::pair<neighbor_range_t, neighbor_range_t>(first, second);
}

template<
		typename vertex_t, typename neighbor_range_map_t,
		typename planar_embedding_t, typename edge_iterator_t
		    = typename boost::property_traits<planar_embedding_t>::value_type
		    	::const_iterator
	>
static inline void initialize(
		vertex_t v, edge_iterator_t start_iter,
		neighbor_range_map_t & neighbor_range_map,
		const planar_embedding_t & planar_embedding
	)
{
	neighbor_range_map[v].first = start_iter;
	neighbor_range_map[v].second = start_iter;
	remove_last_edge(v, neighbor_range_map, planar_embedding);
}

#endif
