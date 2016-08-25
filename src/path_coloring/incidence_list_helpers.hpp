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
		typename graph_t,
		typename vertex_t = typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename edge_t = typename boost::graph_traits<graph_t>::edge_descriptor
	>
inline vertex_t get_incident_vertex(vertex_t vertex, edge_t edge, const graph_t & graph) {
	vertex_t neighbor = boost::source(edge, graph);
	if(neighbor == vertex)
		neighbor = boost::target(edge, graph);
	return neighbor;
}

template<
		typename graph_t, typename planar_embedding_t,
		typename vertex_t = typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename edge_iterator_t = typename boost::property_traits<planar_embedding_t>::value_type::const_iterator
	>
inline edge_iterator_t find_neighbor_iterator(vertex_t vertex, vertex_t target,
	const planar_embedding_t & planar_embedding, const graph_t & graph)
{
	for(edge_iterator_t edge_iter = planar_embedding[vertex].begin();
		edge_iter != planar_embedding[vertex].end(); ++edge_iter)
	{
		vertex_t neighbor = get_incident_vertex(vertex, *edge_iter, graph);
		if(neighbor == target) {
			return edge_iter;
		}
	}
	return planar_embedding[vertex].end();
}

template<
		typename graph_t, typename planar_embedding_t,
		typename vertex_t = typename boost::graph_traits<graph_t>::vertex_descriptor,
		typename edge_iterator_t = typename boost::property_traits<planar_embedding_t>::value_type::const_iterator
	>
inline edge_iterator_t find_neighbor_iterator_restricted(vertex_t vertex, vertex_t target,
	edge_iterator_t begin, edge_iterator_t end, const planar_embedding_t & planar_embedding, const graph_t & graph)
{
	for(edge_iterator_t edge_iter = begin; edge_iter != end; ++edge_iter) {
		if(edge_iter == planar_embedding[vertex].end()) edge_iter = planar_embedding[vertex].begin();
		vertex_t neighbor = get_incident_vertex(vertex, *edge_iter, graph);
		if(neighbor == target) {
			return edge_iter;
		}
	}
	return planar_embedding[vertex].end();
}

#endif
