/*
 * incidence_list_helpers.h
 * Author: Aven Bross
 * 
 * Simple helper functions for working with incidence lists
 */

#ifndef __INCIDENCE_LIST_HELPERS_HPP
#define __INCIDENCE_LIST_HELPERS_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

template<
		typename graph,
		typename vertex_descriptor = typename boost::graph_traits<graph>::vertex_descriptor,
		typename edge_descriptor = typename boost::graph_traits<graph>::edge_descriptor
	>
inline vertex_descriptor get_incident_vertex(vertex_descriptor v, edge_descriptor e, const graph & g) {
	vertex_descriptor n = boost::source(e, g);
	if(n == v)
		n = boost::target(e, g);
	return n;
}

template<
		typename graph, typename planar_embedding,
		typename vertex_descriptor = typename boost::graph_traits<graph>::vertex_descriptor,
		typename edge_iterator = typename boost::property_traits<planar_embedding>::value_type::const_iterator
	>
inline edge_iterator find_neighbor_iterator(vertex_descriptor v, vertex_descriptor u,
	const planar_embedding & e, const graph & g)
{
	for(edge_iterator edge_iter = e[v].begin(); edge_iter != e[v].end(); ++edge_iter) {
		vertex_descriptor n = get_incident_vertex(v, *edge_iter, g);
		if(n == u) {
			return edge_iter;
		}
	}
	return e[v].end();
}

template<
		typename graph, typename planar_embedding,
		typename vertex_descriptor = typename boost::graph_traits<graph>::vertex_descriptor,
		typename edge_iterator = typename boost::property_traits<planar_embedding>::value_type::const_iterator
	>
inline edge_iterator find_neighbor_iterator_restricted(vertex_descriptor v, vertex_descriptor u,
	edge_iterator begin, edge_iterator end, const planar_embedding & e, const graph & g)
{
	for(edge_iterator edge_iter = begin; edge_iter != end; ++edge_iter) {
		if(edge_iter == e[v].end()) edge_iter = e[v].begin();
		vertex_descriptor n = get_incident_vertex(v, *edge_iter, g);
		if(n == u) {
			return edge_iter;
		}
	}
	return e[v].end();
}

#endif