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

namespace boost {
	
	template<
		typename graph,
		typename vertex_descriptor = typename graph_traits<graph>::vertex_descriptor,
		typename edge_descriptor = typename graph_traits<graph>::edge_descriptor
	>
	vertex_descriptor get_incident_vertex(vertex_descriptor v, edge_descriptor e, const graph & g)
	{
		vertex_descriptor n = source(e, g);
		if(n == v)
			n = target(e, g);
		return n;
	}
	
}

#endif