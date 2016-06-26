/*
 * draw_tikz_graph.h
 * Author: Aven Bross
 *
 * Produce tikz drawings of boost graphs.
 */

#ifndef DRAW_TIKZ_GRAPH_HPP
#define DRAW_TIKZ_GRAPH_HPP

#include <string>
#include <utility>
#include <map>
#include <cmath>
#include <vector>
#include <stdexcept>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/tuple/tuple.hpp>

const std::vector<std::string> color_strings =
	{"red", "blue","yellow", "green", "cyan", "magenta", "white", "black"};

// Produce a dot language string to draw the given graph that has been three colored
template<typename index_graph, typename color_map, typename plane_drawing>
std::string draw_tikz_graph(const index_graph & graph, const color_map & coloring,
	const plane_drawing & drawing)
{
	typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
	typedef typename boost::graph_traits<index_graph>::edge_descriptor edge_descriptor;
	typedef typename boost::property_traits<color_map>::value_type color_type;
	
	std::map<color_type, std::string> color_string_map;
	
	typename boost::graph_traits<index_graph>::vertex_iterator v_iter, v_end;

    for (boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; v_iter++) {
		vertex_descriptor curr_vertex = *v_iter;
		color_type curr_color = coloring[curr_vertex];
		
		if(color_string_map.count(curr_color) == 0)
		{
			color_string_map[curr_color] = color_strings[curr_color];
		}
	}
	
	std::string tikz_string;
	
	tikz_string += "\\begin{tikzpicture}[scale=.7, every node/.style={circle,"
		+ std::string(" draw, minimum size=8.5mm, scale=.7}]\n");
	
	for (boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; v_iter++) {
		vertex_descriptor curr_vertex = *v_iter;
		std::string vname = std::to_string(curr_vertex);
		std::string vcolor = color_string_map[coloring[curr_vertex]];
		std::string x = std::to_string(get(drawing, curr_vertex).x);
		std::string y = std::to_string(get(drawing, curr_vertex).y);
		
		tikz_string += "  \\node [fill=" + vcolor + "!35] (" + vname
			+ ") at (" + x + "cm, " + y + "cm) {$" +  vname + "$};\n";
	}
	
	typename boost::graph_traits<index_graph>::edge_iterator e_iter, e_end;
	
	for (boost::tie(e_iter, e_end) = boost::edges(graph); e_iter != e_end; e_iter++) {
		edge_descriptor curr_edge = *e_iter;
		vertex_descriptor v1 = boost::source(curr_edge, graph);
		vertex_descriptor v2 = boost::target(curr_edge, graph);
		std::string v1name = std::to_string(v1);
		std::string v2name = std::to_string(v2);
		
		tikz_string += "  \\draw (" + v1name + ") -- (" + v2name + ");\n";
	}
	
	tikz_string += "\\end{tikzpicture}";
	
	return tikz_string;
}

#endif