/*
 * draw_tikz_graph.hpp
 * Author: Aven Bross
 *
 * Produce tikz drawings of boost graphs.
 */

#ifndef DRAW_TIKZ_GRAPH_HPP
#define DRAW_TIKZ_GRAPH_HPP

// STL headers
#include <string>
#include <utility>
#include <map>
#include <cmath>
#include <vector>
#include <stdexcept>

// Basic graph headers
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Vector of all default tikz colors (therefore use fewer than 8 colors)
const std::vector<std::string> color_strings =
	{"white", "red", "blue","yellow", "green", "cyan", "magenta", "black"};

// Produce string of LaTeX to draw the given graph and coloring
template<typename graph_t, typename color_map_t, typename plane_drawing_t>
std::string draw_tikz_graph(
		const graph_t & graph, const color_map_t & color_map,
		const plane_drawing_t & drawing
	)
{
	typedef typename boost::graph_traits<graph_t>::vertex_descriptor
		vertex_t;
	typedef typename boost::graph_traits<graph_t>::edge_descriptor
		edge_t;
	typedef typename boost::property_traits<color_map>::value_type color_type;
	
	std::map<color_type, std::string> color_string_map;
	
	typename boost::graph_traits<graph_t>::vertex_iterator v_iter, v_end;

    for (boost::tie(v_iter, v_end) = boost::vertices(graph);
    	v_iter != v_end; v_iter++)
	{
		vertex_t curr_vertex = *v_iter;
		color_type curr_color = color_map[curr_vertex];
		
		if(color_string_map.count(curr_color) == 0)
		{
			color_string_map[curr_color] = color_strings[curr_color];
		}
	}
	
	std::string tikz_string;
	
	tikz_string += "\\begin{tikzpicture}[scale=.7, every node/.style={circle,"
		+ std::string(" draw, minimum size=8.5mm, scale=.7}]\n");
	
	for (boost::tie(v_iter, v_end) = boost::vertices(graph);
		v_iter != v_end; v_iter++)
	{
		vertex_t curr_vertex = *v_iter;
		std::string vname = std::to_string(curr_vertex);
		std::string vcolor = color_string_map[color_map[curr_vertex]];
		std::string x = std::to_string(get(drawing, curr_vertex).x);
		std::string y = std::to_string(get(drawing, curr_vertex).y);
		
		tikz_string += "  \\node [fill=" + vcolor + "!50] (" + vname
			+ ") at (" + x + "cm, " + y + "cm) {$" +  vname + "$};\n";
	}
	
	typename boost::graph_traits<graph_t>::edge_iterator e_iter, e_end;
	
	for (boost::tie(e_iter, e_end) = boost::edges(graph);
		e_iter != e_end; e_iter++)
	{
		edge_t curr_edge = *e_iter;
		vertex_t v1 = boost::source(curr_edge, graph);
		vertex_t v2 = boost::target(curr_edge, graph);
		std::string v1name = std::to_string(v1);
		std::string v2name = std::to_string(v2);
		
		if(color_map[v1] == color_map[v2]) {
			std::string vcolor = color_string_map[color_map[v1]];
			tikz_string += "  \\draw (" + v1name + ") -- (" + v2name
				+ ") [ultra thick,color=black];\n";
			tikz_string += "  \\draw (" + v1name + ") -- (" + v2name
				+ ") [very thick,color=" + vcolor + "!50];\n";
		}
		else {
			tikz_string += "  \\draw (" + v1name + ") -- (" + v2name + ");\n";
		}
	}
	
	tikz_string += "\\end{tikzpicture}";
	
	return tikz_string;
}

#endif
