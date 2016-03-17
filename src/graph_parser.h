/*
 * graph_parser.h
 * Author: Aven Bross
 *
 * Recursive descent parser for simplified dot language
 */

#ifndef GRAPH_PARSER_HPP
#define GRAPH_PARSER_HPP

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

const std::vector<std::string> color_strings = {"red", "blue","yellow"};

// Produce a dot language string to draw the given graph that has been three colored
template<typename Graph, typename Coloring, typename Drawing>
std::string drawTikzGraph(const Graph & graph, const Coloring & coloring,
	const Drawing & drawing)
{
	typedef boost::graph_traits<Graph> GraphTraits;
	typedef typename GraphTraits::vertex_descriptor vertex_descriptor;
	typedef typename GraphTraits::edge_descriptor edge_descriptor;
	typedef typename boost::property_traits<Coloring>::value_type color_type;
	
	std::map<color_type, std::string> color_map;
	
	typename GraphTraits::vertex_iterator v_iter, v_end;

    for (boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; v_iter++) {
		vertex_descriptor curr_vertex = *v_iter;
		color_type curr_color = coloring[curr_vertex];
		
		if(color_map.count(curr_color) == 0)
		{
			if(color_map.size() == 3)
			{
				throw std::logic_error("More than 3 colors in coloring");
			}
			if(color_map.count(curr_color) == 0)
			{
				color_map[curr_color] = color_strings[color_map.size()];
			}
		}
	}
	
	std::string tikz_string;
	
	tikz_string += "\\begin{tikzpicture}[scale=.7, every node/.style={circle,"
		+ std::string(" draw, minimum size=8.5mm, scale=.7}]\n");
	
	for (boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; v_iter++) {
		vertex_descriptor curr_vertex = *v_iter;
		std::string vname = std::to_string(curr_vertex);
		std::string vcolor = color_map[coloring[curr_vertex]];
		std::string x = std::to_string(get(drawing, curr_vertex).x);
		std::string y = std::to_string(get(drawing, curr_vertex).y);
		
		tikz_string += "  \\node [fill=" + vcolor + "!35] (" + vname
			+ ") at (" + x + "cm, " + y + "cm) {$" /* +  vname */ + "$};\n";
	}
	
	typename GraphTraits::edge_iterator e_iter, e_end;
	
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