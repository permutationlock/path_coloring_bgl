/*
 * graph_parser.cpp
 * Author: Aven Bross
 *
 * Recursive descent parser for simplified dot language
 */

#include "graph_parser.h"

// Lex input expression into tokens
GraphParser::GraphParser(const std::string & input)
{
	bool concat_node = false;
	std::string value;
	size_t i = 0;

	while(i < input.size())
	{
		char c = input[i++];
		
		if(c == '\t' || c == ' ' || c == '\n')
		{
			continue;
		}
		else if(concat_node)
		{
			if((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')|| (c >= '1' && c <= '9'))
			{
				value += c;
				continue;
			}

			tokens.push_back(Token(NODE,value));
			value = "";
			concat_node = false;
		}
		else if(c == '-' && input[i] == '-')
		{
			tokens.push_back(Token(EDGE,"--"));
		}
		else if(c == '}')
		{
			tokens.push_back(Token(RBRKT,"}"));
		}
		else if(c == '{')
		{
			tokens.push_back(Token(LBRKT,"{"));
		}
		else if(c == 's')
		{
			if(!input.substr(i,5).compare("trict"))
			{
				tokens.push_back(Token(STRCT,"strict"));
				i += 5;
			}
			else
			{
				value += c;
				concat_node = true;
			}
		}
		else if(c == 'g')
		{
			if(!input.substr(i,4).compare("raph"))
			{
				tokens.push_back(Token(GRPH,"graph"));
				i += 4;
			}
			else
			{
				value += c;
				concat_node = true;
			}
		}
		else if((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')|| (c >= '1' && c <= '9'))
		{
			value += c;
			concat_node = true;
		}
		else
		{
			tokens.push_back(Token(INVALID,std::string(c,1)));
		}
	}

	// Extra check to catch single number values
	if(value.compare("")) tokens.push_back(Token(NODE,value));
}

// Check next token and return if it matches the type
// If true it advances index to next token
bool GraphParser::match(TokenType the_type)
{
	if(index == tokens.size()) return false;

	Token t = tokens[index];

	if(t.first == the_type || the_type == ANY)
	{
		index++;
		cur_val = t.second;
		return true;
	}

	return false;
}

// Parse: start -> graph + start | e
// and:   graph -> STRCT GRPH LBRKT edge_list RBRKT
std::vector<Graph<GraphParser::Vertex> > GraphParser::parseGraphs()
{
	std::vector<Graph<Vertex> > graphs;
	
	while(match(STRCT))
	{
		if(match(GRPH))
		{
			if(match(LBRKT))
			{
				std::vector<Edge> edges;
				std::vector<Vertex> vertices;
				
				if(match(NODE))
				{
					Vertex v1(cur_val);
					if(match(EDGE) && match(NODE))
					{
						Vertex v2(cur_val);
						edges.push_back(Edge(v1,v2));
						vertices.push_back(v1);
						vertices.push_back(v2);
					}
					else
					{
						vertices.push_back(v1);
					}
				}
				
				graphs.push_back(Graph<Vertex>(vertices,edges));
				
				if(!match(RBRKT))
				{
					throw std::runtime_error("Missing ending \'}\'.");
				}
			}
			else
			{
				throw std::runtime_error("Missing initial \'{\'.");
			}
		}
		else if(match(ANY))
		{
			throw std::runtime_error("Invalid token \'" + cur_val + "\'.");
		}
	}
	
	if(match(ANY))
	{
		throw std::runtime_error("Invalid token \'" + cur_val + "\'.");
	}

	return graphs;
}

// Produce a dot language string to draw the given graph that has been three colored
template<typename Label, typename Color>
std::string generateGraph(const Graph<Label> & graph, const std::unordered_map<Label,Color> & coloring)
{
	std::map<Color, std::string> color_map;
	std::map<Label, std::string> label_map;
	
	// Generate strings for each label and color
	for(const auto & vertex_pair : graph)
	{
		Label curr_vertex = vertex_pair.first;
		Color curr_color = coloring.at(curr_vertex);
		
		if(label_map.count(curr_vertex) == 0)
		{
			label_map[curr_vertex] = std::to_string(label_map.size());
		}
		
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
	
	std::unordered_set<Label> visited;
	std::queue<Label> bfs_queue;
	std::string dot_string = "strict graph{\nnode [shape=circle]\nsplines=spline\n";
	
	visited.insert(graph.begin() -> first);
	bfs_queue.push(graph.begin() -> first);
	
	while(!bfs_queue.empty())
	{
		Label curr_vertex = bfs_queue.front();
		bfs_queue.pop();
		Color curr_color = coloring.at(curr_vertex);
		
		dot_string += label_map[curr_vertex] + " [color=" + color_map[curr_color] + "]" + "\n";
		
		for(const Label & neighbor : graph[curr_vertex])
		{
			dot_string += label_map[curr_vertex] + " -- " + label_map[neighbor] + "\n";
			
			if(visited.count(neighbor) == 0)
			{
				visited.insert(neighbor);
				bfs_queue.push(neighbor);
			}
		}
	}
	
	dot_string += "}";
	
	return dot_string;
}

