/*
 * augmented_embedding.hpp
 * Author: Aven Bross
 *
 * Implementation of a simple adjacency list embedding structure with reverse lookup.
 */

#ifndef __AUGMENTED_EMBEDDING_HPP
#define __AUGMENTED_EMBEDDING_HPP

#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

// Stores embedding as an adjacency list that allows back lookup
template<typename index_graph>
class augmented_embedding {
	public:
		typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
		
		// Stores the neighbor as well as the current vertex's position in the neighbor's list
		struct neighbor_data {
			neighbor_data(vertex_descriptor neighbor) : neighbor(neighbor) {}
			
			vertex_descriptor neighbor;
			std::size_t back_index;
		};
		
		typedef typename std::vector<neighbor_data>::const_iterator neighbor_data_iterator;
		typedef std::vector<std::vector<neighbor_data> > augmented_embedding_storage;
		typedef boost::iterator_property_map<
				typename augmented_embedding_storage::iterator,
				typename boost::property_map<index_graph, boost::vertex_index_t>::const_type
			> augmented_embedding_map;
		
		// Construct augmented embedding from the given graph and plane embedding
		template<typename planar_embedding>
		augmented_embedding(const index_graph & graph, const planar_embedding & embedding)
		{
			// Initialize our property map
			ndata_storage = augmented_embedding_storage(boost::num_vertices(graph));
			ndata_map = augmented_embedding_map(ndata_storage.begin(), boost::get(boost::vertex_index, graph));
			
			// Temporary mappping from edges to their indices in their endpoint's adjacency lists
			typedef std::vector<std::pair<std::size_t, std::size_t> > edge_index_pair_storage;
			typedef boost::iterator_property_map<
					std::vector<std::pair<std::size_t, std::size_t> >::iterator,
					typename boost::property_map<index_graph, boost::edge_index_t>::const_type
				> edge_index_pair_map;
			
			edge_index_pair_storage edge_order_storage(boost::num_edges(graph));
			edge_index_pair_map edge_order_map(edge_order_storage.begin(),
				boost::get(boost::edge_index, graph));
			
			// Construct embedded adjacency list and remember indices for each edge
			typename boost::graph_traits<index_graph>::vertex_iterator v_iter, v_end; 
			for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; ++v_iter) {
				vertex_descriptor vertex = *v_iter;
				
				auto ordering = embedding[vertex];
				for(auto edge_iter = ordering.begin(); edge_iter != ordering.end(); ++edge_iter) {
					vertex_descriptor source_vertex = boost::source(*edge_iter, graph);
					vertex_descriptor target_vertex = boost::target(*edge_iter, graph);
					if(vertex == source_vertex)
					{
						if(source_vertex > target_vertex)
							edge_order_map[*edge_iter].second = ndata_map[vertex].size();
						else
							edge_order_map[*edge_iter].first = ndata_map[vertex].size();
						ndata_map[vertex].push_back(neighbor_data(target_vertex));
					}
					else {
						if(source_vertex > target_vertex)
							edge_order_map[*edge_iter].first = ndata_map[vertex].size();
						else
							edge_order_map[*edge_iter].second = ndata_map[vertex].size();
						ndata_map[vertex].push_back(neighbor_data(source_vertex));
					}
				}
			}
			
			// Set up back index lookups by iterating through each edge
			typename boost::graph_traits<index_graph>::edge_iterator edge_iter, edge_end; 
			for(boost::tie(edge_iter, edge_end) = boost::edges(graph); edge_iter != edge_end; ++edge_iter) {
				vertex_descriptor source_vertex = boost::source(*edge_iter, graph);
				vertex_descriptor target_vertex = boost::target(*edge_iter, graph);
				
				std::size_t source_index = edge_order_map[*edge_iter].first;
				std::size_t target_index = edge_order_map[*edge_iter].second;
				
				if(source_vertex > target_vertex) {
					vertex_descriptor temp = target_vertex;
					target_vertex = source_vertex;
					source_vertex = temp;
				}
				
				ndata_map[source_vertex][source_index].back_index = target_index;
				ndata_map[target_vertex][target_index].back_index = source_index;
			}
		}
		
		const std::vector<neighbor_data> & operator[](vertex_descriptor vertex) const {
			return ndata_map[vertex];
		}
		
	private:
		augmented_embedding_storage ndata_storage;
		augmented_embedding_map ndata_map;
};

#endif