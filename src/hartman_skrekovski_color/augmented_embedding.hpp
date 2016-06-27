/*
 * augmented_embedding.hpp
 * Author: Aven Bross
 *
 * Implementation of a simple adjacency list embedding structure with reverse lookup.
 */

#ifndef __AUGMENTED_EMBEDDING_HPP
#define __AUGMENTED_EMBEDDING_HPP

// STL headers
#include <vector>

// Basic graph headers
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

/*
 * augmented_embedding
 * This structure takes a planar graph with vertex indices and a valid plane embedding (modelling
 * the boost PlanarEmbedding concept). It then constructs itself as an adjacency list graph, with lists
 * ordered according to the given plane embedding. Furthermore, each entry in a vertex's adjacency list
 * provides the index of its position in the corresponding the neighbor's ajdacency list.
 *
 * NOTE: Construction is done in O(|V|) time as edge lookup is amortized constant in plane graphs.
 */

template<typename index_graph>
class augmented_embedding {
	public:
		typedef typename boost::graph_traits<index_graph>::vertex_descriptor vertex_descriptor;
		
		// Stores the neighbor as well as the current vertex's position in the neighbor's list
		struct neighbor_data {
			neighbor_data(vertex_descriptor neighbor, std::size_t back_index)
				: neighbor(neighbor), back_index(back_index) {}
			
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
			
			// Construct embedded adjacency list and remember indices for each edge
			typename boost::graph_traits<index_graph>::vertex_iterator v_iter, v_end; 
			for(boost::tie(v_iter, v_end) = boost::vertices(graph); v_iter != v_end; ++v_iter) {
				vertex_descriptor v = *v_iter;
				
				for(auto e_iter = embedding[v].begin(); e_iter != embedding[v].end(); ++e_iter) {
					vertex_descriptor n = boost::source(*e_iter, graph);
					if(v == n) n = boost::target(*e_iter, graph);
					
					std::size_t count = 0;
					
					for(auto ne_iter = embedding[n].begin(); ne_iter != embedding[n].end(); ++ne_iter) {
						vertex_descriptor u = boost::source(*ne_iter, graph);
						if(n == u) u = boost::target(*ne_iter, graph);
					
						if(v == u) {
							ndata_map[v].push_back(neighbor_data(n, count));
							break;
						}
						++count;
					}
				}
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