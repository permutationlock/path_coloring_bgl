//=======================================================================
// Copyright 2007 Aaron Windsor
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>

#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/is_straight_line_drawing.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/graphviz.hpp>

#include "path_coloring.h"



using namespace boost;

//a class to hold the coordinates of the straight line embedding
struct coord_t
{
  std::size_t x;
  std::size_t y;
};


int main(int argc, char** argv)
{
  typedef adjacency_list
    < vecS,
      vecS,
      undirectedS,
      property<vertex_index_t, int>
    > graph;

  

  //Define the storage type for the planar embedding
  typedef std::vector< std::vector< graph_traits<graph>::edge_descriptor > > 
    embedding_storage_t;
  typedef boost::iterator_property_map
    < embedding_storage_t::iterator, 
      property_map<graph, vertex_index_t>::type 
    >
    embedding_t;



  // Create the graph - a maximal planar graph on 7 vertices. The functions
  // planar_canonical_ordering and chrobak_payne_straight_line_drawing both
  // require a maximal planar graph. If you start with a graph that isn't
  // maximal planar (or you're not sure), you can use the functions
  // make_connected, make_biconnected_planar, and make_maximal planar in
  // sequence to add a set of edges to any undirected planar graph to make
  // it maximal planar.
  
  graph g(3);
  add_edge(0,1,g);
  add_edge(1,2,g);
  add_edge(0,3,g);

  // Create the planar embedding
  embedding_storage_t embedding_storage(num_vertices(g));
  embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));

  boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                               boyer_myrvold_params::embedding = embedding);
  
  for(int i=0; i<3; i++)
  {                      
  	  std::cout << "Adjacency list for " << i << ".\n";       
	  for(auto iter = embedding[i].begin(); iter != embedding[i].end(); iter++)
	  {
	  	std::cout << "{" << source(*iter, g) << ", " << target(*iter, g)
	  		<< "} adj = " << getIncidentVertex(i,*iter,g) << "\n"; 
	  }
	  std::cout << "\n";
  }

  std::vector<graph::vertex_descriptor> p;
  std::vector<graph::vertex_descriptor> q;
  
  p.push_back(0);
  
  q.push_back(1);
  q.push_back(2);
  
  std::map<graph::vertex_descriptor,int> color_map;
  boost::associative_property_map< std::map<graph::vertex_descriptor, int> >
    color_property_map(color_map);
  
  int color = 0;
  
  poh_path_color(g, embedding, p.begin(), p.end(), q.begin(), q.end(), color_property_map, color);

  // Find a canonical ordering
  /*std::vector<graph_traits<graph>::vertex_descriptor> ordering;
  planar_canonical_ordering(g, embedding, std::back_inserter(ordering));


  //Set up a property map to hold the mapping from vertices to coord_t's
  typedef std::vector< coord_t > straight_line_drawing_storage_t;
  typedef boost::iterator_property_map
    < straight_line_drawing_storage_t::iterator, 
      property_map<graph, vertex_index_t>::type 
    >
    straight_line_drawing_t;

  straight_line_drawing_storage_t straight_line_drawing_storage
    (num_vertices(g));
  straight_line_drawing_t straight_line_drawing
    (straight_line_drawing_storage.begin(), 
     get(vertex_index,g)
     );



  // Compute the straight line drawing
  chrobak_payne_straight_line_drawing(g, 
                                      embedding, 
                                      ordering.begin(),
                                      ordering.end(),
                                      straight_line_drawing
                                      );
  


  std::cout << "The straight line drawing is: " << std::endl;
  graph_traits<graph>::vertex_iterator vi, vi_end;
  for(boost::tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
    {
      coord_t coord(get(straight_line_drawing,*vi));
      std::cout << *vi << " -> (" << coord.x << ", " << coord.y << ")" 
                << std::endl;
    }

  // Verify that the drawing is actually a plane drawing
  if (is_straight_line_drawing(g, straight_line_drawing))
    std::cout << "Is a plane drawing." << std::endl;
  else
    std::cout << "Is not a plane drawing." << std::endl;*/

  return 0;
}
