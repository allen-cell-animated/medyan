//
//  Graph.h
//  Cyto
//
//  Created by James Komianos on 11/11/14.
//  Copyright (c) 2014 University of Maryland. All rights reserved.
//

#ifndef __Cyto__Graph__
#define __Cyto__Graph__

#include <stdio.h>

#include <boost/config.hpp>
#include <boost/version.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/static_assert.hpp>

//using namespace boost;

/// definition of basic boost::graph properties
enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };
namespace boost {
    BOOST_INSTALL_PROPERTY(vertex, properties);
    BOOST_INSTALL_PROPERTY(edge, properties);
}

///Graph class is a template for a undirected graph

template < typename VERTEXPROPERTIES, typename EDGEPROPERTIES >
class Graph
{
public:
    ///Typedefs for easier access
    typedef boost::adjacency_list< boost::listS, boost::hash_setS, boost::undirectedS,
    boost::property<vertex_properties_t, VERTEXPROPERTIES>,
    boost::property<edge_properties_t, EDGEPROPERTIES>> GraphContainer;
    
    /* a bunch of graph-specific typedefs */
    typedef typename boost::graph_traits<GraphContainer>::vertex_descriptor Vertex;
    typedef typename boost::graph_traits<GraphContainer>::edge_descriptor Edge;
    
    typedef typename boost::graph_traits<GraphContainer>::vertex_iterator vertex_iter;
    typedef typename boost::graph_traits<GraphContainer>::edge_iterator edge_iter;
    typedef typename boost::graph_traits<GraphContainer>::adjacency_iterator adjacency_iter;
    
    typedef typename boost::graph_traits<GraphContainer>::degree_size_type degree_t;
    
    typedef std::pair<adjacency_iter, adjacency_iter> adjacency_vertex_range_t;
    typedef std::pair<edge_iter, edge_iter> edge_range_t;
    typedef std::pair<vertex_iter, vertex_iter> vertex_range_t;
    
    ///Default constructor, does nothing
    Graph() {}
    
    ///Copy constructor
    Graph(const Graph& g) : graph(g.graph) {}
    
    /// Assignment operator
    Graph& operator=(const Graph &rhs) {
        graph = rhs.graph;
        return *this;
    }
    
    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
    /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
    /// (as of gcc 4.703), and will presumbaly be fixed in the future.
    virtual ~Graph() noexcept {}
    
    ///clear the graph
    void clear() { graph.clear(); }
    
    ///Add a vertex with given properties
    Vertex addVertex(const VERTEXPROPERTIES& prop) {
        Vertex v = add_vertex(graph);
        properties(v) = prop;
        return v;
    }
    
    ///Remove a vertex
    void removeVertex(const Vertex& v) {
        clear_vertex(v, graph);
        remove_vertex(v, graph);
    }
    
    ///Add an edge to two vertices
    Edge addEdge(const Vertex& v1, const Vertex& v2, const EDGEPROPERTIES& prop) {
        Edge addedEdge = add_edge(v1, v2, graph).first;
        properties(addedEdge) = prop;
        
        return addedEdge;
    }
    
    ///Accecss vertex and edge properties
    VERTEXPROPERTIES& properties(const Vertex& v) {
        typename boost::property_map<GraphContainer, vertex_properties_t>::type param = get(vertex_properties, graph);
        return param[v];
    }
    
    const VERTEXPROPERTIES& properties(const Vertex& v) const {
        typename boost::property_map<GraphContainer, vertex_properties_t>::const_type param = get(vertex_properties, graph);
        return param[v];
    }
    
    EDGEPROPERTIES& properties(const Edge& v) {
        typename boost::property_map<GraphContainer, edge_properties_t>::type param = get(edge_properties, graph);
        return param[v];
    }
    
    const EDGEPROPERTIES& properties(const Edge& v) const {
        typename boost::property_map<GraphContainer, edge_properties_t>::const_type param = get(edge_properties, graph);
        return param[v];
    }
    
    ///getters for graph
    vertex_range_t getVertices() const { return vertices(graph); }
    edge_range_t getEdges() const { return edges(graph); }
    adjacency_vertex_range_t getAdjacentVertices(const Vertex& v) const { return adjacent_vertices(v, graph); }
    
    int getVertexDegree(const Vertex& v) const { return out_degree(v, graph); }
    
    GraphContainer& getGraph() {return graph;}
    
    
protected:
    GraphContainer graph;
};

#endif /* defined(__Cyto__Graph__) */

