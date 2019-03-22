#ifndef GRAPH_COMP_HPP
#define GRAPH_COMP_HPP

#include <iostream>
#include <stack>
#include "Tree.hpp"

using namespace std;

class Graph {

private:
  Vertex *V;      // Array of vertices of the tree
  Edge *E;        // Array of edges of the tree. It is the concatenation of the
	          // adjacency lists of all nodes 
  unsigned int n; // Number of vertices in the tree
  unsigned int m; // Number of edges in the tree

public:
  Graph () {
    this->n = 0;
    this->m = 0;
  }

  Graph (unsigned int n, unsigned int m) {
    this->n = n;
    this->m = m;
    this->V = new Vertex[n];
    this->E = new Edge[2*m];
  }

  unsigned int vertices() {
    return n;
  }

  unsigned int edges() {
    return m;
  }

  // E[i].src = src and E[i].tgt = tgt
  void setEdge(int i, unsigned int src, unsigned int tgt) {
    this->E[i].setSrc(src);
    this->E[i].setTgt(tgt);
    this->E[i].setCmp(-1);
  }

  // V[i].first = first and V[i].last = last
  void setVertex(int i, unsigned int first, unsigned int last) {
    this->V[i].setFirst(first);
    this->V[i].setLast(last);
  }

  // V[i].first = first
  void setVertexFirst(int i, unsigned int first) {
    this->V[i].setFirst(first);
  }

  // V[i].last = last
  void setVertexLast(int i, unsigned int last) {
    this->V[i].setLast(last);
  }

  void setEdgeSrc(int i, unsigned int s) {
    this->E[i].setSrc(s);
  }

  void setEdgeTgt(int i, unsigned int t) {
    this->E[i].setTgt(t);
  }

  void setEdgeCmp(int i, int c) {
    this->E[i].setCmp(c);
  }

  unsigned int getEdgeSrc(int i) {
    return this->E[i].getSrc();
  }

  unsigned int getEdgeTgt(int i) {
    return this->E[i].getTgt();
  }

  int getEdgeCmp(int i) {
    return this->E[i].getCmp();
  }

  Vertex getVertex(int i) {
    return this->V[i];
  }

  Edge getEdge(int i) {
    return this->E[i];
  }

  unsigned int getVertexFirst(int i) {
    return this->V[i].getFirst();
  }

  unsigned int getVertexLast(int i) {
    return this->V[i].getLast();
  }

  Tree dfs_spanning_tree(unsigned int init, int *parent, unsigned int
				*count_edges, unsigned int *references) {
    unsigned int n = this->vertices();
    unsigned int m = this->edges();
    Tree t(n);

    
    char *visited = new char[n](); // TODO: Change to a boolean array
    unsigned int *edges = new unsigned int[2*m]();
    unsigned int num_tree_edges = 2*(n-1);

    stack <unsigned int> s;
    visited[init] = 1;
    s.push(init);
    parent[init] = -1;
  
    while(!s.empty()) {
      unsigned int curr = s.top(); s.pop();
    
      for(unsigned int i = this->V[curr].getFirst(); i <= this->V[curr].getLast(); i++)
	{
	  unsigned int tgt = this->E[i].getTgt();
      
	  if(visited[tgt] == 0) { // Not visited	
	    visited[tgt] = 1;
	    s.push(tgt);
	    parent[tgt] = this->E[i].getCmp(); // Edge child-to-parent
	  }
	}
    }
    
    // Marking the edges of G that are in T
    for(unsigned int i = 0; i < init; i++) {
      edges[(unsigned int)parent[i]] = 1;
      edges[(unsigned int)this->E[parent[i]].getCmp()] = 1;
    }
    for(unsigned int i = init+1; i < n; i++) {
      edges[(unsigned int)parent[i]] = 1;
      edges[(unsigned int)this->E[parent[i]].getCmp()] = 1;
    }

    
    unsigned int mm = 0;
    /* Counting edges */
    for(unsigned int i = 0; i < n; i++) {
      unsigned int first = this->V[i].getFirst();
      unsigned int last = this->V[i].getLast();
      unsigned int zeros = 0;
      int last_one = -1;
      unsigned int carry_last = 0;
    
      for(unsigned int j = first; j <= last; j++) {   
	if(edges[j] == 1) {
	  t.setEdge(mm, this->E[j].getSrc(), this->E[j].getTgt());

	  if(last_one == -1)
	    carry_last = zeros;	
	  else
	    count_edges[last_one] += zeros;

	  references[mm] = j;
	  last_one = mm;
	  zeros = 0;

	  if(mm == 0)
	    t.setNodeFirst(t.getEdgeSrc(mm), mm);
	  else if(t.getEdgeSrc(mm) != t.getEdgeSrc(mm-1)) {	
	    t.setNodeLast(t.getEdgeSrc(mm-1), mm-1);
	    t.setNodeFirst(t.getEdgeSrc(mm), mm);
	  }
	  mm++;
	}
	else
	  zeros++;
	edges[j] = mm;
      }

      count_edges[last_one] += zeros + carry_last;
      t.setNodeLast(t.getEdgeSrc(mm-1), mm-1);
    }

    for(unsigned int i = 0; i < init; i++)
      parent[i] = edges[parent[i]]-1;
    for(unsigned int i = init+1; i < n; i++)
      parent[i] = edges[parent[i]]-1;
  
    for(unsigned int i = 0; i < num_tree_edges; i++) {
      Vertex tgt = t.getNode(t.getEdgeTgt(i));
    
      for(unsigned int j = tgt.getFirst(); j <= tgt.getLast(); j++)
    	if(t.getEdgeTgt(j) == t.getEdgeSrc(i)) {
    	  t.setEdgeCmp(i, j);
    	  break;
    	}
    }
    
    return t;
  }

  void connected_graph() {
    unsigned int n = this->vertices();
    stack <unsigned int> s;
    char *visited = new char[n]();
    int curr = -1;
    unsigned int edge;
    int first = 1;
    int num_vertices = 0;

    while(!s.empty() || first) {
      if(first) { // Root
	curr = 0;
	first = 0;
      }
      else {
	edge = s.top(); s.pop();
	curr = this->E[edge].getTgt(); // current
      }
      visited[curr] = 1;
      for(unsigned int i = this->V[curr].getFirst(); i <= this->V[curr].getLast();
	  i++) {
	if(!visited[this->E[i].getTgt()])
	  s.push(i);      
      }
    }
  
    for(unsigned int i = 0; i < n; i++)
      if(visited[i] == 0) { // There are unvisited vertices
	num_vertices++;
      }
    
    cout << "unvisited vertices: " << num_vertices << ", visited vertices: " <<
      n-num_vertices << endl;
  }
};

#endif
