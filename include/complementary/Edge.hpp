#ifndef EDGE_COMP_HPP
#define EDGE_COMP_HPP

// Used to represent graphs and trees
class Edge {
  
private:
  unsigned int src; // Index of the source vertex of the edge
  unsigned int tgt; // Index of the target vertex of the edge
  int cmp; // Position of the complementary edge (in the adjacency list of
	   // tgt). cmp < 0 means that the field cmp is undefined

public:
  Edge(){}

  ~Edge(){}

  unsigned int getSrc() {
    return this->src;
  }
  
  unsigned int getTgt() {
    return this->tgt;
  }
  
  int getCmp(){
    return this->cmp;
  }

  void setSrc(unsigned int s) {
    this->src = s;
  }
  
  void setTgt(unsigned int t) {
    this->tgt = t;
  }
  
  void setCmp(int c) {
    this->cmp = c;
  }
};

#endif
