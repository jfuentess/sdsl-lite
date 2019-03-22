#include <fstream>
#include <cstdlib>
#include <string.h>
#include "Graph.hpp"

/*
  Format of the expected input/output file:
  
  <number of nodes>
  <number of edges>
  <source vertex> <target vertex>
  ....

  The list of edges (<source vertex> <target vertex>) must be ordered:
  - First, all edges where vertex 0 is the source vertex
  - Second, all edges where vertex 1 is the source vertex
  - and so on

  For a source vertex v, the corresponding target vertices must be in
  counterclockwise order

  Assumption: To support multiple edges, the indices of the multiple edges must
  be always increasing (in other words, the adjacency list of a node cannot
  start in the middle of a list of multiple edges
*/
Graph read_graph_from_file(const char* fn) {

  FILE *fp = fopen(fn, "r");
  char line_buffer[BUFSIZ]; /* BUFSIZ is defined if you include stdio.h */

  if (!fp) {
    fprintf(stderr, "Error opening file \"%s\".\n", fn);
    exit(EXIT_FAILURE);
  }

  unsigned int n = atoi(fgets(line_buffer, sizeof(line_buffer), fp));
  unsigned int m = atoi(fgets(line_buffer, sizeof(line_buffer), fp));

  Graph g(n, m);

  uint source = 0, target = 0, mm = 0;

  while (fgets(line_buffer, sizeof(line_buffer), fp)) {
    source = atoi(strtok(line_buffer, " "));
    target = atoi(strtok(NULL, " "));
    g.setEdge(mm, source, target); // default: E[mm].cmp = -1

    if(mm==0)
      g.setVertexFirst(source, mm);
    else if(source != g.getEdgeSrc(mm-1)) {
      g.setVertexLast(g.getEdgeSrc(mm-1), mm-1);
      g.setVertexFirst(source, mm);
    }
    mm++;
  }
  g.setVertexLast(g.getEdgeSrc(mm-1), mm-1);
  fclose(fp);

  for(uint i = 0; i < 2*m; i++) {
    Vertex target = g.getVertex(g.getEdgeTgt(i));
    int cmp = -1;

    if(g.getEdgeCmp(i) != -1)
      continue;
    
    for(uint j = target.getFirst(); j <= target.getLast(); j++) {
      /* Condition i != j allows to support self-loops */
      if((g.getEdgeCmp(j) == -1) && (g.getEdgeTgt(j) == g.getEdgeSrc(i)) &&
	 (i != j))
	cmp = j; // Choose the last unvisited (e.cmp=-1) edge, not the first one
    }

    if(cmp != -1) {
      g.setEdgeCmp(i, cmp);
      g.setEdgeCmp(cmp, i);
    }
  }

  return g;
}
