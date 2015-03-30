#ifndef DIJKSTRA_SEARCH_HEURISTIC_H
#define DIJKSTRA_SEARCH_HEURISTIC_H

#include "heuristic.h"

class DijkstraSearchHeuristic : public Heuristic {
 
	 int min_operator_cost;
protected:

  	virtual void initialize();
  	virtual int compute_heuristic(const GlobalState &state);
public:

  	DijkstraSearchHeuristic(const Options &options);
  	~DijkstraSearchHeuristic();
  	virtual string get_heur_name() {string temp = "dijkstra"; return temp;}
};

#endif
