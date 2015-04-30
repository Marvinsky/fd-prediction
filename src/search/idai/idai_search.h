#ifndef IDAI_SEARCH_H
#define IDAI_SEARCH_H

#define INT_MAX   2147483647


#include <vector>
#include <map>
#include <stack>
#include <queue>

#include "../evaluator.h"
#include "../global_state.h"
#include "../search_engine.h"
#include "../search_progress.h"
#include "../search_space.h"
#include "../timer.h"

#include "../open_lists/open_list.h"
#include "../state_id.h"
#include "../ss/node2.h"

#include <random>
#include <iostream>

#include <limits>


class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

class Node {
private:
        StateID id;
        double h_value;
        double g_real;
        double level;
public:
        Node() : id(StateID::no_state), h_value(-1),  g_real(-1), level(-1) {}     
	Node(StateID identifier, double h, double g, double l) : id(identifier), h_value(h), g_real(g), level(l) {}
        StateID get_id() const {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        double getHvalue() const {return this->h_value;}
        void setHvalue(double h) {this->h_value = h;}
        double getGreal() const {return this->g_real;}
        void setGreal(double g) {this->g_real = g;}
        double getLevel() const {return this->level;}
        void setLevel(double l) {this->level = l;}
};

struct classcomp {
        bool operator() (const Node& lhs, const Node& rhs) const { 
                return lhs.get_id() < rhs.get_id(); 
        }
};

class IDAISearch : public SearchEngine {
private:

    //DFS
    int depth;
    string heuristic_name;
    Timer ida_timer;

    //IDA*
    int best_soln_sofar;
    int nodes_expanded_for_bound;
    int nodes_generated_for_bound;
    int nodes_expanded_for_start_state;
    int nodes_generated_for_start_state;
    //IDA* bfs
    set<Node, classcomp> check;
    set<Node, classcomp> L;


protected:
    SearchStatus step();
    void print_heuristic_values(const std::vector<int> &values) const;

    std::vector<Heuristic *> heuristics;

    virtual void initialize();
    
public:
    IDAISearch(const Options &opts); 
    int idastar(Node node);
    int dfs_heur(Node node, double bound, double &next_bound);
    void BFS(Node root);
};

#endif
