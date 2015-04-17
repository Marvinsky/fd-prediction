#ifndef IDA_SEARCH_H
#define IDA_SEARCH_H

#define INT_MAX   2147483647


#include <vector>
#include <map>

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
#include <queue>

#include <limits>


class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

class SSNode {
private:
        StateID id;
        double h_value;
        double g_real;
        double level;
public:
        SSNode() : id(StateID::no_state), h_value(-1),  g_real(-1), level(-1) {}     
	SSNode(StateID identifier, double h, double g, double l) : id(identifier), h_value(h), g_real(g), level(l) {}
        StateID get_id() {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        double getHvalue() {return this->h_value;}
        void setHvalue(double h) {this->h_value = h;}
        double getGreal() {return this->g_real;}
        void setGreal(double g) {this->g_real = g;}
        double getLevel() {return this->level;}
        void setLevel(double l) {this->level = l;}
};

class IDASearch : public SearchEngine {
private:

    //DFS
    int depth;
    string heuristic_name;
    queue<SSNode> buffer;
    map<Node2, double> expanded;
    map<Node2, double> generated;
    Timer ida_timer;

    //IDA*
    int best_soln_sofar;
    int nodes_expanded_for_bound;
    int nodes_generated_for_bound;
    int nodes_expanded_for_start_state;
    int nodes_generated_for_start_state;
    //int next_bound;

protected:
    SearchStatus step();
    void print_heuristic_values(const std::vector<int> &values) const;

    std::vector<Heuristic *> heuristics;

    virtual void initialize();
    
public:
    IDASearch(const Options &opts);

    int idastar(SSNode node);
    int dfs_heur(SSNode node, double bound, double &next_bound, double g_real);

    //BFS
    queue<SSNode> BFS(SSNode root);        
};

#endif
