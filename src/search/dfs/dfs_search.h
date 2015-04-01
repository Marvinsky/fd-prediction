#ifndef DFS_SEARCH_H
#define DFS_SEARCH_H

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
#include <stack>

class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

class SSNode {
private:
        StateID id;
	int h_value;
        int g_value;
public:
        SSNode() : id(StateID::no_state), h_value(-1), g_value(-1) {}     
	SSNode(StateID identifier, int h, int g) : id(identifier), h_value(h), g_value(g) {}
        StateID get_id() {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        int get_h_value() {return this->h_value;}
        void set_h_value(int h) {this->h_value = h;}
        int get_g_value() {return this->g_value;}
        void set_g_value(int g) {this->g_value = g;}
};

class DFSSearch : public SearchEngine {
private:
    // Search Behavior parameters
    bool reopen_closed_nodes; // whether to reopen closed nodes upon finding lower g paths
    bool do_pathmax; // whether to use pathmax correction
    bool use_multi_path_dependence;
    bool mark_children_as_finished;

    OpenList<StateID> *open_list;
    ScalarEvaluator *f_evaluator;

    //DFS
    int depth;
    string heuristic_name;
    stack<SSNode> queue;
    map<Node2, double> collector;

protected:
    SearchStatus step();
    std::pair<SearchNode, bool> fetch_next_node();
    void update_jump_statistic(const SearchNode &node);
    void print_heuristic_values(const std::vector<int> &values) const; 
    void reward_progress();

    std::vector<Heuristic *> heuristics;
    std::vector<Heuristic *> preferred_operator_heuristics;
    std::vector<Heuristic *> estimate_heuristics;

    virtual void initialize();
    
public:
    DFSSearch(const Options &opts);
    void statistics() const; 

    void dump_search_space();
};

#endif
