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
#include <stack>

#include <limits>


class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

class SSNode {
private:
        StateID id;
        int h_value;
        int g_real;
        int level;
public:
        SSNode() : id(StateID::no_state), h_value(-1),  g_real(-1), level(-1) {}     
	SSNode(StateID identifier,int h, int g, int l) : id(identifier), h_value(h), g_real(g), level(l) {}
        StateID get_id() {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        int getHvalue() {return this->h_value;}
        void setHvalue(int h) {this->h_value = h;}
        int getGreal() {return this->g_real;}
        void setGreal(int g) {this->g_real = g;}
        int getLevel() {return this->level;}
        void setLevel(int l) {this->level = l;}
};

class IDASearch : public SearchEngine {
private:

    //DFS
    int depth;
    string heuristic_name;
    stack<SSNode> queue;
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

    void generateGeneratedReport(bool flag);
    void generateExpandedReport(bool flag);
    int idastar(SSNode node);
    int dfs_heur(SSNode node, int bound, int &next_bound, int g_real);
};

#endif
