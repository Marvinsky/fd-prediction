#ifndef BFS_SEARCH_H
#define BFS_SEARCH_H

#include <vector>
#include <map>
#include <set>

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
        StateID get_id() const {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        int getHvalue() const {return this->h_value;}
        void setHvalue(int h) {this->h_value = h;}
        int getGreal() const {return this->g_real;}
        void setGreal(int g) {this->g_real = g;}
        int getLevel() const {return this->level;}
        void setLevel(int l) {this->level = l;}
};

struct classcomp {
        bool operator() (const SSNode& lhs, const SSNode& rhs) const {
                if (lhs.get_id().hash() != rhs.get_id().hash()) { 
                        return lhs.get_id().hash() < rhs.get_id().hash();
                }
                return false;
        }
};

class BFSSearch : public SearchEngine {
private:

    //DFS
    int depth;
    string heuristic_name;
    queue<SSNode> buffer;
    set<SSNode, classcomp> check;
    map<Node2, double> expanded;
    map<Node2, double> generated;
    Timer ida_timer;

protected:
    SearchStatus step();
    void print_heuristic_values(const std::vector<int> &values) const;

    std::vector<Heuristic *> heuristics;

    virtual void initialize();
    
public:
    BFSSearch(const Options &opts);

    void generateGeneratedReport(bool flag);
    void generateExpandedReport(bool flag);
};

#endif
