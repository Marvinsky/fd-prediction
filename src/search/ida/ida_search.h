#ifndef IDA_SEARCH_H
#define IDA_SEARCH_H

#define INT_MAX   2147483647


#include <vector>
#include <map>
#include <stack>
#include <set>
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

using namespace std;


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
        StateID get_id() const {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        double getHvalue() const {return this->h_value;}
        void setHvalue(double h) {this->h_value = h;}
        double getGreal() const {return this->g_real;}
        void setGreal(double g) {this->g_real = g;}
        double getLevel() const {return this->level;}
        void setLevel(double l) {this->level = l;}
};

bool fncomp (SSNode lhs, SSNode rhs) {

        if (lhs.get_id().hash() != rhs.get_id().hash()) {
                return lhs.get_id().hash() < rhs.get_id().hash();
        }
        return false;
}

struct classcomp {
        bool operator() (const SSNode& lhs, const SSNode& rhs) const {
		//cout<<"lhs.get_id().hash() = "<<lhs.get_id().hash();
		//cout<<"\trhs.get_id().hash() = "<<rhs.get_id().hash()<<"\n";
                if (lhs.get_id().hash() != rhs.get_id().hash()) {
			//bool flag = lhs.get_id().hash() < rhs.get_id().hash();
			//cout<<"\t\treturn "<<flag<<"\n\n";
                        return lhs.get_id().hash() < rhs.get_id().hash();
                }
                return false;
        }
};

class IDASearch : public SearchEngine {
private:

    //DFS
    int depth;
    string heuristic_name;
    set<SSNode, classcomp> buffer;
    set<SSNode, classcomp> check;
    map<Node2, double> expanded;
    map<Node2, double> generated;
    Timer ida_timer;

    //IDA*
    int best_soln_sofar;
    int nodes_expanded_for_bound;
    int nodes_generated_for_bound;
    int nodes_expanded_for_start_state;
    int nodes_generated_for_start_state;
    int next_bound;
    bool SOLUTION_FOUND;
    std::queue<SSNode> D;
    //int next_bound;

protected:
    SearchStatus step();
    void print_heuristic_values(const std::vector<int> &values) const;

    std::vector<Heuristic *> heuristics;

    virtual void initialize();
    
public:
    IDASearch(const Options &opts); 
    int idastar(SSNode node);
    int dfs_heur(SSNode node, double bound, double next_bound);
    //BFS
    set<SSNode, classcomp> BFS(SSNode root);
};

#endif
