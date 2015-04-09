#ifndef EAGER_DIJKSTRA_SEARCH_H
#define EAGER_DIJKSTRA_SEARCH_H

#include <vector>

#include "../evaluator.h"
#include "../global_state.h"
#include "../search_engine.h"
#include "../search_progress.h"
#include "../search_space.h"
#include "../timer.h"

#include "../open_lists/open_list.h"

//A* prediction
#include "../ss/node2.h"
#include <map>

class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

class EagerDijkstraSearch : public SearchEngine {
    // Search Behavior parameters
    bool reopen_closed_nodes; // whether to reopen closed nodes upon finding lower g paths
    bool do_pathmax; // whether to use pathmax correction
    bool use_multi_path_dependence;

    OpenList<StateID> *open_list;
    ScalarEvaluator *f_evaluator;


//Dijkstra Algorithm
    map<Node2, double> nodes_expanded;
    vector<double> v_timer;
    double count_nodes;
    bool first_time_in_solved;
    int F_boundary;
    double count_last_nodes_generated;
    Timer time_level;
 
    int nodes_expanded_for_start_state;
    int nodes_generated_for_start_state;
protected:
    SearchStatus step();
    std::pair<SearchNode, bool> fetch_next_node();
    void update_jump_statistic(const SearchNode &node);
    void print_heuristic_values(const std::vector<int> &values) const;
    void reward_progress();

    std::vector<Heuristic *> heuristics;
    std::vector<Heuristic *> preferred_operator_heuristics;
    std::vector<Heuristic *> estimate_heuristics;
    // TODO: in the long term this
    // should disappear into the open list

    virtual void initialize();

public:
    EagerDijkstraSearch(const Options &opts);
    void statistics() const;

    void dump_search_space();

//A* prediction
    void generateExpandedReport();
    int returnMaxF(vector<int> levels);
    int returnMinF(vector<int> levels);
};

#endif
