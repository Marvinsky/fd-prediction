#ifndef EAGER_SEARCH_H
#define EAGER_SEARCH_H

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

#include <iostream>
#include <fstream>

class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

class EagerSearch : public SearchEngine {
    // Search Behavior parameters
    bool reopen_closed_nodes; // whether to reopen closed nodes upon finding lower g paths
    bool do_pathmax; // whether to use pathmax correction
    bool use_multi_path_dependence;

    OpenList<StateID> *open_list;
    ScalarEvaluator *f_evaluator;


//A* prediction
    map<Node2, double> nodes_expanded;
    vector<Timer> v_timer;
    double count_nodes;
    bool first_time_in_solved;
    int F_boundary;
    double count_last_nodes_generated;
    Timer time_level;

//Velocity-Based Search Speed Estimator
    
    int initial_value;
    int total_min;
    int nodes_expanded_for_start_state;
    int nodes_generated_for_start_state;
    Timer search_time;
    Timer level_time; //time required to expand an entire level
    double target_search_velocity;
    double V; // Search velocity - it is calcultated based the number of nodes generated
    double search_speed; // Search velocity - it is calculated based the number of nodes expanded
    double SEv;  //Future search effort or search effort estimation
    double VeSP; //Velocity Search Progress Estimator

    ofstream outputFile2;
    
    //Vacillation-Based Search Speed Estimator   
    //TODO: There is no enough information about DAS
    
    //Path-Based Progress Estimator
    //Unit-cost Domains
    double NPBP; //Naive path-based progress estimator

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
    EagerSearch(const Options &opts);
    void statistics() const;

    void dump_search_space();

//A* prediction
    void generateExpandedReport();
    int returnMaxF(vector<int> levels);
    int returnMinF(vector<int> levels);


//Speed Progress
    void reportProgress();
    //int generatedSoFar();
    //int expandedSoFar();
};

#endif
