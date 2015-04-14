#include "bfs_search.h"

#include "../globals.h"
#include "../heuristic.h"
#include "../option_parser.h"
#include "../successor_generator.h"
#include "../g_evaluator.h"
#include "../sum_evaluator.h"
#include "../plugin.h"

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <set>
#include "../global_state.h"

#include <iostream>
#include <fstream>

using namespace std;

BFSSearch::BFSSearch(
    const Options &opts)
    : SearchEngine(opts) {
    
        ScalarEvaluator * evaluator = opts.get<ScalarEvaluator *>("eval");
        std::set<Heuristic *> hset;
        evaluator->get_involved_heuristics(hset);
        for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); it++) {
                heuristics.push_back(*it);
        }
        assert(heuristics.size() == 1);

        int max_h = 0;
        for (size_t i = 0; i < heuristics.size(); i++) {
            heuristics[i]->evaluate(g_initial_state());
            int aux = heuristics[i]->get_heuristic();
            if (max_h < aux) {
                max_h = aux;
            }
        }
        cout<<"max_h(constructor) = "<<max_h<<endl;
}

void BFSSearch::initialize() {
    //DFS
    ida_timer.reset();

}

SearchStatus BFSSearch::step() {

        const GlobalState &initial_state = g_initial_state();

        int hmax = 0; 
        for (size_t i = 0; i < heuristics.size(); i++) {
            heuristics[i]->evaluate(initial_state);
            hmax = max(hmax, heuristics[i]->get_heuristic());
        }
        int h_initial = hmax;
        cout<<"h_initial = "<<h_initial<<endl;
        
        StateID initial_state_id = initial_state.get_id();
        SSNode node(initial_state_id, h_initial,  0, 0); //set the h_value, g_real and the level
    
        depth = 2*h_initial;
	cout<<"depth ="<<depth<<endl;
        buffer.push(node);
        double count_nodes = 0.0;

        while (!buffer.empty()) {
              SSNode nodecp = buffer.front();
              int g_real = nodecp.getLevel();
              int level = nodecp.getGreal();
              cout<<"\tRaiz: h = "<<nodecp.getHvalue()<<", g_real = "<<g_real<<", f = "<<nodecp.getHvalue() + g_real<<", level = "<<level<<"\n";
	      buffer.pop();

              StateID state_id = nodecp.get_id();

             
              std::vector<const GlobalOperator *> applicable_ops;

              GlobalState global_state = g_state_registry->lookup_state(nodecp.get_id());                
              g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
              cout<<"---------------------Begin Childs-----------------------\n";
	      for (size_t i = 0; i < applicable_ops.size(); ++i) {
                  const GlobalOperator *op = applicable_ops[i];
                  GlobalState child =  g_state_registry->get_successor_state(global_state, *op);
                  int hmax_value = 0; 
                  for (size_t i = 0; i < heuristics.size(); i++) {
            	      heuristics[i]->evaluate(child);
            	      hmax_value = max(hmax_value, heuristics[i]->get_heuristic());
        	  }	   

                  int succ_h = hmax_value;

                  search_progress.inc_generated();


                  SSNode succ_node(child.get_id(), succ_h, g_real + get_adjusted_cost(*op), level + 1);

                  cout<<"\t\tChild_"<<(i+1)<<" : h = "<<succ_h<<", g = "<<g_real + get_adjusted_cost(*op)<<", f = "<<succ_h + g_real + get_adjusted_cost(*op)<<"\n";
                  if (succ_h + succ_node.getGreal() <= depth) {
                       cout<<"\t\t\tNode generated: h = "<<succ_h<<", g = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<"\n";
                     buffer.push(succ_node);
                  } else {
                     cout<<"\t\t\tpruned!\n";
                  }
		  cout<<"\t\tEnd child_"<<(i+1)<<"\n";
              }
	      cout<<"------------------------End Childs--------------------\n";
        }
 
   	return SOLVED;
}

void BFSSearch::generateGeneratedReport(bool flag) {
        double nodes_total_generated = 0;
        for (map<Node2, double>::iterator iter = generated.begin(); iter != generated.end(); iter++) {

            double q = iter->second;
            nodes_total_generated += q;
        }

        if (flag) {
                cout<<"Total of nodes generated: "<<nodes_total_generated<<endl;
        } else {
                cout<<", Parcial of nodes generated: "<<nodes_total_generated<<endl;
        }

}

void BFSSearch::generateExpandedReport(bool flag) {
        double nodes_total_expanded = 0;
        for (map<Node2, double>::iterator iter = expanded.begin(); iter != expanded.end(); iter++) {
            double q = iter->second;
            nodes_total_expanded += q;
        }

        if (flag) {
                cout<<"Total of nodes expanded: "<<nodes_total_expanded<<endl;
        } else {
                cout<<", Parcial of nodes expanded: "<<nodes_total_expanded<<endl;
        }
}


void BFSSearch::print_heuristic_values(const vector<int> &values) const {
    for (size_t i = 0; i < values.size(); ++i) {
        cout << values[i];
        if (i != values.size() - 1)
            cout << "/";
    }
}

static SearchEngine *_parse_astar(OptionParser &parser) {
    
    parser.add_option<ScalarEvaluator *>("eval", "evaluator for h-value");
    
    SearchEngine::add_options_to_parser(parser);
    Options opts = parser.parse();

    BFSSearch *engine = 0;
    if (!parser.dry_run()) { 
        engine = new BFSSearch(opts);
    }

    return engine;
}

static Plugin<SearchEngine> _plugin_astar("bfs", _parse_astar);
