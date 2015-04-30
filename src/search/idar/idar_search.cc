#include "idar_search.h"

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

IDARSearch::IDARSearch(
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
            if (!heuristics[i]->is_dead_end()) {
		max_h = max(max_h, heuristics[i]->get_heuristic());
	    } else {
		max_h = INT_MAX/2;
	    }
        }
        cout<<"max_h(constructor) = "<<max_h<<endl;
}

void IDARSearch::initialize() {
    //DFS
    ida_timer.reset();
    best_soln_sofar = INT_MAX/2;
    nodes_expanded_for_start_state = 0;
    nodes_generated_for_start_state = 0;
}

SearchStatus IDARSearch::step() {

        const GlobalState &initial_state = g_initial_state();

        int hmax = 0; 
        for (size_t i = 0; i < heuristics.size(); i++) {
            heuristics[i]->evaluate(initial_state);
	    if (!heuristics[i]->is_dead_end()) {
		hmax = max(hmax, heuristics[i]->get_heuristic());
	    } else {
		hmax = INT_MAX/2;
	    }
        }
        int h_initial = hmax;
        cout<<"h_initial = "<<h_initial<<endl;
        
        StateID initial_state_id = initial_state.get_id();
        Node node(initial_state_id, h_initial,  0, 0); //set the h_value, g_real and the level

        int total_d = 0;
        int total_expanded = 0;
	int total_generated = 0;

        int d = 0;
	d = idastar(node);
        if (d == INT_MAX) {
		//Solution not found
		cout<<"Solution NOT FOUND"<<endl;
		cout<<", expanded nodes = "<<nodes_expanded_for_start_state;
		cout<<", generated nodes = "<<nodes_generated_for_start_state<<"\n";
	} else {
                //Solution found
		cout<<"SOLUTION FOUND"<<endl;
		cout<<"\tcost = "<<d;
		cout<<", expanded nodes = "<<nodes_expanded_for_start_state;
		cout<<", generated nodes = "<<nodes_generated_for_start_state<<"\n";
	}
	total_d += d;
	total_expanded += nodes_expanded_for_start_state;
	total_generated += nodes_generated_for_start_state;


	cout<<"\n\tTotal depth: "<<total_d;
	cout<<", expansion: "<<total_expanded;
	cout<<", generation: "<<total_generated;
	cout<<"\n";



	return SOLVED;
}

int IDARSearch::idastar(Node node) {
	double bound, next_bound,  done;

        GlobalState global_state = g_state_registry->lookup_state(node.get_id());
	if (check_goal_and_set_plan(global_state)) {
		cout<<"line solution-found."<<endl;
                return 0; 
        }
	best_soln_sofar = INT_MAX;
	bound =  node.getHvalue();
        int count_bound = 1;
	while (1) {
		next_bound = INT_MAX;
		nodes_expanded_for_bound = 0;
		nodes_generated_for_bound = 0;
		done = dfs_heur(node, bound, next_bound, 0);
                cout<<"\t time_"<<count_bound<<" = "<<g_timer;
		cout<<", bound_"<<count_bound<<" = "<<bound;
		cout<<", nodes_expanded_for_bound = "<<nodes_expanded_for_bound;
		cout<<", nodes_generated_for_bound = "<<nodes_generated_for_bound;
		cout<<"\n";
		nodes_expanded_for_start_state += nodes_expanded_for_bound;
		nodes_generated_for_start_state += nodes_generated_for_bound;
		cout<<"done = "<<done<<endl;
		if (done) {
			cout<<"break the application because done = "<<done<<endl;
			break;
		}
		bound = next_bound;
		cout<<"the new bound is = "<<bound<<endl;
                cout<<"best_soln_sofar = "<<best_soln_sofar<<endl;
		if (best_soln_sofar <= bound) {
			cout<<"break the application because best_soln_sofar <= bound"<<endl;
			break;
		}
		count_bound++;
	}
	return best_soln_sofar;
}

int IDARSearch::dfs_heur(Node node, double bound, double &next_bound, double g_real) {
	nodes_expanded_for_bound++;
        cout<<"node expanded: h = "<<node.getHvalue()<<", g_real = "<<node.getGreal()<<", f = "<<node.getHvalue() + node.getGreal()<<"\n";
	std::vector<const GlobalOperator *> applicable_ops;
        GlobalState global_state = g_state_registry->lookup_state(node.get_id());
        g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
        cout<<"applicable_ops.size() = "<<applicable_ops.size()<<endl;
	cout<<"--------------childs-------------\n";
        for (size_t i = 0; i < applicable_ops.size(); i++) {
            nodes_generated_for_bound++;

	    const GlobalOperator *op = applicable_ops[i];
            GlobalState child =  g_state_registry->get_successor_state(global_state, *op);
	                
	    int hmax_value = 0; 
            for (size_t i = 0; i < heuristics.size(); i++) {
                heuristics[i]->evaluate(child);
		if (!heuristics[i]->is_dead_end()) {
			hmax_value = max(hmax_value, heuristics[i]->get_heuristic());
		} else {
			hmax_value = INT_MAX/2;
		}
            }

            int succ_h = hmax_value;
            search_progress.inc_generated();
            Node succ_node(child.get_id(), succ_h, g_real + get_adjusted_cost(*op), node.getLevel()+ 1);
	    cout<<"\tChild_"<<(i+1)<<" : h = "<<succ_h<<", g_real = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<"\n";

	    if (check_goal_and_set_plan(child)) {
		cout<<"\tSolution-found in dfs_heur."<<endl;
		if (best_soln_sofar > g_real + get_adjusted_cost(*op)) {
			best_soln_sofar = g_real + get_adjusted_cost(*op);
		}
		
		cout<<"\tbest_soln_sofar = "<<best_soln_sofar<<endl;
		if (best_soln_sofar <= bound) {
	           cout<<"\tbest_soln_sofar <= bound => return 1;"<<endl;
                   return  1;
		} else {
		   continue;
		}
            } else {
		cout<<"\t\tthe soluton WAS NOT found."<<endl;
		if (g_real + get_adjusted_cost(*op) + succ_h > bound) {
			if (next_bound > g_real + get_adjusted_cost(*op) + succ_h) {
				next_bound = g_real + get_adjusted_cost(*op) + succ_h;
			}

			cout<<"\t\tnext_bound = "<<next_bound<<endl;
		} else {
			cout<<"\t\tcall dfs again."<<endl;
			if (dfs_heur(succ_node, bound, next_bound, g_real + get_adjusted_cost(*op))) {
				cout<<"\t\tdfs_heur is executed again and return 1;"<<endl;
				return 1;
			} else {
				cout<<"\t\tdfs_heur is not returning true."<<endl;
			}
		}
	    }
	}
	cout<<"-------------end Childs-----------\n";
	cout<<"return 0;"<<endl;
	return 0;
}

void IDARSearch::print_heuristic_values(const vector<int> &values) const {
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

    IDARSearch *engine = 0;
    if (!parser.dry_run()) { 
        engine = new IDARSearch(opts);
    }

    return engine;
}

static Plugin<SearchEngine> _plugin_astar("idar", _parse_astar);
