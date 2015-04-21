#include "ida_search.h"

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

IDASearch::IDASearch(
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
		break;
	    }
        }

        cout<<"max_h(constructor) = "<<max_h<<endl;
}

void IDASearch::initialize() {
    //DFS
    ida_timer.reset();
    best_soln_sofar = INT_MAX/2;
    nodes_expanded_for_start_state = 0;
    nodes_generated_for_start_state = 0;
    next_bound = 0;
    SOLUTION_FOUND = false; 
}

SearchStatus IDASearch::step() {

        const GlobalState &initial_state = g_initial_state();

        int hmax = 0; 
        for (size_t i = 0; i < heuristics.size(); i++) {
            heuristics[i]->evaluate(initial_state);
	    if (!heuristics[i]->is_dead_end()) {
		hmax = max(hmax, heuristics[i]->get_heuristic());
	    } else {
		hmax = INT_MAX/2;
		break;
	    }
        }

        int h_initial = hmax;
        cout<<"h_initial = "<<h_initial<<endl;
        
        StateID initial_state_id = initial_state.get_id();
        SSNode node(initial_state_id, h_initial,  0, 0); //set the h_value, g_real and the level

        int total_d = 0;
        int total_expanded = 0;
	int total_generated = 0;

        int d;
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

int IDASearch::idastar(SSNode node) {
	double bound;

        GlobalState global_state = g_state_registry->lookup_state(node.get_id());
	if (check_goal_and_set_plan(global_state)) {
		cout<<"line 90 solution-found."<<endl;
                return 0; 
        }
	best_soln_sofar = INT_MAX;
	bound =  node.getHvalue();
        int count_bound = 1;
	
	while (!SOLUTION_FOUND) {
		next_bound = INT_MAX;
		nodes_expanded_for_bound = 0;
		nodes_generated_for_bound = 0;
		next_bound = dfs_heur(node, bound, next_bound);
                cout<<"\t time_"<<count_bound<<" = "<<g_timer;
		cout<<", bound_"<<count_bound<<" = "<<bound;
		cout<<", nodes_expanded_for_bound = "<<nodes_expanded_for_bound;
		cout<<", nodes_generated_for_bound = "<<nodes_generated_for_bound;
		cout<<"\n";
		nodes_expanded_for_start_state += nodes_expanded_for_bound;
		nodes_generated_for_start_state += nodes_generated_for_bound;
		bound = next_bound;	
		count_bound++;
		//g_state_registry = new StateRegistry;
	}
	return best_soln_sofar;
}

int IDASearch::dfs_heur(SSNode node, double bound, double next_bound) {
	stack<SSNode> queue;
	queue.push(node);
	cout<<"\n\t\tbound = "<<bound;
	cout<<"\t\tnext_bound = "<<next_bound<<"\n";

	double g_real = 0, h_value = 0, level = 0;
	//GlobalState global_state;
	std::vector<const GlobalOperator *> applicable_ops;

	while (!queue.empty()) {
		nodes_expanded_for_bound++;
		SSNode nodecp = queue.top();
		g_real = nodecp.getGreal();
		h_value = nodecp.getHvalue();
		level = nodecp.getLevel();
		queue.pop();
		
		StateID state_id = nodecp.get_id();
				
		applicable_ops.clear();
        	GlobalState global_state = g_state_registry->lookup_state(state_id);
        	g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
		
		cout<<"Node expanded: h = "<<h_value<<", g_real = "<<g_real<<", f = "<<h_value + g_real<<", level = "<<level<<", StateID = "<<state_id<<"\n";

		cout<<"applicable_ops.size() = "<<applicable_ops.size()<<endl;
		cout<<"------------------Child----------------\n\n";
		
		for (size_t i = 0; i < applicable_ops.size(); ++i) {
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
					break;
				}
                  	}

			int succ_h = hmax_value;
			int cost_op = get_adjusted_cost(*op);
            		search_progress.inc_generated();
            		
			SSNode succ_node(child.get_id(), succ_h, g_real + cost_op, level+ 1);
            		cout<<"\t\tChild_"<<(i+1)<<" : h = "<<succ_h<<", g_real = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<", StateID = "<<child.get_id()<<"\n";
			cout<<"\t\tChild_state:(";
                	child.dump_inline();
                	cout<<"),state_id:"<<child.get_id()<<"\n";
			fflush(NULL);

			if (cost_op == 0) {
				cout<<"\tget_adjusted_cost(*op) == 0\n";
                		buffer = BFS(succ_node);
                		std::set<SSNode>::iterator it;
			
				double new_g_real = 0, new_level = 0;
				std::vector<const GlobalOperator *> applicable_ops2;
				for (it = buffer.begin(); it != buffer.end(); it++) {
                        		SSNode ncp = *it;
                        		StateID new_state_id = ncp.get_id();
                        		new_g_real = ncp.getGreal();
                        		new_level = ncp.getLevel();

                        		cout<<"\t\tExpanded node that comes from BFS: h = "<<ncp.getHvalue()<<", g_real = "<<ncp.getGreal()<<", f  = "<<ncp.getHvalue() + ncp.getGreal()<<", level = "<<ncp.getLevel()<<", StateID = "<<new_state_id<<"\n";

                        		applicable_ops2.clear();
                        		GlobalState global_state2 = g_state_registry->lookup_state(new_state_id);
                        		g_successor_generator->generate_applicable_ops(global_state2, applicable_ops2);
                        		cout<<"\t\tExpanded node that comes from BFS:(";
                			global_state2.dump_inline();
                			cout<<"),state_id:"<<global_state2.get_id()<<"\n";
					fflush(NULL);

					cout<<"\t\tapplicable_ops.size() = "<<applicable_ops2.size()<<endl;
					cout<<"\t\t--------------childs-------------\n";
					
                        		for (size_t j = 0; j < applicable_ops2.size(); j++) {
                                		const GlobalOperator *op2 = applicable_ops2[j];
                                		GlobalState new_child =  g_state_registry->get_successor_state(global_state2, *op2);

                                		int hmax_value2 = 0;
                                		for (size_t i = 0; i < heuristics.size(); i++) {
                                        		heuristics[i]->evaluate(new_child);
                                        		if (!heuristics[i]->is_dead_end()) {
								hmax_value2 = max(hmax_value2, heuristics[i]->get_heuristic());
							} else {
								hmax_value2 = INT_MAX/2;
								break;
							}
                                		}
                                		int succ_h2 = hmax_value2;
						int cost_op2 = get_adjusted_cost(*op2);

                                		SSNode succ_node2(new_child.get_id(), succ_h2, new_g_real + cost_op2, new_level + 1);

						if (check_goal_and_set_plan(new_child)) {
							cout<<"\t\tSolution Found!"<<endl;
							SOLUTION_FOUND = true;
		 					if (best_soln_sofar > succ_h2 + succ_node2.getGreal()) {
								best_soln_sofar = succ_h2 + succ_node2.getGreal(); 
							}
							return best_soln_sofar;
						} else {
							cout<<"\t\tThe solution was not found."<<endl;
							if (new_g_real + cost_op2 + succ_h2 <= bound) {
								cout<<"\t\tInserting to the queue f <= bound"<<endl;
								queue.push(succ_node2);
							} else {
								cout<<"\t\tFinding the next_bound"<<endl;
								if (next_bound >  new_g_real + cost_op2 + succ_h2) {
									next_bound =  new_g_real + cost_op2 + succ_h2;
									cout<<"\t\tnext_bound = "<<next_bound<<endl;
								}
							}
						} //End check goal
					} //End for applicable_ops2
				} //End set buffer
				
				buffer.clear();
			} else { //else cost_op != 0
				if (check_goal_and_set_plan(child)) {
					cout<<"\t\tSolution Found!"<<endl;
					SOLUTION_FOUND = true;
		 			if (best_soln_sofar > succ_h + succ_node.getGreal()) {
						best_soln_sofar = succ_h + succ_node.getGreal(); 
					}
					return best_soln_sofar;
				} else {
					cout<<"\t\tThe solution was not found."<<endl;
					if (g_real + cost_op + succ_h <= bound) {
						cout<<"\t\tInserting to the queue f <= bound"<<endl;
						queue.push(succ_node);
					} else {
						cout<<"\t\tFinding the next_bound"<<endl;
						if (next_bound >  g_real + cost_op + succ_h) {
							next_bound =  g_real + cost_op + succ_h;
							cout<<"\t\tnext_bound = "<<next_bound<<endl;
						}
					}
				}
			}
			cout<<"End Child_"<<(i+1)<<"\n";
		}//End for
		cout<<"-----------------End Childs------------------\n\n";
	}//End while

	cout<<"return_next_bound = "<<next_bound<<"\n";
	return next_bound;	
}

set<SSNode, classcomp> IDASearch::BFS(SSNode root) {
        std::set<SSNode, classcomp> L;
	std::queue<SSNode> empty;
	std::swap(D, empty);
        D.push(root);
        cout<<"\t\t\tD.size() = "<<D.size()<<endl;
	int g_real = 0, level = 0;
	std::vector<const GlobalOperator *> applicable_ops;	
        while (!D.empty()) {
                cout<<"\t\t\tD.size() = "<<D.size()<<endl;
                SSNode nodecp = D.front();
                g_real = nodecp.getGreal();
                StateID state_id = nodecp.get_id();
                level = nodecp.getLevel();
                

		cout<<"\t\t\tBFSNode expanded: h = "<<nodecp.getHvalue()<<", g_real = "<<nodecp.getGreal()<<", f = "<<nodecp.getHvalue() + nodecp.getGreal()<<", level = "<<level<<", stateID,: "<<state_id<<"\n";

                D.pop();
		applicable_ops.clear();
                //Recover the global_state
                GlobalState global_state = g_state_registry->lookup_state(state_id);
                g_successor_generator->generate_applicable_ops(global_state, applicable_ops);

                cout<<"\nBFSExpanded state:(";
                global_state.dump_inline();
                cout<<"),state_id:"<<global_state.get_id()<<"\n";
		fflush(NULL);


                cout<<"\t\t\tBFS: applicable_ops.size() = "<<applicable_ops.size()<<endl;
                cout<<"\t\t\t--------------childs-------------\n";
                for (size_t i = 0; i < applicable_ops.size(); i++) {
			const GlobalOperator *op = applicable_ops[i];
                        GlobalState child =  g_state_registry->get_successor_state(global_state, *op);

                        int hmax_value = 0;
                        for (size_t i = 0; i < heuristics.size(); i++) {
                                heuristics[i]->evaluate(child);
				if (!heuristics[i]->is_dead_end()) {
			 		hmax_value = max(hmax_value, heuristics[i]->get_heuristic());	
				} else {
					hmax_value = INT_MAX/2;
					break;
				}
                        }
                        int succ_h = hmax_value;
			int cost_op = get_adjusted_cost(*op);
                        SSNode succ_node(child.get_id(), succ_h, g_real + cost_op, level + 1);

                        cout<<"\t\t\tBFSChild state:(";
                        child.dump_inline();
                        cout<<"),state_id:"<<child.get_id()<<"\n";
			fflush(NULL);

                        cout<<"\t\t\tBFSChild_"<<(i+1)<<" : h = "<<succ_h<<", g_real = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<", level = "<<succ_node.getLevel()<<", stateID,: "<<child.get_id()<<"\n";

			if (cost_op == 0) {
                                cout<<"\t\t\tcost = 0\n";
                                D.push(succ_node);
                        } else {
                                cout<<"\t\t\tcost != 0\n";
                                L.insert(succ_node);
                        }
                }
                cout<<"\t\t\t-------------End childs------------\n";
        }
        cout<<"\t\t\tD.empty() == "<<D.empty()<<endl;
        return L;
}


void IDASearch::print_heuristic_values(const vector<int> &values) const {
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

    IDASearch *engine = 0;
    if (!parser.dry_run()) { 
        engine = new IDASearch(opts);
    }

    return engine;
}

static Plugin<SearchEngine> _plugin_astar("ida", _parse_astar);
