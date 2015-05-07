#include "idai_search.h"

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

IDAISearch::IDAISearch(
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

void IDAISearch::initialize() {
    //DFS
    ida_timer.reset();
    best_soln_sofar = INT_MAX/2;
    nodes_expanded_for_start_state = 0;
    nodes_generated_for_start_state = 0;
    found_solution = false;
}

SearchStatus IDAISearch::step() {

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
        Node node(initial_state_id, h_initial,  0, 0); //set the h_value, g_real and the level

        int total_d = 0;
        int total_expanded = 0;
	int total_generated = 0;

        int d = 0;
	d = idastar(node);
	cout<<"d = "<<d<<"\n";
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
	cout<<"total_d = "<<total_d<<"\n";
	total_expanded += nodes_expanded_for_start_state;
	total_generated += nodes_generated_for_start_state;


	cout<<"\n\tTotal depth: "<<total_d;
	cout<<", expansion: "<<total_expanded;
	cout<<", generation: "<<total_generated;
	cout<<"\n";



	return SOLVED;
}

int IDAISearch::idastar(Node node) {
	double bound, next_bound,  done;

        GlobalState global_state = g_state_registry->lookup_state(node.get_id());
	if (check_goal_and_set_plan(global_state)) {
		cout<<"solution_found_1"<<endl;
                return 0; 
        }
	best_soln_sofar = INT_MAX;
	bound =  node.getHvalue();
        int count_bound = 1;
	while (1) {
		next_bound = INT_MAX;
		nodes_expanded_for_bound = 0;
		nodes_generated_for_bound = 0;
		done = dfs_heur(node, bound, next_bound);
                cout<<"\t time_"<<count_bound<<": "<<g_timer;
		cout<<", bound_"<<count_bound<<": "<<bound;
		cout<<", nodes_expanded_for_bound: "<<nodes_expanded_for_bound;
		cout<<", nodes_generated_for_bound: "<<nodes_generated_for_bound;
		cout<<"\n";
		nodes_expanded_for_start_state += nodes_expanded_for_bound;
		nodes_generated_for_start_state += nodes_generated_for_bound;
		//cout<<"done = "<<done<<endl;
		if (done) {
			//cout<<"break the application because done = "<<done<<endl;
			break;
		}
		bound = next_bound;
		//cout<<"the new bound is = "<<bound<<endl;
                //cout<<"best_soln_sofar = "<<best_soln_sofar<<endl;
		if (best_soln_sofar <= bound) {
			//cout<<"break the application because best_soln_sofar <= bound"<<endl;
			break;
		}
		count_bound++;
	}
	return best_soln_sofar;
}

int IDAISearch::dfs_heur(Node node, double bound, double &next_bound) {
	std::stack<Node> buffer;
	buffer.push(node);

	while (!buffer.empty()) {
		nodes_expanded_for_bound++;
		Node ncp = buffer.top();
		double g_real = ncp.getGreal();
		buffer.pop();

		//cout<<"node expanded: h = "<<ncp.getHvalue()<<", g_real = "<<ncp.getGreal()<<", f = "<<ncp.getHvalue() + ncp.getGreal()<<"\n";
		std::vector<const GlobalOperator *> applicable_ops;
        	GlobalState global_state = g_state_registry->lookup_state(ncp.get_id());
        	g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
        	//cout<<"applicable_ops.size() = "<<applicable_ops.size()<<endl;
		//cout<<"--------------childs-------------\n";
		//cout<<"before checking.."<<endl;
		if (check_goal_and_set_plan(global_state)) {
			//cout<<"\tsolution_found_4"<<endl;
			//cout<<"best_soln_sofar = "<<best_soln_sofar<<"\n";
			//cout<<"g_real = "<<g_real<<"\n";
			//cout<<"h_value = "<<ncp.getHvalue()<<"\n";
			if (best_soln_sofar > g_real + ncp.getHvalue()) {
				best_soln_sofar = g_real + ncp.getHvalue();
			}
			//cout<<"\tbest_soln_sofar = "<<best_soln_sofar<<endl;
			/*if (best_soln_sofar <= bound) {
	           		//cout<<"\tbest_soln_sofar <= bound => return 1;"<<endl;
                   		return  1;
			} else {
		   		continue;
			}*/
			return 1;
            	}
		//cout<<"after checking.."<<endl;
		L.clear();
		check.clear();
        	for (size_t i = 0; i < applicable_ops.size(); i++) {
			nodes_generated_for_bound++;
			
			const GlobalOperator *op = applicable_ops[i];
            		GlobalState child =  g_state_registry->get_successor_state(global_state, *op);
	                
	    		int hmax_value = 0;
			int cost_op = get_adjusted_cost(*op); 
            		for (size_t i = 0; i < heuristics.size(); i++) {
                		heuristics[i]->evaluate(child);
				if (!heuristics[i]->is_dead_end()) {
					hmax_value = max(hmax_value, heuristics[i]->get_heuristic());
				} else {
					hmax_value = INT_MAX/2;
				}
            		}
            		int succ_h = hmax_value;
			Node succ_node(child.get_id(), succ_h, g_real + cost_op, node.getLevel()+ 1);
	    		//cout<<"\tChild_"<<(i+1)<<" : h = "<<succ_h<<", g_real = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<"\n";

			if (cost_op == 0) {
				BFS(succ_node);
				if (found_solution) {
					return 1;
				}

				std::set<Node, classcomp>::iterator iter;
				double new_g_real = 0, new_h_value = 0;
				for (iter = L.begin(); iter != L.end(); ++iter) {
					Node n = *iter;
					StateID new_state_id = n.get_id();
					new_g_real = n.getGreal();
					new_h_value = n.getHvalue();
		
					//cout<<"\t\tNode from BFS: h "<<new_h_value<<", g = "<<new_g_real<<", f = "<<new_h_value + new_g_real<<", level = "<<n.getLevel()<<", stateID = "<<new_state_id<<"\n";

					GlobalState global_state2 = g_state_registry->lookup_state(new_state_id);

					/*if (check_goal_and_set_plan(global_state2)) {
						cout<<"\tsolution_found_2"<<endl;
						cout<<"best_soln_sofar = "<<best_soln_sofar<<"\n";
						cout<<"new_g_real = "<<new_g_real<<"\n";
						cout<<"new_h_value = "<<new_h_value<<"\n";
						if (best_soln_sofar > bound) {
							best_soln_sofar = bound;
						}
						
						return 1;
            				} else {*/
						//cout<<"\t\tthe soluton WAS NOT found."<<endl;
						if (new_g_real + new_h_value > bound) {
							if (next_bound > new_g_real + new_h_value) {
								next_bound = new_g_real + new_h_value;
							}
							//cout<<"\t\tnext_bound = "<<next_bound<<endl;
						} else {
							buffer.push(n);
						}
	    				//}//End check goal
				} //End for set L
			} else {
				/*if (check_goal_and_set_plan(child)) {
					cout<<"\tsolution_found_3"<<endl;
					cout<<"best_soln_sofar = "<<best_soln_sofar<<"\n";
					cout<<"g_real = "<<g_real<<"\n";
					cout<<"cost_op = "<<cost_op<<"\n";
					cout<<"succ_h = "<<succ_h<<"\n";
					if (best_soln_sofar > g_real + cost_op) {
						best_soln_sofar = g_real + cost_op;
					}
					cout<<"\tbest_soln_sofar = "<<best_soln_sofar<<endl;

					return 1;
            			} else {*/
					//cout<<"\t\tthe soluton WAS NOT found."<<endl;
					if (g_real + cost_op + succ_h > bound) {
						if (next_bound > g_real + cost_op + succ_h) {
							next_bound = g_real + cost_op + succ_h;
						}
						//cout<<"\t\tnext_bound = "<<next_bound<<endl;
					} else {
						buffer.push(succ_node);
					}
	    			//}//End check goal
			}
		}//End for applicable ops
		//cout<<"-------------end Childs-----------\n";
	}
	//cout<<"return 0;"<<endl;
	return 0;
}

void IDAISearch::BFS(Node root) {
	std::queue<Node> expand;
	check.insert(root);
	expand.push(root);
	std::vector<const GlobalOperator *> applicable_ops;

	int g_real = 0, level = 0;
	//int counter = 0;
	while (!expand.empty()) {
		nodes_expanded_for_bound++; 
		Node nodecp = expand.front();
		g_real = nodecp.getGreal();
                StateID state_id = nodecp.get_id();
                level = nodecp.getLevel();
                //cout<<"\t\tBFSNode expanded: h = "<<nodecp.getHvalue()<<", g_real = "<<nodecp.getGreal()<<", f = "<<nodecp.getHvalue() + nodecp.getGreal()<<", level = "<<level<<", stateID,: "<<state_id<<"\n";
		expand.pop();
		//counter++;
		applicable_ops.clear();
                //Recover the global_state
                GlobalState global_state = g_state_registry->lookup_state(state_id);
                g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
		//cout<<"checing goal in BFS."<<endl;
		if (check_goal_and_set_plan(global_state)) {
			//cout<<"\tsolution_found_5"<<endl;
			//cout<<"best_soln_sofar = "<<best_soln_sofar<<"\n";
			//cout<<"g_real = "<<g_real<<"\n";
			//cout<<"h_value = "<<nodecp.getHvalue()<<"\n";
			if (best_soln_sofar > g_real + nodecp.getHvalue()) {
				best_soln_sofar = g_real + nodecp.getHvalue();
			}
			//cout<<"\tbest_soln_sofar = "<<best_soln_sofar<<endl;
			/*if (best_soln_sofar <= bound) {
	           		//cout<<"\tbest_soln_sofar <= bound => return 1;"<<endl;
                   		return  1;
			} else {
		   		continue;
			}*/
			//return 1;
			found_solution = true;
            	}

		if (found_solution) {
			break;
		}

                //cout<<"\t\tBFSExpanded state:";
                //cout<<": state_id:"<<global_state.get_id()<<":";
                //global_state.dump_inline();
                //fflush(NULL);

                //cout<<"\n\t\tBFS: applicable_ops.size() = "<<applicable_ops.size()<<endl;
                //cout<<"\t\t--------------childs-------------\n";

		for (size_t i = 0; i < applicable_ops.size(); i++) {
                        const GlobalOperator *op = applicable_ops[i];
                        GlobalState child =  g_state_registry->get_successor_state(global_state, *op);
                        //cout<<"\n\t\t\tBFSChild_"<<(i+1);
                        //cout<<": state_id:"<<child.get_id()<<":";
                        //child.dump_inline();
                        //fflush(NULL);

                        int hmax_value = 0;
			int cost_op = get_adjusted_cost(*op);

                        Node succ_node(child.get_id(), hmax_value, g_real + cost_op, level + 1);
                        std::pair<std::set<Node, classcomp>::iterator, bool> p;
                        p = check.insert(succ_node);
			if (p.second) {
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
                                succ_node.setHvalue(succ_h);				

                                //cout<<"\t\t\tNew node : h = "<<succ_h<<", g_real = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<", level = "<<succ_node.getLevel()<<", stateID,: "<<child.get_id()<<"\n";

				if (cost_op == 0) {
					expand.push(succ_node);
				} else {
					L.insert(succ_node);
				}
			} else {
				//cout<<"\t\t\tnode with stateID = "<<child.get_id()<<"\n";
			}
		}
	}//End while
	//nodes_expanded_for_bound += counter;
}


void IDAISearch::print_heuristic_values(const vector<int> &values) const {
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

    IDAISearch *engine = 0;
    if (!parser.dry_run()) { 
        engine = new IDAISearch(opts);
    }

    return engine;
}

static Plugin<SearchEngine> _plugin_astar("idai", _parse_astar);
