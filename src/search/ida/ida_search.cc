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
            int aux = heuristics[i]->get_heuristic();
            if (max_h < aux) {
                max_h = aux;
            }
        }
        cout<<"max_h(constructor) = "<<max_h<<endl;
}

void IDASearch::initialize() {
    //DFS
    ida_timer.reset();

}

SearchStatus IDASearch::step() {

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
        queue.push(node);
        double count_nodes = 0.0;

        while (!queue.empty()) {
              SSNode nodecp = queue.top();
              int g_real = nodecp.getLevel();
              int level = nodecp.getGreal();
              //cout<<"Raiz: h = "<<nodecp.getHvalue()<<", g_real = "<<g_real<<", f = "<<nodecp.getHvalue() + g_real<<"\n";
	      queue.pop();

              StateID state_id = nodecp.get_id();

              Node2 node2(nodecp.getHvalue() + g_real, level);
              int new_f_value = nodecp.getHvalue() + g_real;

              
	      if (search_progress.updated_lastjump_f_value_sscc(new_f_value)) {

                 bool flag = search_progress.showReportLastjump(new_f_value);
                 if (flag)  {
                        cout<<"f_boundary = "<<new_f_value;
                        generateExpandedReport(false);                 
                        generateGeneratedReport(false);
                 }
              }

	      //count nodes expanded
              if (new_f_value <= depth) {
		std::pair<std::map<Node2, double>::iterator, bool> ret0;
              	std::map<Node2, double>::iterator it0;

              	ret0 = expanded.insert(pair<Node2, double>(node2, 1));
              	it0 = ret0.first;

              	if (ret0.second) {

	      	} else {
		  it0->second += 1;
	      	} 

              	count_nodes++;
              }
           
              

              std::vector<const GlobalOperator *> applicable_ops;

              GlobalState global_state = g_state_registry->lookup_state(nodecp.get_id());                
              g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
              


              double amount = (double)applicable_ops.size();

              //count nodes generated
	      std::pair<std::map<Node2, double>::iterator, bool> ret;

              std::map<Node2, double>::iterator it;
              
              ret = generated.insert(std::pair<Node2, double>(node2, amount));
              it = ret.first;

              if (ret.second) {
                 //cout<<"new node is added."<<endl;
              }  else {
                 //cout<<"node is updated."<<endl;
                 it->second += amount;
                 //cout<<"new = "<<it->second<<endl;
              }
               
	      /*if (search_progress.updated_lastjump_f_value_sscc(new_f_value)) {

                 bool flag = search_progress.showReportLastjump(new_f_value);
                 if (flag)  {
                        cout<<"f_boundary = "<<new_f_value;
                        generateExpandedReport(false);
                        generateGeneratedReport(false);
                 }
              }*/

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
                  //cout<<"\tNode generated: h = "<<succ_h<<", g = "<<succ_node.get_g_value()<<", f = "<<succ_h + succ_node.get_g_value()<<"\n";

                  if (succ_h + succ_node.getGreal() <= depth) {
                     //cout<<"\tNode generated: h = "<<succ_h<<", g = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<"\n";
                     queue.push(succ_node);
                  } else {
                     //cout<<"\tpruned!\n";
                  }
              }
        }

        double ida_timer_value = ida_timer();
        cout<<"ida_timer: "<<ida_timer()<<endl; 
        cout<<"end expansion of nodes finished."<<endl;
        cout<<"Total of nodes expanded by counter marvinsky: "<<count_nodes<<endl;
       
        double nodes_total_expanded = 0;
	
	for (map<Node2, double>::iterator iter = expanded.begin(); iter != expanded.end(); iter++) {
            double q = iter->second;
            nodes_total_expanded += q;
        }
        cout<<"Total of nodes expanded: "<<nodes_total_expanded<<endl;
 
        generateGeneratedReport(true); 
        generated.clear();

        ofstream output;
    	string dominio = domain_name;
    	string tarefa =  problem_name2;
    	string heuristica = heuristic_name2;
    	cout<<"changing the code."<<endl;
    	cout<<"dominio = "<<dominio<<endl;
    	cout<<"tarefa = "<<tarefa<<endl;
    	cout<<"heuristica = "<<heuristica<<endl;

    	string directoryDomain = "mkdir /home/marvin/marvin/testdfs/"+heuristica+"/reportdfs/"+dominio;
    	if (system(directoryDomain.c_str())) {
           cout<<dominio.c_str()<<" directory created."<<endl;
        }

    	string directoryFdist = "mkdir /home/marvin/marvin/testdfs/"+heuristica+"/reportdfs/"+dominio+"/fdist/";
    	if (system(directoryFdist.c_str())) {
           cout<<"fdist created."<<endl;
        }

    	string outputFile = "/home/marvin/marvin/testdfs/"+heuristica+"/reportdfs/"+dominio+"/fdist/"+tarefa;
    	cout<<"outputFile = "<<outputFile.c_str()<<endl;
    	output.open(outputFile.c_str());
    	output<<"\t"<<outputFile.c_str()<<"\n";
    	output<<"nodes_expanded: "<<nodes_total_expanded<<"\n";
    	output<<"ida_timer: "<<ida_timer_value<<"\n\n";

    	for (int i = 0; i <= depth; i++) {
       		int k = 0;
       		vector<int> f;
       		vector<double> q;
       		for (map<Node2, double>::iterator iter = expanded.begin(); iter != expanded.end(); iter++) {
           		Node2 n = iter->first;
           		if (i == n.getL()) {
              		    k++;
              		    f.push_back(n.getF());
              		    q.push_back(iter->second);
           	        }      
                }
       		cout<<"g:"<<i<<"\n";
       		output<<"g:"<<i<<"\n";

       		cout<<"size: "<<k<<"\n";
       		output<<"size: "<<k<<"\n";

       		for (size_t j = 0; j < f.size(); j++) {
           		cout<<"\tf: "<<f.at(j)<<"\tq: "<<q.at(j)<<"\n";
           		output<<"\tf: "<<f.at(j)<<"\tq: "<<q.at(j)<<"\n";
       		}
       		output<<"\n";
       		cout<<"\n";
    	}
    	output.close();
        expanded.clear(); 
   	return SOLVED;
}

void IDASearch::generateGeneratedReport(bool flag) {
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

void IDASearch::generateExpandedReport(bool flag) {
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
