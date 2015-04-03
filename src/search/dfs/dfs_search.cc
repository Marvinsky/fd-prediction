#include "dfs_search.h"

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

DFSSearch::DFSSearch(
    const Options &opts)
    : SearchEngine(opts),
      reopen_closed_nodes(opts.get<bool>("reopen_closed")),
      do_pathmax(opts.get<bool>("pathmax")),
      use_multi_path_dependence(opts.get<bool>("mpd")),
      open_list(opts.get<OpenList<StateID> *>("open")) {
    if (opts.contains("f_eval")) {
        f_evaluator = opts.get<ScalarEvaluator *>("f_eval");
    } else {
        f_evaluator = 0;
    }
    if (opts.contains("preferred")) {
        preferred_operator_heuristics =
            opts.get_list<Heuristic *>("preferred");
    } 
}

void DFSSearch::initialize() {
    //TODO children classes should output which kind of search
    cout << "Conducting best first search"
         << (reopen_closed_nodes ? " with" : " without")
         << " reopening closed nodes, (real) bound = " << bound
         << endl;
    if (do_pathmax)
        cout << "Using pathmax correction" << endl;
    if (use_multi_path_dependence)
        cout << "Using multi-path dependence (LM-A*)" << endl;
    assert(open_list != NULL);

    set<Heuristic *> hset;
    open_list->get_involved_heuristics(hset);

    for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); ++it) {
        estimate_heuristics.push_back(*it);
        search_progress.add_heuristic(*it);
    }

    // add heuristics that are used for preferred operators (in case they are
    // not also used in the open list)
    hset.insert(preferred_operator_heuristics.begin(),
                preferred_operator_heuristics.end());

    // add heuristics that are used in the f_evaluator. They are usually also
    // used in the open list and hence already be included, but we want to be
    // sure.
    if (f_evaluator) {
        f_evaluator->get_involved_heuristics(hset);
    }

    for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); ++it) {
        heuristics.push_back(*it);
    }

    assert(!heuristics.empty());

    //DFS
    ida_timer.reset();

    const GlobalState &initial_state = g_initial_state();
    for (size_t i = 0; i < heuristics.size(); ++i)
        heuristics[i]->evaluate(initial_state);
    open_list->evaluate(0, false);
    search_progress.inc_evaluated_states();
    search_progress.inc_evaluations(heuristics.size());

    if (open_list->is_dead_end()) {
        cout << "Initial state is a dead end." << endl;
    } else {
        search_progress.get_initial_h_values();
        if (f_evaluator) {
            f_evaluator->evaluate(0, false);
            search_progress.report_f_value(f_evaluator->get_value());
        }
        search_progress.check_h_progress(0);
        SearchNode node = search_space.get_node(initial_state);
        node.open_initial(heuristics[0]->get_value());

        open_list->insert(initial_state.get_id());
    }
}


void DFSSearch::statistics() const {
    search_progress.print_statistics();
    search_space.statistics();
}

SearchStatus DFSSearch::step() {

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
              //cout<<"new_f_value = "<<new_f_value<<endl;
	      //if (search_progress.updated_lastjump_f_value(new_f_value)) {
                 //search_progress.report_f_value(new_f_value);
              //}

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
 
        double nodes_total_generated = 0;
        for (map<Node2, double>::iterator iter = generated.begin(); iter != generated.end(); iter++) {
            
  	    double q = iter->second;
            nodes_total_generated += q;  
        }
        cout<<"Total of nodes generated: "<<nodes_total_generated<<endl;
        cout<<"generated.size() = "<<generated.size()<<endl;
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

pair<SearchNode, bool> DFSSearch::fetch_next_node() {
    /* TODO: The bulk of this code deals with multi-path dependence,
       which is a bit unfortunate since that is a special case that
       makes the common case look more complicated than it would need
       to be. We could refactor this by implementing multi-path
       dependence as a separate search algorithm that wraps the "usual"
       search algorithm and adds the extra processing in the desired
       places. I think this would lead to much cleaner code. */

    while (true) {
        if (open_list->empty()) {
            cout << "Completely explored state space -- no solution!" << endl;
            // HACK! HACK! we do this because SearchNode has no default/copy constructor
            SearchNode dummy_node = search_space.get_node(g_initial_state());
            return make_pair(dummy_node, false);
        }
        vector<int> last_key_removed;
        StateID id = open_list->remove_min(
            use_multi_path_dependence ? &last_key_removed : 0);
        // TODO is there a way we can avoid creating the state here and then
        //      recreate it outside of this function with node.get_state()?
        //      One way would be to store GlobalState objects inside SearchNodes
        //      instead of StateIDs
        GlobalState s = g_state_registry->lookup_state(id);
        SearchNode node = search_space.get_node(s);

        if (node.is_closed())
            continue;

        if (use_multi_path_dependence) {
            assert(last_key_removed.size() == 2);
            if (node.is_dead_end())
                continue;
            int pushed_h = last_key_removed[1];
            assert(node.get_h() >= pushed_h);
            if (node.get_h() > pushed_h) {
                // cout << "LM-A* skip h" << endl;
                continue;
            }
            assert(node.get_h() == pushed_h);
            if (!node.is_closed() && node.is_h_dirty()) {
                for (size_t i = 0; i < heuristics.size(); ++i)
                    heuristics[i]->evaluate(node.get_state());
                node.clear_h_dirty();
                search_progress.inc_evaluations(heuristics.size());

                open_list->evaluate(node.get_g(), false);
                bool dead_end = open_list->is_dead_end();
                if (dead_end) {
                    node.mark_as_dead_end();
                    search_progress.inc_dead_ends();
                    continue;
                }
                int new_h = heuristics[0]->get_heuristic();
                if (new_h > node.get_h()) {
                    assert(node.is_open());
                    node.increase_h(new_h);
                    open_list->insert(node.get_state_id());
                    continue;
                }
            }
        }

        node.close();
        assert(!node.is_dead_end());
        update_jump_statistic(node);
        search_progress.inc_expanded();
        return make_pair(node, true);
    }
}

void DFSSearch::reward_progress() {
    // Boost the "preferred operator" open lists somewhat whenever
    // one of the heuristics finds a state with a new best h value.
    open_list->boost_preferred();
}

void DFSSearch::dump_search_space() {
    search_space.dump();
}

void DFSSearch::update_jump_statistic(const SearchNode &node) {
    if (f_evaluator) {
        heuristics[0]->set_evaluator_value(node.get_h());
        //f_evaluator->evaluate(node.get_g(), false);
        //int new_f_value = f_evaluator->get_value();
	int new_f_value = node.get_g()+node.get_h();
	
        search_progress.report_f_value(new_f_value);
    }
}

void DFSSearch::print_heuristic_values(const vector<int> &values) const {
    for (size_t i = 0; i < values.size(); ++i) {
        cout << values[i];
        if (i != values.size() - 1)
            cout << "/";
    }
}

static SearchEngine *_parse(OptionParser &parser) {
    //open lists are currently registered with the parser on demand,
    //because for templated classes the usual method of registering
    //does not work:
    Plugin<OpenList<StateID> >::register_open_lists();

    parser.document_synopsis("Eager best first search", "");

    parser.add_option<OpenList<StateID> *>("open", "open list");
    parser.add_option<bool>("reopen_closed",
                            "reopen closed nodes", "false");
    parser.add_option<bool>("pathmax",
                            "use pathmax correction", "false");
    parser.add_option<ScalarEvaluator *>(
        "f_eval",
        "set evaluator for jump statistics. "
        "(Optional; if no evaluator is used, jump statistics will not be displayed.)",
        "",
        OptionFlags(false));
    parser.add_list_option<Heuristic *>
        ("preferred",
        "use preferred operators of these heuristics", "[]");
    SearchEngine::add_options_to_parser(parser);
    Options opts = parser.parse();

    DFSSearch *engine = 0;
    if (!parser.dry_run()) {
        opts.set<bool>("mpd", false);
        engine = new DFSSearch(opts);
    }

    return engine;
}

static SearchEngine *_parse_astar(OptionParser &parser) {
    parser.document_synopsis(
        "A* search (eager)",
        "A* is a special case of eager best first search that uses g+h "
        "as f-function. "
        "We break ties using the evaluator. Closed nodes are re-opened.");
    parser.document_note(
        "mpd option",
        "This option is currently only present for the A* algorithm and not "
        "for the more general eager search, "
        "because the current implementation of multi-path depedence "
        "does not support general open lists.");
    parser.document_note(
        "Equivalent statements using general eager search",
        "\n```\n--search astar(evaluator)\n```\n"
        "is equivalent to\n"
        "```\n--heuristic h=evaluator\n"
        "--search eager(tiebreaking([sum([g(), h]), h], unsafe_pruning=false),\n"
        "               reopen_closed=true, pathmax=false, progress_evaluator=sum([g(), h]))\n"
        "```\n", true);
    parser.add_option<ScalarEvaluator *>("eval", "evaluator for h-value");
    parser.add_option<bool>("pathmax",
                            "use pathmax correction", "false");
    parser.add_option<bool>("mpd",
                            "use multi-path dependence (LM-A*)", "false");
    SearchEngine::add_options_to_parser(parser);
    Options opts = parser.parse();

    DFSSearch *engine = 0;
    if (!parser.dry_run()) {
        GEvaluator *g = new GEvaluator();
        vector<ScalarEvaluator *> sum_evals;
        sum_evals.push_back(g);
        ScalarEvaluator *eval = opts.get<ScalarEvaluator *>("eval");
        sum_evals.push_back(eval);
        ScalarEvaluator *f_eval = new SumEvaluator(sum_evals);

        // use eval for tiebreaking
        std::vector<ScalarEvaluator *> evals;
        evals.push_back(f_eval);
        evals.push_back(eval);
        OpenList<StateID> *open = \
            new TieBreakingOpenList<StateID>(evals, false, false);

        opts.set("open", open);
        opts.set("f_eval", f_eval);
        opts.set("reopen_closed", true);
        engine = new DFSSearch(opts);
    }

    return engine;
}

static SearchEngine *_parse_greedy(OptionParser &parser) {
    parser.document_synopsis("Greedy search (eager)", "");
    parser.document_note(
        "Open list",
        "In most cases, eager greedy best first search uses "
        "an alternation open list with one queue for each evaluator. "
        "If preferred operator heuristics are used, it adds an extra queue "
        "for each of these evaluators that includes only the nodes that "
        "are generated with a preferred operator. "
        "If only one evaluator and no preferred operator heuristic is used, "
        "the search does not use an alternation open list but a "
        "standard open list with only one queue.");
    parser.document_note(
        "Closed nodes",
        "Closed node are not re-opened");
    parser.document_note(
        "Equivalent statements using general eager search",
        "\n```\n--heuristic h2=eval2\n"
        "--search eager_greedy([eval1, h2], preferred=h2, boost=100)\n```\n"
        "is equivalent to\n"
        "```\n--heuristic h1=eval1 --heuristic h2=eval2\n"
        "--search eager(alt([single(h1), single(h1, pref_only=true), single(h2), \n"
        "                    single(h2, pref_only=true)], boost=100),\n"
        "               preferred=h2)\n```\n"
        "------------------------------------------------------------\n"
        "```\n--search eager_greedy([eval1, eval2])\n```\n"
        "is equivalent to\n"
        "```\n--search eager(alt([single(eval1), single(eval2)]))\n```\n"
        "------------------------------------------------------------\n"
        "```\n--heuristic h1=eval1\n"
        "--search eager_greedy(h1, preferred=h1)\n```\n"
        "is equivalent to\n"
        "```\n--heuristic h1=eval1\n"
        "--search eager(alt([single(h1), single(h1, pref_only=true)]),\n"
        "               preferred=h1)\n```\n"
        "------------------------------------------------------------\n"
        "```\n--search eager_greedy(eval1)\n```\n"
        "is equivalent to\n"
        "```\n--search eager(single(eval1))\n```\n", true);

    parser.add_list_option<ScalarEvaluator *>("evals", "scalar evaluators");
    parser.add_list_option<Heuristic *>(
        "preferred",
        "use preferred operators of these heuristics", "[]");
    parser.add_option<int>(
        "boost",
        "boost value for preferred operator open lists", "0");
    SearchEngine::add_options_to_parser(parser);


    Options opts = parser.parse();
    opts.verify_list_non_empty<ScalarEvaluator *>("evals");

    DFSSearch *engine = 0;
    if (!parser.dry_run()) {
        vector<ScalarEvaluator *> evals =
            opts.get_list<ScalarEvaluator *>("evals");
        vector<Heuristic *> preferred_list =
            opts.get_list<Heuristic *>("preferred");
        OpenList<StateID> *open;
        if ((evals.size() == 1) && preferred_list.empty()) {
            open = new StandardScalarOpenList<StateID>(evals[0], false);
        } else {
            vector<OpenList<StateID> *> inner_lists;
            for (size_t i = 0; i < evals.size(); ++i) {
                inner_lists.push_back(
                    new StandardScalarOpenList<StateID>(evals[i], false));
                if (!preferred_list.empty()) {
                    inner_lists.push_back(
                        new StandardScalarOpenList<StateID>(evals[i], true));
                }
            }
            open = new AlternationOpenList<StateID>(
                inner_lists, opts.get<int>("boost"));
        }

        opts.set("open", open);
        opts.set("reopen_closed", false);
        opts.set("pathmax", false);
        opts.set("mpd", false);
        ScalarEvaluator *sep = 0;
        opts.set("f_eval", sep);
        opts.set("preferred", preferred_list);
        engine = new DFSSearch(opts);
    }
    return engine;
}

static Plugin<SearchEngine> _plugin("eager", _parse);
static Plugin<SearchEngine> _plugin_astar("dfs", _parse_astar);
static Plugin<SearchEngine> _plugin_greedy("eager_greedy", _parse_greedy);
