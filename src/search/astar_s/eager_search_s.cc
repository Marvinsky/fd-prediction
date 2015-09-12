#include "eager_search_s.h"

#include "../globals.h"
#include "../heuristic.h"
#include "../option_parser.h"
#include "../successor_generator.h"
#include "../g_evaluator.h"
#include "../sum_evaluator.h"
#include "../max_evaluator.h"
#include "../plugin.h"

#include <cassert>
#include <cstdlib>
#include <set>


//file
#include <iostream>
#include <fstream>

//ss+culprits
#include "../ext/boost/dynamic_bitset.hpp"
#include <boost/lexical_cast.hpp>
#include <stdlib.h>

using namespace std;

EagerSearchS::EagerSearchS(
    const Options &opts)
    : SearchEngine(opts),
      reopen_closed_nodes(opts.get<bool>("reopen_closed")),
      do_pathmax(opts.get<bool>("pathmax")),
      use_multi_path_dependence(opts.get<bool>("mpd")),
      open_list(opts.get<OpenList<StateID> *>("open")) {
	f_evaluator = opts.get<ScalarEvaluator *>("f_eval");
    /*if (opts.contains("f_eval")) {
        f_evaluator = opts.get<ScalarEvaluator *>("f_eval");
    } else {
        f_evaluator = 0;
    }*/
    if (opts.contains("preferred")) {
        preferred_operator_heuristics =
            opts.get_list<Heuristic *>("preferred");
    }
}

void EagerSearchS::initialize() {
    if(use_saved_pdbs){
      stored_GA_patterns.clear();
      cout<<"cleared stored_GA_patterns"<<endl;
    }
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
	
    for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); it++) {
	  //Eliminate any heuristics which were not generated because we ran out of time
	  //currently this is hacked to return not using heuristics
      if((*it)->is_using()){
	heuristics.push_back(*it);
        estimate_heuristics.push_back(*it);
        search_progress.add_heuristic(*it);
      }
    }
    cout<<"initial heuristics size:"<<heuristics.size()<<endl;

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
    else{
      cout<<"f_evaluator is false"<<endl;
      exit(0);
    }

    const GlobalState &initial_state = g_initial_state();
    cout<<"# heuristics before eliminating those not supporting conditional effects:"<<heuristics.size()<<endl;

    orig_heuristics=heuristics;
    //remove any heuristics which have been set to stop_using
    heuristics.clear();
    for (size_t i = 0; i < orig_heuristics.size(); i++){
	    if(!heuristics[i]->is_using()){
	      //cout<<"removing heur";heuristics[i]->print_heur_name();cout<<endl;
	      //heuristics[i]->free_up_memory(search_space);
	      //string heur_name=heuristics[i]->get_heur_name();
	      //if(heur_name.substr(0,13)=="heur:,lp_pdb,"){
		//exit(0);
	      //}
	      continue;
	    }
	    else{
	      heuristics[i]->evaluate(initial_state);
	      heuristics.push_back(orig_heuristics.at(i));
	    }
    }
    cout<<"# heuristics after eliminating those not supporting conditional effects:"<<heuristics.size()<<endl;
    heuristics=orig_heuristics;
    
    assert(!heuristics.empty());
    

    cout<<"active heuristics size:"<<heuristics.size()<<endl;
    
    int max_h=0;

    bool dead_end=false;
    for (size_t i = 0; i < heuristics.size(); i++){
	heuristics[i]->evaluate(initial_state);
	dead_end=heuristics[i]->is_dead_end();
	if(dead_end){
	  //cout<<"dead end found"<<endl;fflush(NULL);
	  break;
	}
	max_h = max(max_h,heuristics[i]->get_heuristic());
    }
    cout<<"Initial max_h:"<<max_h<<endl;
    open_list->evaluate2(0, max_h);

    //open_list->evaluate(0, false);
    search_progress.inc_evaluated_states();
    search_progress.inc_evaluations(heuristics.size());

    if (dead_end) {
        cout << "Initial state is a dead end." << endl;
    } else {
        search_progress.get_initial_h_values();
        if (f_evaluator) {
            f_evaluator->evaluate(0, false);
            search_progress.report_f_value(max_h);
        }
        search_progress.check_h_progress(0);
        SearchNode node = search_space.get_node(initial_state);
        node.open_initial(heuristics[0]->get_value());

        open_list->insert(initial_state.get_id());
    }
}


void EagerSearchS::statistics() const {
    search_progress.print_statistics();
    search_space.statistics();
}

SearchStatus EagerSearchS::step() {
    pair<SearchNode, bool> n = fetch_next_node();
    if (!n.second) {
        return FAILED;
    }
    SearchNode node = n.first;

    GlobalState s = node.get_state();
    if (check_goal_and_set_plan(s))
        return SOLVED;

    vector<const GlobalOperator *> applicable_ops;
    set<const GlobalOperator *> preferred_ops;

    g_successor_generator->generate_applicable_ops(s, applicable_ops);
    // This evaluates the expanded state (again) to get preferred ops
    for (size_t i = 0; i < preferred_operator_heuristics.size(); ++i) {
        Heuristic *h = preferred_operator_heuristics[i];
        h->evaluate(s);
        if (!h->is_dead_end()) {
            // In an alternation search with unreliable heuristics, it is
            // possible that this heuristic considers the state a dead end.
            vector<const GlobalOperator *> preferred;
            h->get_preferred_operators(preferred);
            preferred_ops.insert(preferred.begin(), preferred.end());
        }
    }
    search_progress.inc_evaluations(preferred_operator_heuristics.size());

    for (size_t i = 0; i < applicable_ops.size(); ++i) {
      //cout<<"op:"<<i<<"out of"<<applicable_ops.size()<<endl;fflush(stdout);
        const GlobalOperator *op = applicable_ops[i];

        if ((node.get_real_g() + op->get_cost()) >= bound)
            continue;

        GlobalState succ_state = g_state_registry->get_successor_state(s, *op);
	//cout<<"Succ_state:";succ_state.dump_inline();fflush(NULL);
        search_progress.inc_generated();
        bool is_preferred = (preferred_ops.find(op) != preferred_ops.end());

        SearchNode succ_node = search_space.get_node(succ_state);

        // Previously encountered dead end. Don't re-evaluate.
        if (succ_node.is_dead_end()){
	  //cout<<"dead end found"<<endl;fflush(NULL);
            continue;
	}

        // update new path
        if (use_multi_path_dependence || succ_node.is_new()) {
            bool h_is_dirty = false;
            for (size_t j = 0; j < heuristics.size(); ++j) {
                /*
                  Note that we can't break out of the loop when
                  h_is_dirty is set to true or use short-circuit
                  evaluation here. We must call reach_state for each
                  heuristic for its side effects.
                */
                if (heuristics[j]->reach_state(s, *op, succ_state))
                    h_is_dirty = true;
            }
            if (h_is_dirty && use_multi_path_dependence)
                succ_node.set_h_dirty();
        }

	int succ_h=0;
	bool dead_end=false;
        if (succ_node.is_new()) {
	  //cout<<"succ_node is new"<<endl;fflush(stdout);
            // We have not seen this state before.
            // Evaluate and create a new node.
            for (size_t j = 0; j < heuristics.size(); ++j){
                heuristics[j]->evaluate(succ_state);
		dead_end=heuristics[j]->is_dead_end();
	      if (dead_end) {
		//cout<<"dead end found"<<endl;fflush(NULL);
		  succ_node.mark_as_dead_end();
		  search_progress.inc_dead_ends();
		  succ_h = INT_MAX/2;
		  break;
	      }
	      else{
		succ_h = max(succ_h,heuristics[j]->get_heuristic());
	      }
	    }
	    //cout<<"succ_h:"<<succ_h<<endl;
            succ_node.clear_h_dirty();
            search_progress.inc_evaluated_states();
            search_progress.inc_evaluations(heuristics.size());

            // Note that we cannot use succ_node.get_g() here as the
            // node is not yet open. Furthermore, we cannot open it
            // before having checked that we're not in a dead end. The
            // division of responsibilities is a bit tricky here -- we
            // may want to refactor this later.
            //open_list->evaluate(node.get_g() + get_adjusted_cost(*op), is_preferred);
            open_list->evaluate2(node.get_g() + get_adjusted_cost(*op),succ_h);
	    //cout<<"After evaluate2"<<endl;
            //bool dead_end = open_list->is_dead_end();

            //TODO:CR - add an ID to each state, and then we can use a vector to save per-state information
            if (do_pathmax) {
                if ((node.get_h() - get_adjusted_cost(*op)) > succ_h) {
                    //cout << "Pathmax correction: " << succ_h << " -> " << node.get_h() - get_adjusted_cost(*op) << endl;
                    succ_h = node.get_h() - get_adjusted_cost(*op);
                    heuristics[0]->set_evaluator_value(succ_h);
                    open_list->evaluate(node.get_g() + get_adjusted_cost(*op), is_preferred);
                    search_progress.inc_pathmax_corrections();
                }
            }
            succ_node.open(succ_h, node, op);

            open_list->insert(succ_state.get_id());
	    //cout<<"After insert"<<endl;fflush(NULL);
            if (search_progress.check_h_progress(succ_node.get_g())) {
                reward_progress();
            }
        } else if (succ_node.get_g() > node.get_g() + get_adjusted_cost(*op)) {
	  //cout<<"succ_node is not new"<<endl;fflush(stdout);
            // We found a new cheapest path to an open or closed state.
            if (reopen_closed_nodes) {
                //TODO:CR - test if we should add a reevaluate flag and if it helps
                // if we reopen closed nodes, do that
                if (succ_node.is_closed()) {
		  //cout<<"Reopening closed node"<<endl;fflush(NULL);
                    /* TODO: Verify that the heuristic is inconsistent.
                     * Otherwise, this is a bug. This is a serious
                     * assertion because it can show that a heuristic that
                     * was thought to be consistent isn't. Therefore, it
                     * should be present also in release builds, so don't
                     * use a plain assert. */
                    //TODO:CR - add a consistent flag to heuristics, and add an assert here based on it
                    search_progress.inc_reopened();
                }
		//cout<<"Revalued opened node"<<endl;fflush(NULL);
                succ_node.reopen(node, op);
                //heuristics[0]->set_evaluator_value(succ_node.get_h());
                // TODO: this appears fishy to me. Why is here only heuristic[0]
                // involved? Is this still feasible in the current version?
                //open_list->evaluate(succ_node.get_g(), is_preferred);
                open_list->evaluate2(succ_node.get_g(), succ_node.get_h());

                open_list->insert(succ_state.get_id());
            } else {
                // if we do not reopen closed nodes, we just update the parent pointers
                // Note that this could cause an incompatibility between
                // the g-value and the actual path that is traced back
                succ_node.update_parent(node, op);
            }
        }else{
	  //cout<<"succ_node is not new but g value is higher, so no worries"<<endl;fflush(stdout);
	}
    }
    //cout<<"returning IN_PROGRESS:"<<endl;fflush(stdout);

    return IN_PROGRESS;
}

void EagerSearchS::generateSSCCReport() {
	string dominio = domain_name;
        string tarefa = problem_name2;
        string heuristica = heuristic_name2;
	int threshold = 0;

	//get the deep_F_boundary
        if (deep_F_boundary != 0) {
        	threshold = deep_F_boundary;
        }

	cout<<"dominio = "<<dominio<<endl;
        cout<<"tarefa = "<<tarefa<<endl;
        cout<<"heuristica = "<<heuristica<<endl;
        cout<<"problem_name_gapdb = "<<problem_name_gapdb<<"\n";

	string name, dirDomain, dirSSCC, dirSSCCFile, outputFile;
	if (problem_name_gapdb != "temp") {
                string t0 = tarefa;
                size_t found0 = t0.find(".");
                string fileDir = t0.substr(0, found0);

                string t = problem_name_gapdb;
                size_t found = t.find(".");
                name = t.substr(0, found);
                name += "_F_";

                name+=boost::lexical_cast<std::string>(threshold);
                name += ".csv";

                dirDomain = "mkdir /home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio;
                dirSSCC = "mkdir /home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio+"/bc";
                dirSSCCFile = "mkdir /home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio+"/bc/"+fileDir;
                outputFile = "/home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio+"/bc/"+fileDir+"/"+name;

                if (system(dirDomain.c_str())) {
                        cout<<"Directory: "<<heuristica<<" created."<<endl;
                }

                if (system(dirSSCC.c_str())) {
                        cout<<"Directory: SSCC created."<<endl;
                }

                if (system(dirSSCCFile.c_str())) {
                        cout<<"Directory: SSCCFile created."<<endl;
                }
        } else {
		size_t found = tarefa.find(".");
                cout<<"found = "<<found<<endl;
                name = tarefa.substr(0, found);
                name+="_F_";
                name+=boost::lexical_cast<std::string>(threshold);
                name += ".csv";

                dirDomain = "mkdir /home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio;
                dirSSCC = "mkdir /home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio+"/bc";
                outputFile = "/home/marvin/marvin/astar/"+heuristica+"/reportastar/"+dominio+"/bc/"+name;


                if (system(dirDomain.c_str())) {
                        cout<<"Directory: "<<heuristica<<" created."<<endl;
                }

                if (system(dirSSCC.c_str())) {
                        cout<<"Directory: SSCC created."<<endl;
                }
        }
        cout<<"name = "<<name<<endl;
        ofstream output;
        output.open(outputFile.c_str());
}



pair<SearchNode, bool> EagerSearchS::fetch_next_node() {
    /* TODO: The bulk of this code deals with multi-path dependence,
       which is a bit unfortunate since that is a special case that
       makes the common case look more complicated than it would need
       to be. We could refactor this by implementing multi-path
       dependence as a separate search algorithm that wraps the "usual"
       search algorithm and adds the extra processing in the desired
       places. I think this would lead to much cleaner code. */
  //cout<<"Fetching next_node"<<endl;fflush(stdout);

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
	//cout<<"Fetching next_node finished"<<endl;fflush(stdout);
        return make_pair(node, true);
    }
}

void EagerSearchS::reward_progress() {
    // Boost the "preferred operator" open lists somewhat whenever
    // one of the heuristics finds a state with a new best h value.
    open_list->boost_preferred();
}

void EagerSearchS::dump_search_space() {
    search_space.dump();
}

void EagerSearchS::update_jump_statistic(const SearchNode &node) {
    //cout<<"Started update_jump_statistic"<<endl;fflush(stdout);
    if (f_evaluator) {
      int new_f_value = node.get_g()+node.get_h();
        //heuristics[0]->set_evaluator_value(node.get_h());
        //f_evaluator->evaluate(node.get_g(), false);
        //int new_f_value = f_evaluator->get_value();
        search_progress.report_f_value(new_f_value);
    }
    //cout<<"Finished update_jump_statistic"<<endl;fflush(stdout);
}

void EagerSearchS::print_heuristic_values(const vector<int> &values) const {
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

    EagerSearchS *engine = 0;
    if (!parser.dry_run()) {
        opts.set<bool>("mpd", false);
        engine = new EagerSearchS(opts);
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

    EagerSearchS *engine = 0;
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
        engine = new EagerSearchS(opts);
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

    EagerSearchS *engine = 0;
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
        engine = new EagerSearchS(opts);
    }
    return engine;
}

static Plugin<SearchEngine> _plugin("eager", _parse);
static Plugin<SearchEngine> _plugin_astar("astar_s", _parse_astar);
static Plugin<SearchEngine> _plugin_greedy("eager_greedy", _parse_greedy);
