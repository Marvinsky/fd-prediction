#include "pattern_generation_edelkamp.h"

#include "pdb_heuristic.h"
#include "zero_one_pdbs_heuristic.h"

#include "../blind_search_heuristic.h"
#include "../causal_graph.h"
#include "../globals.h"
#include "../plugin.h"
#include "../rng.h"
#include "../timer.h"
#include "../utilities.h"
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <vector>
#include <ext/hash_set>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace __gnu_cxx;
using namespace std;

//Root and fd information
string _HOME_INFO_PDB = "/home";
string _FD_INFO_PDB = "/fd";

PatternGenerationEdelkamp::PatternGenerationEdelkamp(const Options &opts)
    : task(get_task_from_options(opts)),
      pdb_max_size(opts.get<int>("size")),
      colls(opts.get<int>("colls")),
      num_episodes(opts.get<int>("eps")),
      mutation_probability(opts.get<double>("mp")),
      disjoint_patterns(opts.get<bool>("disjoint")),
      complementary(opts.get<bool>("complementary")),
      cost_type(OperatorCost(opts.get<int>("cost_type"))) {
	static bool no_more_pdb_gen_print=true;
      if(pdb_max_size>=20000000){
	time_limit=600;
      }
      else if(pdb_max_size>=2000000){
	time_limit=40;
      }
      else if(pdb_max_size>=200000){
	time_limit=60; //time_limit=20
      }
      else if(pdb_max_size>=20000){
	time_limit=10;
      }
      else {
	time_limit=5;
      }
      //cout<<"overall GA's pdb_gen_time_limit, currently:"<<pdb_gen_time_limit<<endl;
      if(g_timer()>pdb_gen_time_limit){
	if(no_more_pdb_gen_print==true){
	  cout<<"no more PDB generation, overall time("<<g_timer<<")>pdb_gen_time_limit("<<pdb_gen_time_limit<<")"<<endl;
	  no_more_pdb_gen_print=false;
	}
	//best_heuristic->set_stop_using(true);
	no_more_ga_pdbs=true;
	return;
      }
      if((double(get_peak_memory_in_kb())/1024)>1000){
	cout<<"no more PDB generation, Peak memory above 1 GB max"<<endl;
	no_more_ga_pdbs=true;
	return;
      }
      //to add more variety to patterns
      g_random_seed+=1;
      srand(g_random_seed);
      g_rng.seed(g_random_seed);
      
      ///cout<<"GAPDB,mutation_probability:"<<mutation_probability<<".Disj:"<<disjoint_patterns<<",restarted random sequence for g_random_seed"<<g_random_seed<<",complementary:"<<endl;
      timer.reset();
      //As we are running several methods for comparison, we ran PDB generation only once per problem and then add the generation time(minus the new generation time of just the selected pdbs) and the selected patterns
      cout<<"use_saved_pdbs = "<<use_saved_pdbs<<endl;
      if(use_saved_pdbs){
	//problem_name=g_plan_filename;
      /*for (int i = 0; i < argc_copy; ++i) {
	//puts(argv_copy[i]);
	cout<<"i:"<<i<<","<<argv_copy[i]<<endl;
      }*/
//	vector<int> v1 {9,17,18,31,50,56}; 
//	vector<int> v2 {24,26,43,45,54};
//	vector<int> v3 {2, 16, 19, 25, 46, 51};
//	vector<int> v4 {3, 6, 10, 12, 14, 15, 28, 33, 41, 49};
//	vector<int> v5 {1, 8, 11, 21, 32, 35, 38, 40, 47, 52};
//	vector<int> v6 {7, 23, 27, 34, 39, 44, 55};
//	vector<int> v7 {20, 22, 37, 42, 48, 53};
//	vector<int> v8 {};
	vector<vector<int> > pattern_collection;
//	pattern_collection.push_back(v1);
//	pattern_collection.push_back(v2);
//	pattern_collection.push_back(v3);
//	pattern_collection.push_back(v4);
//	pattern_collection.push_back(v5);
//	pattern_collection.push_back(v6);
//	pattern_collection.push_back(v7);
//	pattern_collection.push_back(v8);

	if(!get_GA_patterns_from_file(pattern_collection,disjoint_patterns,mutation_probability,pdb_max_size)){
	  //cout<<"Cant find at least one previous GA, so returning dummy heuristic from now on"<<endl;
	  no_more_ga_pdbs=true;
	  return;
	}
	//cout<<"returned from get_GA_patterns_from_file"<<endl;fflush(NULL);
	      
	Options opts2;

	opts2.set<TaskProxy *>("task", task);

	opts2.set<int>("cost_type", cost_type);
	opts2.set<bool>("disjoint", disjoint_patterns);
	opts2.set<vector<vector<int> > >("patterns", pattern_collection);
	opts2.set<double>("mp", mutation_probability);
	opts2.set<bool>("complementary", complementary);
	opts2.set<int>("size", pdb_max_size);
	//cout<<"complementary:"<<complementary<<endl;
	ZeroOnePDBsHeuristic *zoppch =
	  new ZeroOnePDBsHeuristic(opts2);
	best_heuristic = zoppch;
	best_fitness = best_heuristic->get_approx_mean_finite_h();
	//dump_best_heuristic();
        } else {
          genetic_algorithm();
	  cout << "Pattern generation (Edelkamp) time: " << timer() << endl;
        }
    //Timer timer;
    //genetic_algorithm();
    //cout << "Pattern generation (Edelkamp) time: " << timer << endl;
}

PatternGenerationEdelkamp::~PatternGenerationEdelkamp() {
}

void PatternGenerationEdelkamp::dump_best_heuristic() const {
	static int i = 0;
	if (best_fitness <0) {
		cout<<"GAPDB dit not generate any complete PDB collection, so stop using!"<<endl;
		best_heuristic->set_stop_using(true);
		return;
	}
	if ((!use_saved_pdbs) && best_fitness_was_duplicate) {
		cout<<"last best_fitness_was_duplicate,";
	}
	cout<<"returning best heuristic(GAPDB)[,"<<i<<",]:";
	best_heuristic->dump();
	cout<<",mp:,"<<mutation_probability;
	cout<<",disjoint_patterns:,"<<disjoint_patterns;
	cout<<",size:,"<<pdb_max_size;
	cout<<"-best_fitness:"<<best_fitness<<",";
	cout<<",initial value:"<<best_heuristic->compute_heuristic(g_initial_state());
	cout<<",GAPDB generation time:"<<timer()<<endl;
	i++;
}

void PatternGenerationEdelkamp::dump_file() const {
	ofstream outputFile;
        cout<<"callind dump_file() with ";
	static int i = 0;
	if (best_fitness<0){
		//No GAPDB was chosen
		return;
	}

	string datstr = "dat_";
	if (ss_probes == 0) {
        	datstr += "info";
	} else if (run_n_heuristics == 0) {
		datstr += boost::lexical_cast<std::string>(ss_probes);
		datstr += "_probes";
  	} else {
		datstr += boost::lexical_cast<std::string>(ss_probes);
		datstr += "_probes_" + boost::lexical_cast<std::string>(run_n_heuristics);
	}
  	cout<<"datstr = "<<datstr<<"\n";
        cout<<"pdb_dump_counter = "<<pdb_dump_counter<<endl;
	if (pdb_dump_counter == 0) {
                //Create directory dat
                string datDirectory = "mkdir "+_HOME_INFO_PDB+"/marvin"+_FD_INFO_PDB+"/"+datstr;
                string domainDirectory = "mkdir "+_HOME_INFO_PDB+"/marvin"+_FD_INFO_PDB+"/"+datstr+"/"+domain_name;
		if (system(datDirectory.c_str())) {
		   cout<<"dat directory created."<<endl;
		} else {
		   cout<<"dat directory already exists."<<endl;
		}

		if (system(domainDirectory.c_str())) {
		   cout<<"domain directory created."<<endl;
		} else {
		   cout<<"domain directory already exists."<<endl;
		}

		//string system_call = "/bin/rm dat/"+domain_name+"/";
		string system_call = "rm "+_HOME_INFO_PDB+"/marvin"+_FD_INFO_PDB+"/"+datstr+"/"+domain_name+"/";
                string task2 = problem_name2;
                size_t found2 = task2.find(".");
                string task2_final = task2.substr(0, found2);
		system_call += task2_final; 
		system_call += ".dat";
		cout<<"First call, removing system_call to avoid duplicate pdbs:"<<system_call<<endl;

		int temp = system(system_call.c_str());
		cout<<"rm status:"<<temp<<endl;
		if (temp == 0) {
			
		}
	}
	pdb_dump_counter++;

        string task3 = problem_name2;
	
	size_t found3 = task3.find(".");
	
	string file_name =  task3.substr(0, found3);
        
	file_name += ".dat";
        file_name = "/" + file_name;
        file_name = domain_name + file_name;
	file_name = _HOME_INFO_PDB+"/marvin"+_FD_INFO_PDB+"/"+datstr+"/" + file_name;
        cout<<"file_name: "<<file_name<<endl;

	outputFile.open(file_name.c_str(), ios::app);
	problem_name = problem_name2;
	outputFile<<problem_name<<":";
	outputFile<<"returning best heuristic(GAPDB)[,"<<i<<",]:";
	i++;
	string patterns;
	best_heuristic->get_patterns(patterns);
	outputFile<<patterns;


	std::ostringstream ss;
	ss<<std::fixed<<std::setprecision(7);
	ss<<mutation_probability;
	outputFile<<",mp:,"<<ss.str();
	outputFile<<",size:,"<<pdb_max_size;
	outputFile<<",disjoint_patterns:,"<<disjoint_patterns;
	outputFile<<"-best_fitness:"<<best_fitness<<",";
	outputFile<<",initial value:"<<best_heuristic->compute_heuristic(g_initial_state());
	outputFile<<",GAPDB generation time:"<<timer()<<endl;
	outputFile.flush();
	outputFile.close();
} 

void PatternGenerationEdelkamp::select(const vector<double> &fitness_values) {
    vector<double> cumulative_fitness;
    cumulative_fitness.reserve(fitness_values.size());
    double total_so_far = 0;
    for (size_t i = 0; i < fitness_values.size(); ++i) {
        total_so_far += fitness_values[i];
        cumulative_fitness.push_back(total_so_far);
    }
    // total_so_far is now sum over all fitness values

    vector<vector<vector<bool> > > new_pattern_collections;
    new_pattern_collections.reserve(colls);
    for (int i = 0; i <colls; ++i) {
        int selected;
        if (total_so_far == 0) {
            // All fitness values are 0 => choose uniformly.
            selected = g_rng(fitness_values.size());
        } else {
            double random = g_rng() * total_so_far; // [0..total_so_far)
            // Find first entry which is strictly greater than random.
            selected = upper_bound(cumulative_fitness.begin(),
                                   cumulative_fitness.end(), random) -
                       cumulative_fitness.begin();
        }
        new_pattern_collections.push_back(pattern_collections[selected]);
    }
    pattern_collections.swap(new_pattern_collections);
}

void PatternGenerationEdelkamp::mutate() {
    for (size_t i = 0; i < pattern_collections.size(); ++i) {
        for (size_t j = 0; j < pattern_collections[i].size(); ++j) {
            vector<bool> &pattern = pattern_collections[i][j];
            for (size_t k = 0; k < pattern.size(); ++k) {
                double random = g_rng(); // [0..1)
                if (random < mutation_probability) {
                    pattern[k].flip();
                }
            }
        }
    }
}

void PatternGenerationEdelkamp::transform_to_pattern_normal_form(const vector<bool> &bitvector,
                                                                 vector<int> &pattern) const {
    for (size_t i = 0; i < bitvector.size(); ++i) {
        if (bitvector[i])
            pattern.push_back(i);
    }
}

void PatternGenerationEdelkamp::remove_irrelevant_variables(
    vector<int> &pattern) const {
    hash_set<int> in_original_pattern(pattern.begin(), pattern.end());
    hash_set<int> in_pruned_pattern;

    vector<int> vars_to_check;
    for (size_t i = 0; i < g_goal.size(); ++i) {
        int var_no = g_goal[i].first;
        if (in_original_pattern.count(var_no)) {
            // Goals are causally relevant.
            vars_to_check.push_back(var_no);
            in_pruned_pattern.insert(var_no);
        }
    }

    while (!vars_to_check.empty()) {
        int var = vars_to_check.back();
        vars_to_check.pop_back();
        // A variable is relevant to the pattern if it is a goal variable or if
        // there is a pre->eff arc from the variable to a relevant variable.
        // Note that there is no point in considering eff->eff arcs here.
        const vector<int> &rel = g_causal_graph->get_eff_to_pre(var);
        for (size_t i = 0; i < rel.size(); ++i) {
            int var_no = rel[i];
            if (in_original_pattern.count(var_no) &&
                !in_pruned_pattern.count(var_no)) {
                // Parents of relevant variables are causally relevant.
                vars_to_check.push_back(var_no);
                in_pruned_pattern.insert(var_no);
            }
        }
    }

    pattern.assign(in_pruned_pattern.begin(), in_pruned_pattern.end());
    sort(pattern.begin(), pattern.end());
}

bool PatternGenerationEdelkamp::is_pattern_too_large(
    const vector<int> &pattern) const {
    // test if the pattern respects the memory limit
    int mem = 1;
    for (size_t i = 0; i < pattern.size(); ++i) {
        int domain_size = g_variable_domain[pattern[i]];
        if (!is_product_within_limit(mem, domain_size, pdb_max_size))
            return true;
        mem *= domain_size;
    }
    return false;
}

bool PatternGenerationEdelkamp::mark_used_variables(
    const vector<int> &pattern, vector<bool> &variables_used) const {
    for (size_t i = 0; i < pattern.size(); ++i) {
        int var_no = pattern[i];
        if (variables_used[var_no])
            return true;
        variables_used[var_no] = true;
    }
    return false;
}

void PatternGenerationEdelkamp::evaluate(vector<double> &fitness_values) {
    for (size_t i = 0; i < pattern_collections.size(); ++i) {
        //cout << "evaluate pattern collection " << (i + 1) << " of " << pattern_collections.size() << endl;
        double fitness = 0;
        bool pattern_valid = true;
        vector<bool> variables_used(g_variable_domain.size(), false);
        vector<vector<int> > pattern_collection;
        pattern_collection.reserve(pattern_collections[i].size());
        for (size_t j = 0; j < pattern_collections[i].size(); ++j) {
            const vector<bool> &bitvector = pattern_collections[i][j];
            vector<int> pattern;
            transform_to_pattern_normal_form(bitvector, pattern);

            if (is_pattern_too_large(pattern)) {
                //cout << "pattern " << j << " exceeds the memory limit!" << endl;
                pattern_valid = false;
                break;
            }

            if (disjoint_patterns) {
                if (mark_used_variables(pattern, variables_used)) {
                    //cout << "patterns are not disjoint anymore!" << endl;
                    pattern_valid = false;
                    break;
                }
            }

            remove_irrelevant_variables(pattern);
            pattern_collection.push_back(pattern);
        }
        if (!pattern_valid) {
            // set fitness to a very small value to cover cases in which all patterns are invalid
            fitness = 0.001;
        } else {
            // generate the pattern collection heuristic and get its fitness value.
	    if(timer()>time_limit){
	      //cout<<"breaking-1 out of GA Algortihm, current gen_time:"<<timer<<" bigger than time_limit:"<<time_limit<<endl;
	      timer.stop();
	      break;
	    }
            Options opts;
            opts.set<TaskProxy *>("task_proxy", task);
            opts.set<int>("cost_type", cost_type);
            opts.set<bool>("disjoint", disjoint_patterns);
            opts.set<vector<vector<int> > >("patterns", pattern_collection);
	    opts.set<double>("mp", mutation_probability);
	    opts.set<bool>("complementary", complementary);
	    opts.set<int>("size", pdb_max_size);
            ZeroOnePDBsHeuristic *zoppch =
                new ZeroOnePDBsHeuristic(opts);
            fitness = zoppch->get_approx_mean_finite_h();
            // update the best heuristic found so far.
            if (fitness > best_fitness) {
	      std::pair<std::set<vector<vector<int> > >::iterator,bool> ret;

	      ret=chosen_pattern_collections.insert(pattern_collection); // all current pattern collections
	      if(ret.second==false){
		//cout<<"pattern collection already selected"<<endl;
		//zoppch->dump();
                delete zoppch;
		best_fitness_was_duplicate=true;
	      }
	      else{
                best_fitness = fitness;
                cout << "best_fitness = " << best_fitness << endl;
                delete best_heuristic;
                best_heuristic = zoppch;
		best_fitness_was_duplicate=false;
	      }
            } else {
                delete zoppch;
            }
        }
        fitness_values.push_back(fitness);
    }
}

void PatternGenerationEdelkamp::bin_packing() {
    vector<int> variables;
    variables.reserve(g_variable_domain.size());
    for (size_t i = 0; i < g_variable_domain.size(); ++i) {
        variables.push_back(i);
    }

    for (int i = 0; i < colls; ++i) {
        // random variable ordering for all pattern collections
        g_rng.shuffle(variables);
        vector<vector<bool> > pattern_collection;
        vector<bool> pattern(g_variable_name.size(), false);
        int current_size = 1;
        for (size_t j = 0; j < variables.size(); ++j) {
            int var = variables[j];
            int next_var_size = g_variable_domain[var];
            if (next_var_size > pdb_max_size) // var never fits into a bin
                continue;
            if (!is_product_within_limit(current_size, next_var_size, pdb_max_size)) {
                // Open a new bin for var.
                pattern_collection.push_back(pattern);
                pattern.clear();
                pattern.resize(g_variable_name.size(), false);
                current_size = 1;
            }
            current_size *= next_var_size;
            pattern[var] = true;
        }
        // the last bin has not bin inserted into pattern_collection, do so now.
        // We test current_size against 1 because this is cheaper than
        // testing if pattern is an all-zero bitvector. current_size
        // can only be 1 if *all* variables have a domain larger than
        // pdb_max_size.
        if (current_size > 1) {
            pattern_collection.push_back(pattern);
        }
        pattern_collections.push_back(pattern_collection);
    }
}

void PatternGenerationEdelkamp::genetic_algorithm() {
    best_fitness = -1;
    best_heuristic = 0;
    bin_packing();
    //cout << "initial pattern collections:" << endl;
    //dump();
    vector<double> initial_fitness_values;
    evaluate(initial_fitness_values);
    for (int i = 0; i < num_episodes; ++i) {
	if(timer()>time_limit){
	  //cout<<"breaking-3 out of GA Algortihm, current gen time:"<<timer()<<" bigger than time_limit:"<<time_limit<<endl;
	  timer.stop();
	  break;
        }
        //cout << endl;
        //cout << "--------- episode no " << (i + 1) << " ---------" << endl;
        mutate();
        //cout << "current pattern_collections after mutation" << endl;
        //dump();
        vector<double> fitness_values;
        evaluate(fitness_values);
	if(timer()>time_limit){
	  //cout<<"breaking-4 out of GA Algortihm, current gen time:"<<timer()<<" bigger than time_limit:"<<time_limit<<endl;
	  timer.stop();
	  break;
	}
        select(fitness_values); // we allow to select invalid pattern collections
        //cout << "current pattern collections (after selection):" << endl;
        //dump();
    timer.stop();//no need to keep this one ticking after pattern generation finished!
    }
}

void PatternGenerationEdelkamp::dump() const {
    for (size_t i = 0; i < pattern_collections.size(); ++i) {
        cout << "pattern collection no " << (i + 1) << endl;
        for (size_t j = 0; j < pattern_collections[i].size(); ++j) {
            cout << pattern_collections[i][j] << endl;
        }
    }
}

static Heuristic *_parse(OptionParser &parser) {
    parser.document_synopsis(
        "Genetic Algorithm PDB",
        "The following paper describes the automated creation of pattern databases "
        "with a genetic algorithm. Pattern collections are initially created with a "
        "bin-packing algorithm. The genetic algorithm is used to optimize the pattern "
        "collections with an objective function that estimates the mean heuristic "
        "value of the the pattern collections. Pattern collections with higher mean "
        "heuristic estimates are more likely selected for the next generation.\n\n"
        " * Stefan Edelkamp<<BR>>"
        " [Automated Creation of Pattern Database Search Heuristics "
        "http://www.springerlink.com/content/20613345434608x1/].<<BR>>"
        "In //Proceedings of the 4th Workshop on Model Checking and Artificial Intelligence "
        "(!MoChArt 2006)//, pp. 35-50, 2007.");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "not supported");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes");
    parser.document_property("consistent", "yes");
    parser.document_property("safe", "yes");
    parser.document_property("preferred operators", "no");
    parser.document_note(
        "Note",
        "This pattern generation method uses the zero/one pattern database heuristic.");
    parser.document_note(
        "Implementation Notes",
        "The standard genetic algorithm procedure as described in the paper is "
        "implemented in Fast Downward. The implementation is close to the paper.\n\n"
        "+ Initialization<<BR>>"
        "In Fast Downward bin-packing with the next-fit strategy is used. A bin "
        "corresponds to a pattern which contains variables up to ``pdb_max_size``. "
        "With this method each variable occurs exactly in one pattern of a collection. "
        "There are ``colls`` collections created.\n"
        "+ Mutation<<BR>>"
        "With probability ``mutation_probability`` a bit is flipped meaning that "
        "either a variable is added to a pattern or deleted from a pattern.\n"
        "+ Recombination<<BR>>"
        "Recombination isn't implemented in Fast Downward. In the paper recombination "
        "is described but not used.\n"
        "+ Evaluation<<BR>>"
        "For each pattern collection the mean heuristic value is computed. For a "
        "single pattern database the mean heuristic value is the sum of all pattern "
        "database entries divided through the number of entries. Entries with infinite "
        "heuristic values are ignored in this calculation. The sum of these individual "
        "mean heuristic values yield the mean heuristic value of the collection.\n"
        "+ Selection<<BR>>"
        "The higher the mean heuristic value of a pattern collection is, the more "
        "likely this pattern collection should be selected for the next generation. "
        "Therefore the mean heuristic values are normalized and converted into "
        "probabilities and Roulette Wheel Selection is used.\n"
        "+\n\n", true);

    parser.add_option<int>("size", "maximal number of states per pattern database ", "50000");
    parser.add_option<int>("colls", "number of pattern collections to maintain in the genetic algorithm (population size)", "5");
    parser.add_option<int>("eps", "number of episodes for the genetic algorithm", "30");
    parser.add_option<double>("mp", "probability between 0 and 1 for flipping a bit in the genetic algorithm", "0.01");
    parser.add_option<bool>("disjoint","using disjoint variables in the patterns of a collection", "false");
    parser.add_option<bool>("disjoint", "consider a pattern collection invalid (giving it very low fitness) if its patterns are not disjoint", "false");
    parser.add_option<bool>("complementary", "If restarting sampling, was the heuristic already defined as complementary", "false");


    Heuristic::add_options_to_parser(parser);
    Options opts = parser.parse();
    if (parser.help_mode())
        return 0;

    if (opts.get<int>("size") < 1)
        parser.error("size per pdb must be at least 1");
    if (opts.get<int>("colls") < 1)
        parser.error("number of pattern collections must be at least 1");
    if (opts.get<int>("eps") < 0)
        parser.error("number of episodes must be a non negative number");
    if (opts.get<double>("mp") < 0 || opts.get<double>("mp") > 1)
        parser.error("mutation probability must be in [0..1]");

    if (parser.dry_run())
        return 0;
    PatternGenerationEdelkamp pge(opts);
  if(no_more_ga_pdbs||pge.get_fitness()<0){
    //cout<<"We must have reached maximum time or memory limit, returning blind heuristic instead"<<endl;
    if(pge.get_fitness()<0){
      cout<<"No unique pattenrs were found!"<<endl;
    }
    Heuristic::add_options_to_parser(parser);
    Options opts = parser.parse();
    BlindSearchHeuristic *dummy_heur =
                new BlindSearchHeuristic(opts);
    dummy_heur->set_stop_using(true);
    return dummy_heur; 
  }
    //we copy the patterns to be able to reuse them if we crast out before time is up
    if(!use_saved_pdbs){//only if we are not reading the pdbs!
      pge.dump_file();
      pge.dump_best_heuristic();
    }
    return pge.get_pattern_collection_heuristic();
}

static Plugin<Heuristic> _plugin("gapdb", _parse);
