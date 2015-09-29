#include "ss_search.h"

#include "../globals.h"
#include "../heuristic.h"
#include "../option_parser.h"
#include "../successor_generator.h"
#include "../g_evaluator.h"
#include "../sum_evaluator.h"
#include "../plugin.h"
#include "../rng.h"

#include <iostream>
#include <fstream>
#include "../ext/boost/dynamic_bitset.hpp"
#include <boost/lexical_cast.hpp>
#include "../global_state.h"
#include <stdlib.h>

bool ss_debug=false;
int next_f_bound=INT_MAX/2;
bool domination_check=false;
set<vector<int> > F_culprits;
double sampling_time_limit=150;
double overall_time_limit=1400;
int N_default = 8;
//run experiments
string run_method = "both";

//Root and fd information
string _HOME_INFO = "/home";
string _FD_INFO = "/fd";
bool run_min_heuristic = false;//true=select from all the heuristics/false=select just gapdb
bool run_min_eval_time_approach = false;//true run min_eval_time in order to use time to apply meta-reasoning

//#define _SS_DEBUG
//#define _LMCUT_EARLY_TERM

SSSearch::SSSearch(const Options &opts) : SearchEngine(opts), current_state(g_initial_state()) {
//first determine if f_boundary was manually set
if(f_boundary!=0){
  cout<<"f_boundary manually set to "<<f_boundary<<endl;
}  
	ScalarEvaluator * evaluator = opts.get<ScalarEvaluator *>("eval");
	std::set<Heuristic *> hset;
	evaluator->get_involved_heuristics(hset);
	int heuristics_active=0;
	for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); it++) {
#ifdef _LMCUT_EARLY_TERM
	  size_t found_lmcut = (*it)->get_heur_name().find("lmcut");
#endif
	  //Eliminate any heuristics which were not generated because we ran out of time
	  //currently this is hacked to return not using heuristics
	  if((*it)->is_using()){
#ifdef _LMCUT_EARLY_TERM
	    if(found_lmcut!=std::string::npos){
	      cout<<"lmcut is present!, doing early term sampling"<<endl;fflush(stdout);
	      lmcut_heuristic.push_back(*it);
	    }
	    else{
#endif
	      heuristics.push_back(*it);
#ifdef _LMCUT_EARLY_TERM
	    }
#endif
            all_heuristics.push_back(*it);
	    heuristics_active++;
	  }
	}
	assert(heuristics.size() == 1);

	int min_h = INT_MAX/2;
	int max_h = 0;
        for (size_t i = 0; i < heuristics.size(); i++) {
            heuristics[i]->evaluate(g_initial_state());
            int aux = heuristics[i]->get_heuristic();
            if (min_h > aux) {
		min_h = aux;
                heuristic = heuristics[i];
	    }
	    max_h=max(aux,max_h);
	}

	//Beginning timing node generation time
	Timer heur_timings;
	double node_counter=0;
        std::vector<const GlobalOperator*> applicable_ops;
	const GlobalOperator *op;
	double HUST_TPN=0;
        //const GlobalState &s = g_initial_state();
        StateID initial_state_id = g_initial_state().get_id();
	GlobalState s= g_state_registry->lookup_state(initial_state_id);
	std::cout << "scientific:\n" << std::scientific;
	while(heur_timings()<1.0){//||node_counter<1000)//in situ measure node generation costs
	  node_counter++;
	  //cout<<"node_counter:"<<node_counter<<endl;
	  applicable_ops.clear();
	  g_successor_generator->generate_applicable_ops(s, applicable_ops);
	  if(applicable_ops.size()==0){//dead end so start again with root node
	    s= g_state_registry->lookup_state(initial_state_id);
	    //cout<<"restarting random node generation, just hit a dead end"<<endl;
	    continue;
	  }
	  op = applicable_ops[rand()%applicable_ops.size()];//choose operator at random
	  GlobalState succ_state =  g_state_registry->get_successor_state(s, *op);
	  s= g_state_registry->lookup_state(succ_state.get_id());
	}
	node_gen_and_exp_cost=heur_timings.stop()/node_counter;
	cout<<"node_gen_and_exp_cost:"<<node_gen_and_exp_cost<<endl;
	node_time_adjusted_reval=0.5/node_gen_and_exp_cost;
	cout<<"node_time_adjusted_reval="<<node_time_adjusted_reval<<"\n";
        cout<<"starting timing individual heuristics"<<endl;fflush(stdout);
        //Time each heuristic
	double max_TPN=0;
	//vector<pair<int,int> > lp_pdbs;//we greedily remove any lp_pdb with a higher pattern_size if the initial f-value is not increased as a result
	//int lp_pdbs_max_h=0;
	//int max_lp_pattern_size=0;
	vector<int> initial_h_values;
	s= g_state_registry->lookup_state(initial_state_id);
	std::cout << "scientific:\n" << std::scientific;
	double time_limit=0;
	if(heuristics.size()>10000){
	  time_limit=0.0001;
	}
	else if(heuristics.size()>1000){
	  time_limit=0.001;
	}
	else if(heuristics.size()>100){
	  time_limit=0.01;
	}
	else{
	  time_limit=0.1;
	}
	
	for (size_t i = 0; i < heuristics.size(); i++){
	   // heuristics[i]->evaluate(*g_initial_state);
	    //cout<<"h[,"<<i<<",] time_cost is,"<<heuristics[i]->get_time_cost()<<endl;
	      heur_timings.reset();heur_timings.resume();
	      double counter=0;
	      double counter_limit=1000;
	      while(heur_timings()<time_limit){
		while(counter<counter_limit){//in situ measure heuristic evaluation costs
		  counter++;
		  heuristics[i]->evaluate(s);
		}
		counter_limit+=1000;
	      }
	    
	      double TPN=double(heur_timings())/counter;
	      heur_timings.stop();
	      cout<<"heur_timings:"<<heur_timings()<<",counter:"<<counter<<",TPN:"<<TPN<<endl;
	      max_TPN=max(TPN,max_TPN);
	      aggr_TPN+=TPN;
	      HUST_TPN+=TPN;
	        
	      heuristics[i]->set_measured_TPN(TPN);
	      cout<<"h[,"<<i<<",] is:,"<<heuristics[i]->get_heur_name()<<",measured time cost:"<<heuristics[i]->get_measured_TPN()<<",h:"<<heuristics[i]->get_heuristic()<<endl;
	}
	//Ending timing node generation time
	      
	cout<<"max_h before lmcut:"<<max_h<<endl;fflush(stdout);
#ifdef _LMCUT_EARLY_TERM
	cout<<"_LMCUT_EARLY_TERM:"<<endl;fflush(stdout);
	if(lmcut_heuristic.size()>0){
	  cout<<"Before evaluate"<<endl;fflush(stdout);
	  lmcut_heuristic[0]->evaluate(g_initial_state());
	  cout<<"After evaluate"<<endl;fflush(stdout);
	  if (min_h > lmcut_heuristic[0]->get_heuristic()){
	    min_h = lmcut_heuristic[0]->get_heuristic();
	  }
	  max_h=max(lmcut_heuristic[0]->get_heuristic(),max_h);
	}
#endif
        cout<<"initial min_h(constructor) = "<<min_h<<",max_h:"<<max_h<<",heuristics_active:"<<heuristics_active<<",overall input heuristics:"<<hset.size()<<endl;
	sampler = new TypeSystem(heuristic);
	//this->RanGen2 = new CRandomMersenne((unsigned)time(NULL));        
	this->RanGen2 = new CRandomMersenne(1);        
	cout<<"random seed set to 1"<<endl;
	if(!domination_check){
	  if(heuristics.size()>1000){
	    ss_probes=20;
	  }
	  else if(heuristics.size()>500){
	    ss_probes=50;
	  }
	  else if(heuristics.size()>250){
	    ss_probes=100;
	  }
	  else if(heuristics.size()>100){
	    ss_probes=250;
	  }
	  else{
	    ss_probes=500;
	  }
	  cout<<"Not doing domination_check, setting probes to :"<<ss_probes<<endl;
	}
}

SSSearch::~SSSearch() {
}


SearchStatus SSSearch::step() {
  updateGlobalVariables();//update global variables
  bool at_least_one_dominated=false;
  vector<int> demotted_heurs;
  int original_threshold=f_boundary;
  while(search_time()<sampling_time_limit){
	this->RanGen2 = new CRandomMersenne(1);
        	
	//clear all data from previous iteration
	queue.clear();expanded.clear();generated.clear();
	if(!domination_check){
	  //full_collector=collector;
	  //collector.clear();
	  //cout<<"Cleared F_culprits&b_culprits"<<endl;F_culprits.clear();
	}
        predict(ss_probes);
	if(next_f_bound==INT_MAX/2){
	  cout<<"next_f_bound was not updated!, check code!"<<endl;
	  cout<<"For now just selecting best_heuristics_greedy"<<endl;
	  runReports(true); 
	  exit(1);
	}

	threshold=2*next_f_bound;
	f_boundary=threshold;//want to sample the biggest areas of search space
	next_f_bound=INT_MAX/2;
	cout<<"Restarting sampling for Thresshold:"<<threshold<<",overall search_time:"<<search_time()<<endl;
	if(search_time()>5.0&&domination_check){
	  domination_check=false;
	  bool stop_checking_for_dominatedness=false;
	  cout<<"One time Checking dominatedness of heuristics, one time check at 10 seconds sampling time to get rid of any heuristic who is dominated by the rest, only do this when all heurs have similar evaluation time!"<<endl;fflush(stdout);
	  cout<<"domination check, F_culprits:"<<F_culprits.size()<<",b_culprits:"<<collector.size()<<endl;
	  cout<<"domination checks start, time:"<<search_time()<<endl;double start_time_dom_check=search_time();
	  set<int> undominated_heuristics;
	  for (size_t i = 0; i < heuristics.size(); i++){
	    undominated_heuristics.insert(i);
	  }
	  set<int> remaining_heurs_to_check=undominated_heuristics;
	  set<pair<int,int> >skip_checking;
	  //vector<int> dominator_counter(undominated_heuristics.size());

	  //Quick check using b_culprits to eliminate obvious undominated checks
	  cout<<"search time before doing quick b_culprits undomination check"<<search_time()<<endl;

	  for (size_t i = 0; i < remaining_heurs_to_check.size(); i++){
	    cout<<"Checking b_culprits for h:"<<i<<endl;
	      for (size_t j = i+1; j < remaining_heurs_to_check.size(); j++){
		
		bool found_pair=false;
		for (std::map<boost::dynamic_bitset<>,double>::iterator it=collector.begin(); it!=collector.end(); it++){
		  if((!it->first.test(i))&&it->first.test(j)){//we know h(i) can not be dominated by h(j) because i is pruning and j is expanding
		      skip_checking.insert(make_pair(i,j));
		      //cout<<"h("<<i<<") is prunning but h("<<j<<") is expanding, so undominated pair"<<endl;
		      found_pair=true;
		      //dominator_counter[i]++;
		  }//and viceversa, should halve the loop size
		  else if((!it->first.test(j))&&it->first.test(i)){//we know h(i) can not be dominated by h(j) because i is pruning and j is expanding
		      skip_checking.insert(make_pair(j,i));;
		      //cout<<"inverse h("<<j<<") is prunning but h("<<i<<") is expanding, so undominated pair"<<endl;
		      //dominator_counter[j]++;
		      break;
		  }
		  if(found_pair){
		    break;
		  }
	      }
	    }
	  }
	  cout<<"search time after doing quick b_culprits undomination check"<<search_time()<<",found "<<skip_checking.size()<<" undominating pairs"<<endl;

	  cout<<"F_boundary"<<threshold<<",checking for dominated heuri) stics,search_time:"<<search_time();
	  cout<<"F_culprits:"<<F_culprits.size()<<",b_culprits:"<<collector.size()<<endl;

	  bool just_demoted=false;
	  while(true){
	    at_least_one_dominated=false;
	    int last_heuristic=-1;
	  for (std::set<int>::iterator it2=undominated_heuristics.begin(); it2!=undominated_heuristics.end();){
	    just_demoted=false;
	    ///starting with true condition, in case 2 identical heurs get deleted as in openstacks problems 2-15
	    //in that case all heuristics give the same f-values, they are useless but we still have to keep one
	    //if all are identical
	    for (std::set<int>::iterator it3=remaining_heurs_to_check.begin(); it3!=remaining_heurs_to_check.end();it3++){
	      stop_checking_for_dominatedness=false;
	      just_demoted=true;
	      if(undominated_heuristics.find(*it3)==undominated_heuristics.end()){//cant compare to already demoted heurs!
		//cout<<"h("<<*it3<<" is already demoted"<<endl;
		just_demoted=false;
		continue;
	      }
	      /*  if(dominator_counter[*it2]==heuristics.size()-1){
		cout<<"skipping h("<<*it2<<") it is undominated from b_counters"<<endl;
	      }
	      else{
		//cout<<"h is undominated,"<<dominator_counter[*it2]<<",heuristics:"<<heuristics.size()<<endl;
	      }*/
	      if(skip_checking.find(make_pair(*it2,*it3))!=skip_checking.end()){
		  //cout<<"skipping pair("<<*it2<<","<<*it3<<")"<<endl;
		  just_demoted=false;
		  continue;
	      }
	      if(*it2==*it3){
		just_demoted=false;
		continue;//Do not check heuristic against self!!!
	      }
	      last_heuristic=*it3;
	    //cout<<"Working on heuristics "<<*it2<<","<<*it3<<endl;
	      //cout<<"F_culprits.size:"<<F_culprits.size()<<endl;
    
	      for (std::set<vector<int> >::iterator it=F_culprits.begin(); it!=F_culprits.end(); it++){
		//check that h1 (*it2) is pruning and it3 is expanding, then h2 is undominated
		//cout<<"it:"<<*it<<",it2:"<<*it2;fflush(stdout);cout<<"it3:"<<*it3<<endl;fflush(stdout);
		if((*it)[*it2]>(*it)[*it3]){
		  stop_checking_for_dominatedness=true;
		  just_demoted=false;
		  //cout<<"h("<<*it2<<") f:"<<(*it)[*it2]<<">h("<<*it3<<"):"<<(*it)[*it3]<<",so it is undominated by h("<<*it3<<endl;
		  break;
		}
		  //cout<<"h("<<*it2<<") f:"<<(*it)[*it2]<<">h("<<*it3<<"):"<<(*it)[*it3]<<",so it is still undominated"<<endl;

	      }
	      if(stop_checking_for_dominatedness){
		break;
	      }
	    }
	    if(just_demoted){
	      //skip_checking.insert(make_pair(last_heuristic,*it2));//so reverse not possible, in case heuristics have same f values all over!!
	      heuristics[*it2]->set_stop_using(true);
	      demotted_heurs.push_back(*it2);
	      cout<<"h("<<*it2<<" was dominated by heuristic "<<last_heuristic<<", demoting it to weak,time;"<<search_time()<<endl;
	      at_least_one_dominated=true;
	      undominated_heuristics.erase(it2++);//removing heur from checking, this heur is already demoted
	      //also we now now that heur is not undominated to any of our strong heurisitcs
	    }
	    else{
		//cout<<"h("<<*it2<<") is not dominated by any of the other heuristics"<<endl;
	      it2++;
	    }
	    /*std::set<int>::iterator it_to_erase;
	    cout<<"also removing all dominated heuristics from remaining heurs,just_demoted size:"<<just_demotted.size()<<endl;
	    for (unsigned j=0;j<just_demotted.size();j++){
	      it_to_erase=remaining_heurs_to_check.find(just_demotted[j]);
	      if(it_to_erase!=remaining_heurs_to_check.end()){
		remaining_heurs_to_check.erase(it_to_erase);
		cout<<"removed "<<just_demotted[j]<<"from remaining_heurs_to_check"<<endl;
	      }
	      just_demotted.clear();
	    }*/
	  }
	  if(at_least_one_dominated){
	    cout<<"Running test again because at least one heuristic was dominated"<<endl;
	  }
	  else{
	    cout<<"Finished running domination check, no more heuristics are dominated"<<endl;
	    break;
	  }
	  if(search_time()-start_time_dom_check>60.0){
	    cout<<"Domination check limited to 60 secs, need to improve to enable full check!"<<endl;
	    break;
	  }
	  
	}
	  cout<<"Cleared F_culprits"<<endl;F_culprits.clear();
	  vector<Heuristic *> temp_heuristics;
	  for(size_t i=0;i<heuristics.size();i++){
	    if(heuristics[i]->is_using()){
	      temp_heuristics.push_back(heuristics[i]);
	    }
	  }
	  cout<<"Finished checking for dominated_heuristics, heuristics left:"<<temp_heuristics.size()<<",demotted_heurs:"<<demotted_heurs.size()<<",search_time:"<<search_time()<<endl;
	  //update list of heuristics to the remaining ones
	  cout<<"NEED TO FREE MEMORY OF UNUSED HEURISTICS!!!"<<endl;
	  heuristics=temp_heuristics;
	  cout<<"domination checks finished, time:"<<search_time()<<endl;
	  cout<<"domination checks overall time:"<<search_time()-start_time_dom_check<<endl;
	  cout<<"setting threshold back to "<<original_threshold<<endl;
	  threshold=original_threshold;f_boundary=threshold;
	  if(heuristics.size()>1000){
	    ss_probes=20;
	  }
	  else if(heuristics.size()>500){
	    ss_probes=50;
	  }
	  else if(heuristics.size()>250){
	    ss_probes=100;
	  }
	  else if(heuristics.size()>100){
	    ss_probes=250;
	  }
	  else{
	    ss_probes=500;
	  }
	  cout<<"setting probes to "<<ss_probes<<endl;
      }
      if(!domination_check){
	boost::dynamic_bitset<> b_zeros(heuristics.size());b_zeros.reset();
	if(collector.find(b_zeros)!=collector.end()){
	    if(collector[b_zeros]>1000000000){
	      cout<<"minimum expanded nodes==1,000,000,000,stopping sampling"<<endl;
	      break;
	    }
	    else{
	      cout<<"collector["<<b_zeros<<"]:"<<collector[b_zeros]<<endl;
	    }
	}
      }
  }
  cout<<"selecting best heuristic after,"<<search_time()<<", seconds"<<endl;
  runReports(true); 
        return SOLVED;
}

void SSSearch::predict(int probes) {
        totalPrediction = 0;
	totalPredictionMean = 0;
        cout<<"#probes : "<<probes<<",g_timer:"<<g_timer()<<endl;
	cout<<"input heuristics:"<<heuristics.size()<<endl;
	last_probe=0;
	last_n_expanded=0;
        for (int i = 0; i < probes; i++) {
	  cout<<"#probe:"<<i<<",g_timer:"<<g_timer()<<",search_time:"<<search_time()<<endl;
	  if(search_time()>sampling_time_limit||g_timer()>overall_time_limit){
	    cout<<"Search_timer past maximum sampling_time"<<endl;
	    cout<<"selecting best heuristic after search_time: "<<search_time()<<", seconds,g_timer:"<<g_timer()<<endl;
	    runReports(true); 
	    exit(0);
	  }
	  else if(search_time()>5.0&&domination_check){
	    cout<<"Search_time:"<<search_time()<<",returning for domination_check"<<endl;
	    return;
	  }
            vweight.clear();
	    vmeanheur.clear();
	    //Validate that the number of probes do not exceed the order of 150
	    if (last_n_expanded > 1*pow(10,150)) {
		runReports(true);	
		exit(0);
	    } 

            probe();
            double p = getProbingResult();
	    double mean = getMeanHeurResult();
            totalPrediction = totalPrediction + (p - totalPrediction)/(i + 1);
	    totalPredictionMean = totalPredictionMean + (mean - totalPredictionMean)/(i + 1);
            cout<<"**********"<<endl;
            cout<<"p = "<<p<<endl;
            cout<<"prePre_"<<(i+1)<<" = "<<totalPrediction<<endl;
	    cout<<"mean = "<<mean<<endl;
            cout<<"preMean_"<<(i+1)<<" = "<<totalPredictionMean<<endl;
            cout<<"**********"<<endl;
	    last_probe=i;
	    last_n_expanded=p;
        }
	cout<<"\tnext_f_bound:"<<next_f_bound<<"\n";
	ss_timer_value = ss_timer();
        cout<<"\ttotalPrediction : "<<totalPrediction<<"\n";
	cout<<"\tss_timer: "<<ss_timer_value<<"\n";
	cout<<"\tprobes: "<<probes<<"\n";

	runReports(false);
	//generateGeneratedReport();
        //generateExpandedReport(false);
}

void SSSearch::probe()
{
  static int call_number=0;
  call_number++;
	/*
	 * Probing is done based on the types of the children of the root state
	 */

        queue.clear();
	// evaluating the initial state
	boost::dynamic_bitset<> b_initial_v(heuristics.size()+lmcut_heuristic.size()); 
        int min_h = INT_MAX/2;
	int max_h=0;
	vector<int> h_initial_v;
        for (size_t i = 0; i < heuristics.size(); i++) {
	  heuristics[i]->evaluate(g_initial_state());
	  if(heuristics[i]->is_dead_end()){
	    cout<<"initial state is dead_end, no solution possible!"<<endl;exit(0);
	      min_h=min(min_h,INT_MAX/2);
	      h_initial_v.push_back(INT_MAX/2);
	      }
	  else{
	    min_h = min(min_h, heuristics[i]->get_heuristic());
	    max_h = max(max_h, heuristics[i]->get_heuristic());
	    h_initial_v.push_back(heuristics[i]->get_heuristic());            	
	  }
        }

#ifdef _LMCUT_EARLY_TERM
        if(lmcut_heuristic.size()>0){
	  lmcut_heuristic[0]->evaluate(g_initial_state());
	  if(lmcut_heuristic[0]->is_dead_end()){
	    cout<<"initial state is dead_end(lmcut), no solution possible!"<<endl;exit(0);
	      min_h=min(min_h,INT_MAX/2);
	      h_initial_v.push_back(INT_MAX/2);
	      }
	  else{
	    min_h = min(min_h, lmcut_heuristic[0]->get_heuristic());
	    max_h = max(max_h, lmcut_heuristic[0]->get_heuristic());
	    h_initial_v.push_back(lmcut_heuristic[0]->get_heuristic());            	
	  }
        }
#endif 

	//initial_value = min_h;
	initial_value = max_h;
	cout<<"probe:"<<call_number<<",search_time:"<<search_time()<<",initial min_h;"<<min_h<<",initial max_h:"<<max_h<<",max_h:"<<max_h<<endl;
	if(call_number%100==0){
	  cout<<",search_time:"<<search_time()<<endl;
	}
	else{
	  cout<<endl;
	}	

        cout<<"f_boundary = "<<f_boundary<<endl;
        if (f_boundary) {
	   threshold = f_boundary;
	   cout<<"threshold set manually to "<<threshold<<endl;
	} else {
	   threshold = 2 * max_h;
	   cout<<"default threshold twice max_h:"<<threshold<<endl;
	   /*  if(threshold==0){
	     threshold=2;
	     cout<<"max_h is 0 so we arbitrarily set the prediction f value to twice the initial value"<<endl;
	   }*/
	}
        cout<<"threshold = "<<threshold<<endl;

        //const GlobalState &initial_state = g_initial_state();
        int max_h_initial_value = 0;
	double sum_all_h_initial_values = 0;
        for (size_t i = 0; i < h_initial_v.size(); i++) {
            int a = h_initial_v.at(i);
	    sum_all_h_initial_values += a;
            if (a > max_h_initial_value) {
                max_h_initial_value = a;
            }
        }
	cout<<"max_h_initial_value:"<<max_h_initial_value<<endl;
        
        cout<<"\tthreshold: "<<threshold<<endl;
        //cout<<"\nprint h_initial_v\n";
	for (size_t i = 0; i < h_initial_v.size(); i++) {
            int h_value = h_initial_v.at(i);
            //cout<<h_value;
            //if (i != h_initial_v.size() - 1) {
               //cout<<"/";
            //}
            if (h_value <= threshold) {
               b_initial_v.set(i);
            }
        }
	cout<<"b_initial_v,size()"<<b_initial_v.size()<<endl;fflush(NULL);
         
        std::vector<const GlobalOperator*> applicable_ops0;
        GlobalState global_state0 = g_state_registry->lookup_state(g_initial_state().get_id()); 
        g_successor_generator->generate_applicable_ops(global_state0, applicable_ops0);

        //count nodes generated
        double amount_initial = (double)applicable_ops0.size();


	collector.insert(std::pair<boost::dynamic_bitset<>, double>(b_initial_v, 1 + amount_initial));
	
	collector_heur.insert(std::pair<std::vector<int>, double>(h_initial_v, 1 + amount_initial));
	//collector_heur.insert(std::pair<std::vector<int>, double>(h_initial_v, sum_all_h_initial_values));
        //cout<<"\n";
        
        const GlobalState &initial_state2 = g_initial_state();
        StateID initial_state_id = initial_state2.get_id();
	
	SSNode node;
#ifdef _LMCUT_EARLY_TERM
        if(lmcut_heuristic.size()>0){
	  node.set_lmcut_expanding(true);
	}
#endif
	node.setId(initial_state_id);
        node.setWeight(1.0);
        node.setGreal(0);  //g_real value of the node
        //node.setHC(h_initial_v);
        node.setH(min_h);
	/*
	 * Seeding the prediction with the children of the start state
	 *
	 */

	//Only evaluating those heuristics which still have not pruned the current path
	Type type = sampler->getType(node.get_id(), initial_value, 1);

	type.setLevel( 0 ); // level where the node is located

	queue.insert( pair<Type, SSNode>( type, node ) );

        int nraiz = 0;
  
	long queue_counter=0;
	while( !queue.empty() )
	{
	  queue_counter++;
	  if(queue_counter%1000==0){
	    if(search_time()>sampling_time_limit||g_timer()>overall_time_limit){
	      cout<<"Search_timer past maximum sampling_time"<<endl;
	      cout<<"selecting best heuristic after search_time: "<<search_time()<<", seconds,g_timer:"<<g_timer()<<endl;
	      runReports(true); 
	      exit(0);
	    }
	  }
#ifdef _SS_DEBUG
	  cout<<"queue.size:"<<queue.size()<<",search_time:"<<search_time<<endl;
#endif
		Type out = queue.begin()->first;
		SSNode s = queue.begin()->second;

               	int g_real =  s.getGreal();
                int level = out.getLevel();
		double w = s.getWeight(); 
                //std::vector<int> h_global_v = s.getHC();
		int min_h=s.getH();


		/*  
                boost::dynamic_bitset<> b_raiz_v(heuristics.size());
                std::vector<int> f_global_v;

                for (size_t i = 0; i < h_global_v.size(); i++) {
                    int h_value = h_global_v.at(i);
                    f_global_v.push_back(h_value + g_real);
		}

                for (size_t i = 0; i < f_global_v.size(); i++) {
                    int f_value = f_global_v.at(i);
                    if (f_value <= threshold) {
                     	b_raiz_v.set(i);
                    }
                }*/
                //printQueue();

                std::map<Type, SSNode>::iterator rt;
                rt = queue.find(out);


		queue.erase( rt );
                   
                nraiz++;                
#ifdef _SS_DEBUG
		  cout<<"Raiz: "<<nraiz<<" h  = "<<min_h;fflush(stdout);
	
		  cout<<", g_real = "<<g_real<<", f = ";fflush(stdout);
		  cout<<min_h + g_real;fflush(stdout);
		  cout<<", level = "<<level;fflush(stdout);
		  cout<<", w = "<<w<<"\n";fflush(stdout);
#endif

		vweight.push_back(s.getWeight());
		//add the SSNode and Type
                SSQueue smean;
                smean.setNode(s);
                smean.setT(out);
                vmeanheur.push_back(smean);
	        
                //Insert each node.
                //Node2 node(getMinHeur(h_global_v) + g_real, level);
                Node2 node(min_h + g_real, level);
                //node.setFs(f_global_v);
                //node.setL(level);
                //count nodes expanded
                if ( (min_h + g_real) <= threshold) {
			std::pair<std::map<Node2, double>::iterator, bool> ret0;

                	std::map<Node2, double>::iterator it0;

                	ret0 = expanded.insert(std::pair<Node2, double>(node, s.getWeight()));
                	it0 = ret0.first;

                	if (ret0.second) {
                    		//cout<<"new node expanded is added."<<endl;
                	} else {
                    		//cout<<"node expanded is being updated."<<endl;
                    		it0->second += s.getWeight();
                    		//cout<<"it0->second = "<<it0->second<<endl;
                	}
                }
		else{//Nodes could be added through BFS, need to update F bound accordingly
#ifdef _SS_DEBUG
		  int prev_f_bound=next_f_bound;
#endif
		  next_f_bound=min(next_f_bound,min_h + g_real);
#ifdef _SS_DEBUG
		  if(next_f_bound!=prev_f_bound){
		    cout<<"next_f_bound:"<<next_f_bound<<",prev_f_bound:"<<prev_f_bound<<endl;
		  }
#endif
		}
                //end count nodes expanded
 
                std::vector<const GlobalOperator*> applicable_ops;
                GlobalState global_state = g_state_registry->lookup_state(s.get_id()); 
		g_successor_generator->generate_applicable_ops(global_state, applicable_ops);

                //count nodes generated
                double amount = (double)applicable_ops.size();

                std::pair<std::map<Node2, double>::iterator, bool> ret;
                std::map<Node2, double>::iterator it;

                ret = generated.insert(std::pair<Node2, double>(node, amount*w));
                it = ret.first;


                if (ret.second) {
                   //cout<<"new node is added."<<endl;
                } else {
                   //cout<<"old is being updated."<<endl;
                   it->second += amount*w;
                   //cout<<"new = "<<it->second<<endl;
                }
                 //end count nodes generated
#ifdef _SS_DEBUG
		  cout<<"\t_____________________begin Childs________________________\n";fflush(stdout);
#endif
                int h =  INT_MAX/2;
		L.clear();
		check.clear();
		for (size_t i = 0; i < applicable_ops.size(); ++i)
		{
                        const GlobalOperator *op = applicable_ops[i];
                        GlobalState child = g_state_registry->get_successor_state(global_state, *op);

			vector<int> h_child_v;
                  	boost::dynamic_bitset<> b_child_v(heuristics.size()+lmcut_heuristic.size());b_child_v.set();
                  	vector<int> F_culprit;
#ifdef _SS_DEBUG
			cout<<"initial b_child_v:"<<b_child_v<<endl;fflush(NULL);
#endif
                  	
                  	//vector<string> heur_name_v;
			//int max_h=0;
			
                  	for (size_t i = 0; i < heuristics.size(); ++i){
                      		heuristics[i]->evaluate(child);
				if(heuristics[i]->is_dead_end()){
				  b_child_v.reset(i);
				  if(domination_check){
				    F_culprit.push_back(INT_MAX/2);
				  }
				  //max_h=INT_MAX/2;
				  //cout<<"dead_end found"<<endl;
				  h_child_v.push_back(INT_MAX/2);
				}
				else{
				  h = min(h, heuristics[i]->get_heuristic());
				  //max_h = max(h, heuristics[*it]->get_heuristic());
#ifdef _SS_DEBUG
				  //cout<<"h("<<i<<"):"<<heuristics[i]->get_heuristic()<<",f:"<< heuristics[i]->get_heuristic() + g_real + get_adjusted_cost(*op)<<",thresh:"<<threshold<<endl;
				  int prev_f_bound=next_f_bound;
#endif
				  if ( heuristics[i]->get_heuristic() + g_real + get_adjusted_cost(*op)  > threshold) {
				    next_f_bound=min(next_f_bound,heuristics[i]->get_heuristic() + g_real + get_adjusted_cost(*op));
#ifdef _SS_DEBUG
				    if(next_f_bound!=prev_f_bound){
				      cout<<"prev_f_bound:"<<prev_f_bound<<",next_f_bound:"<<next_f_bound<<endl;
				    }
#endif
				  //cout<<"h("<<i<<"):"<<heuristics[i]->get_heuristic()<<",f:"<< heuristics[i]->get_heuristic() + g_real + get_adjusted_cost(*op)<<",thresh:"<<threshold<<endl;
				    b_child_v.reset(i);
 	
					//Record F values until domination check finished
				  } 
				  int new_heur = heuristics[i]->get_heuristic() + g_real + get_adjusted_cost(*op);
				  if(domination_check){
				    F_culprit.push_back(new_heur);
				  }
				  h_child_v.push_back(new_heur);
				  //string heur_name = heuristics[i]->get_heur_name();
				  //heur_name_v.push_back(heur_name);
				}
			}
#ifdef _LMCUT_EARLY_TERM
	bool lmcut_expanding=true;
        if(lmcut_heuristic.size()>0&&s.get_lmcut_expanding()){
	  lmcut_heuristic[0]->evaluate(child);
	  if(lmcut_heuristic[0]->is_dead_end()){
	    b_child_v.reset(heuristics.size());
	    lmcut_expanding=false;
	  }
	  else if ( lmcut_heuristic[0]->get_heuristic() + g_real + get_adjusted_cost(*op)  > threshold) {
	    lmcut_expanding=false;
	    b_child_v.reset(heuristics.size());
	    next_f_bound=min(next_f_bound,lmcut_heuristic[0]->get_heuristic() + g_real + get_adjusted_cost(*op));
	  }
	}
	else{//automatically mark as prunning, ancestor already pruned path with LMCUT
	    b_child_v.reset(heuristics.size());
	}

	//sum heuristics values
        if(lmcut_heuristic.size()>0){
	  lmcut_heuristic[0]->evaluate(child);
	  if(lmcut_heuristic[0]->is_dead_end()){
	    h_child_v.push_back(INT_MAX/2);
	  }
	  else if ( lmcut_heuristic[0]->get_heuristic() + g_real + get_adjusted_cost(*op)  > threshold) {
	    h_child_v.push_back(lmcut_heuristic[0]->get_heuristic() + g_real + get_adjusted_cost(*op));
	  }
	}
	else{//automatically mark as prunning, ancestor already pruned path with LMCUT
	    h_child_v.clear();
	}
#endif

#ifdef _SS_DEBUG
	  if(next_f_bound!=prev_f_bound){
	    cout<<"prev_f_bound:"<<prev_f_bound<<",next_f_bound:"<<next_f_bound<<endl;
	  }
#endif
#ifdef _SS_DEBUG
			  cout<<", g_real = "<<g_real + get_adjusted_cost(*op)<<" f_min = "<< h + g_real + get_adjusted_cost(*op)<<",b_child_v.count:"<<b_child_v.count()<<endl;
			  cout<<"\tget_adjusted_cost(*op) = "<<get_adjusted_cost(*op)<<"\n";
			  cout<<"\tChild_"<<(i+1)<<" : h = "<<h<<",b_child_v:"<<b_child_v<<endl;
			 
			  cout<<h + g_real + get_adjusted_cost(*op)<<endl;
			  cout<<", level = "<<(level + 1);
			  cout<<", w = "<<w<<"\n";
#endif

	     		std::vector<const GlobalOperator *> applicable_ops_2;
             

             		GlobalState global_state_2 = g_state_registry->lookup_state(child.get_id());
			//cout<<"S:"<<endl;global_state_2.dump_inline();
             		g_successor_generator->generate_applicable_ops(global_state_2, applicable_ops_2);
             
             		int amount = applicable_ops_2.size();
			//number of nodes expanded                
             		std::pair<std::map<boost::dynamic_bitset<>, double>::iterator, bool> ret2; 
             		std::map<boost::dynamic_bitset<>, double>::iterator it2;
			//sum the heuristic values
			std::pair<std::map<vector<int>, double>::iterator, bool> ret4;
			std::map<vector<int>, double>::iterator it4;

			double sum_all_h_child_values = 0;
                        for (size_t i = 0; i < h_child_v.size(); i++) {
                            int h_value = h_child_v.at(i);
			    sum_all_h_child_values += h_value;
			    //cout<<h_value + g_real + get_adjusted_cost(*op);
 		            /*if (i != h_child_v.size() -1) {
                         	cout<<"/";
                      	    }*/
                        }
			if(domination_check){
			  //cout<<"inserting F_culprit:"<<F_culprit<<endl;
			  F_culprits.insert(F_culprit);
			}
      
                        //add to the collector
             		if (b_child_v.count() > 0) {
 	     			ret2 = collector.insert(std::pair<boost::dynamic_bitset<>, double>(b_child_v, amount*w));
             			it2 = ret2.first;
				//predict the number of nodes expanded in the search tree.
             			if (ret2.second) {
					//cout<<"raiz bc new is added"<<endl;
             			} else {
                			//cout<<"raiz bc old is being updated"<<endl; 
                			it2->second += amount*w;
                			//cout<<", newcc : "<<it2->second<<"\n"; 
             			}
				//sum the heuristic values in the search tree.
				ret4 = collector_heur.insert(std::pair<vector<int>, double>(h_child_v, amount*w));
				//ret4 = collector_heur.insert(std::pair<vector<int>, double>(h_child_v, sum_all_h_child_values*w));
				it4 = ret4.first;

				if (ret4.second) {

				} else {
					it4->second += amount*w;
					//it4->second += sum_all_h_child_values*w;
				}

                        //Make pruning
			   Type object = sampler->getType(child.get_id(), h, 1);
			   
                           object.setLevel( level + 1 );
                           
                           SSNode child_node;
#ifdef _LMCUT_EARLY_TERM
			   if(lmcut_expanding){
			     child_node.set_lmcut_expanding();
			   }
			   else{
			     child_node.set_lmcut_expanding(false);
			   }
#endif
                           StateID child_state_id = child.get_id();
                           child_node.setId(child_state_id);
                           child_node.setWeight(w);
		           child_node.setGreal(g_real + get_adjusted_cost(*op)); 
 			   //child_node.setHC(h_child_v);
 			   child_node.setH(h);

				
#ifdef _SS_DEBUG
			     cout<<"\t\tChild f<=threshold: h = "<<h; 
			     cout<<", g_real = "<<g_real + get_adjusted_cost(*op)<<" f = ";
				 cout<<h + g_real  +  get_adjusted_cost(*op);
			     cout<<", level = "<<level + 1<<"\n";
#endif
			   if (get_adjusted_cost(*op) == 0) {
#ifdef _SS_DEBUG
				cout<<"\t\tget_adjusted_cost(*op) == 0\n";fflush(NULL);
#endif
			   	BFS(child_node, object);
#ifdef _SS_DEBUG
				cout<<"after BFS, L size:"<<L.size()<<endl;
#endif
				std::set<SSQueue>::iterator it;
				for (it = L.begin(); it != L.end(); it++) {	
					SSQueue sst = *it;
					SSNode node = sst.getNode();
					Type t = sst.getT();
					//double new_g_real = node.getGreal();
					StateID new_state_id = node.get_id();
					double w2 = node.getWeight();

					//cout<<"\n\t\tNode restored: h = "<<t.getH()<<", g_real = "<<node.getGreal()<<", f = "<<t.getH() + node.getGreal()<<", level = "<<t.getLevel()<<", w = "<<w2<<"\n";
                        		GlobalState global_state2 = g_state_registry->lookup_state(new_state_id);
                           			
					map<Type, SSNode>::iterator queueIt = queue.find( t );
#ifdef _SS_DEBUG
					cout<<"state_id:"<<new_state_id<<endl;
#endif
			   		if( queueIt != queue.end() )
			   		{
                                	SSNode snode = queueIt->second;

#ifdef _SS_DEBUG
                                		cout<<"\t\t\tzc: The duplicate node is: h = "<<queueIt->first.getH()<<", g = "<<snode.getGreal()<<", f = "<< queueIt->first.getH() + snode.getGreal()<<", w = "<<snode.getWeight()<<", level = "<<queueIt->first.getLevel()<<"\n";
#endif
                                
						double wa = (double)snode.getWeight();
						//snode.setWeight( wa + w);
                                		queueIt->second.setWeight(wa + w2); // set w the node that already exists
                                		//cout<<"\t\t\tzc: before ss process starts, the w of the duplicate node is updated to: "<<queueIt->second.getWeight()<<endl; 
                                		//std::pair<std::map<Type, SSNode>::iterator, bool> ret0;
                                		//ret0 = queue.insert(pair<Type, SSNode>(object, snode));
                                		//cout<<"\tsnode.getWeight() = "<<snode.getWeight()<<endl;
                                		//queueIt->second.setWeight(snode.getWeight());
						double prob = ( double )w2 / (double)( wa + w2);
						int rand_100 =  RanGen2->IRandom(0, 99);  //(int)g_rng.next(100);
                          	 
                                		double a = (( double )rand_100) / 100;
                                		//cout<<"a = "<<a<<" prob = "<<prob<<endl;
                                
						if (a < prob) 
						{
#ifdef _SS_DEBUG
                                        		cout<<"\t\t\tzc: Added even though is duplicate.\n";                               
#endif
				        		node.setWeight( wa + w2);
#ifdef _SS_DEBUG
                                        		cout<<"\t\t\tzc: the w is updated to = "<<node.getWeight()<<endl;
#endif
                                        		std::pair<std::map<Type, SSNode>::iterator, bool> ret3;
                                     			queue.erase(t); 
                                        
                                        		ret3 = queue.insert( pair<Type, SSNode>( t, node ));      
                                        		queueIt = ret3.first;
                                        		queueIt->second.setWeight(node.getWeight());
						} else {
#ifdef _SS_DEBUG
                                        		cout<<"\t\t\tzc: Not added.\n";
                                        		cout<<"\t\t\tbut the w is updated for the node that already exists to: "<<queueIt->second.getWeight()<<endl;
#endif
                                		}
			   		} 
			   		else
			   		{
#ifdef _SS_DEBUG
                                		cout<<"\t\t\tzc: New L node added.\n";
#endif
						queue.insert( pair<Type, SSNode>( t, node ) );
                                		//cout<<"\t\tsucc_node2.getWeight() = "<<succ_node2.getWeight()<<"\n";
                                
                                		//cout<<"\t\t\tzc: Child: h = "<< t.getH() <<", g_real = "<< new_g_real <<", f = "<< t.getH() + new_g_real << " threshold: " << threshold <<" w = "<<node.getWeight()<<endl;
                           		}// End queueIt != queue.end()
				}   //End for set lopp
#ifdef _SS_DEBUG
				cout<<"Finished with L"<<endl;
#endif
			   } else {

                                                     
                           map<Type, SSNode>::iterator queueIt = queue.find( object );
			   if( queueIt != queue.end() )
			   {
       	                        SSNode snode = queueIt->second;

#ifdef _SS_DEBUG
				cout<<"\t\tThe duplicate node is: h = "<<h;
				    
				cout<<", g_real = "<<g_real + get_adjusted_cost(*op)<<" f = ";
				cout<<h + g_real  +  get_adjusted_cost(*op);
				cout<<", w = "<<snode.getWeight();
				cout<<", level = "<<level + 1<<"\n";
#endif

				double wa = (double)snode.getWeight();
				//snode.setWeight( wa + w);
                                queueIt->second.setWeight(wa + w);
#ifdef _SS_DEBUG
				  cout<<"\t\tbefore ss process starts, the w of the duplicate node is updated to: "<<queueIt->second.getWeight()<<endl;
#endif
                                //std::pair<std::map<Type, SSNode>::iterator, bool> ret0;

                                //ret0 = queue.insert(pair<Type, SSNode>(object, snode));
                                //cout<<"\tsnode.getWeight() = "<<snode.getWeight()<<endl;
                                //queueIt->second.setWeight(snode.getWeight());
 
 
				double prob = ( double )w / (double)( wa + w );
				int rand_100 =  RanGen2->IRandom(0, 99);  //(int)g_rng.next(100);
                          	 
                                double a = (( double )rand_100) / 100;
                                //cout<<"a = "<<a<<" prob = "<<prob<<endl; 
                                
				if( (a < prob))
				{
#ifdef _SS_DEBUG
				    cout<<"\t\tAdded even though is duplicate.\n";
#endif
                                        
				        child_node.setWeight( wa + w);
#ifdef _SS_DEBUG
					  cout<<"\t\tthe w is updated to = "<<child_node.getWeight()<<endl;
#endif
                                        std::pair<std::map<Type, SSNode>::iterator, bool> ret;
                                     	queue.erase(object); 

                                        ret = queue.insert( pair<Type, SSNode>( object, child_node ));      

                                        queueIt = ret.first;
                                        queueIt->second.setWeight(child_node.getWeight());
                                        
                                      	
				} else {
#ifdef _SS_DEBUG
					  cout<<"\t\tNot added.\n";
					  cout<<"\t\tbut the w is updated for the node that already exists to: "<<queueIt->second.getWeight()<<endl;
#endif
                                }
			   } 
			   else
			   {
#ifdef _SS_DEBUG
			       cout<<"\t\tNew node added\n";
			       //Now update the non-prunning set of heuristics for the node
#endif
				    queue.insert( pair<Type, SSNode>( object, child_node ) );
			       }
			   }
                        }
			else 
			{
#ifdef _SS_DEBUG
				cout << "\tNode was pruned!" << endl;
#endif
			}
#ifdef _SS_DEBUG
                        cout<<"\tend Child_"<<(i+1)<<"\n";
#endif
		}
#ifdef _SS_DEBUG
		  cout<<"\t______________________end Childs_____________________\n";
#endif

	}
	//cout<<"end queue"<<endl;
}

int SSSearch::getMinHeur(vector<int> v) {
	int h_min = INT_MAX/2;
	for (size_t i = 0; i < v.size(); i++) {
		int aux = v.at(i);
		if (h_min > aux) {
			h_min = aux;
		}
	}
	return h_min;
}


int SSSearch::getMaxHeur(vector<int> v) {
	int h_max = 0;
	//cout<<"print v\n";
	for (size_t i = 0; i < v.size(); i++) {
		int aux = v.at(i);
		//cout<<aux<<", ";
		if (h_max < aux) {
			h_max = aux;
		}
	}
	//cout<<"\n";
	return h_max;
}

void SSSearch::runReports(bool cmd) {
	if (run_method == "grhs") {
		select_random_greedy(cmd);
	} else if (run_method == "sscc") {
	  	generateSSCCReport(cmd);
        } else if (run_method == "both") {
	  	select_random_greedy(cmd);	
	  	generateSSCCReport(cmd);
	}
	generateGeneratedReport();
        generateExpandedReport(false);
	///select_best_heuristics_greedy();
}

void SSSearch::generateGeneratedReport() {
       double count_nodes = 0;
       for (map<Node2, double>::iterator iter = generated.begin(); iter != generated.end(); iter++) {
           double n = iter->second;
           count_nodes += n;
       }
       number_nodes_generated = count_nodes/(double)ss_probes;
       cout<<"count nodes generates : "<<number_nodes_generated<<endl;
       generated.clear();
}

void SSSearch::generateExpandedReport(bool create_report) {
        double total_nodes = 0.0;
        for (map<Node2, double>::iterator iter = expanded.begin(); iter != expanded.end(); iter++) {
            double n = iter->second;
            total_nodes += n;
        }
	double predictionExpanded = (double)total_nodes/(double)n_probes_global;
	number_nodes_expanded = predictionExpanded; 
        cout<<"count nodes expanded : "<<predictionExpanded <<endl;

	if (create_report) {
		string dirDomain, dirfDist, outputFile;

		if (is_mov_bound) {
                	string nameProbes = "reportss_bounds";
                	nameProbes += "_probes_";
                	nameProbes += boost::lexical_cast<std::string>(ss_probes);
                	cout<<"nameProbes = "<<nameProbes<<"\n";
                	dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes  +"/"+dominio_global;
                	dirfDist = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes  +"/"+dominio_global+"/fdist";
                	outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes  +"/"+dominio_global+"/fdist/"+tarefa_global;
		} else {
			string nameProbes = "reportss_";
                	nameProbes += boost::lexical_cast<std::string>(ss_probes);
                	nameProbes += "_probes";
                	cout<<"nameProbes = "<<nameProbes<<"\n";

        		dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global;
        		dirfDist = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/fdist";
        		outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/fdist/"+tarefa_global;
		}

        	ofstream output;

        	output.open(outputFile.c_str());
        	output<<"\t"<<outputFile.c_str()<<"\n";
        	//output<<"predictionSS: "<<totalPrediction<<"\n";
        	output<<"predictionSS: "<<predictionExpanded<<"\n";fflush(stdout);
        	output<<"ss_timer: "<<ss_timer_value<<"\n";

        	if (system(dirDomain.c_str())) {
           		cout<<"Directory: "<<dirDomain.c_str()<<" created."<<endl;
        	}

        	if (system(dirfDist.c_str())) {
           		cout<<"Directory: "<<dirfDist.c_str()<<" created."<<endl;
        	}
        	cout<<"print."<<endl;
        	for (int i = 0; i <= threshold; i++) {
            	int k = 0;
            	vector<long> f;
            	vector<double> q;
            	for (map<Node2, double>::iterator iter = expanded.begin(); iter != expanded.end(); iter++) {
                 	Node2 n = iter->first;
                
                 	if (i == n.getL()) {
                    		k++;
                    		f.push_back(n.getF());
                    		q.push_back(((double)iter->second)/(double)n_probes_global);

                    		//cout<<"l = "<<n.getL()<<" f = "<<n.getF()<<" q = "<<(iter->second)/ss_probes<<endl;
                    		//output<<"l = "<<n.getL()<<" f = "<<n.getF()<<" q = "<<(iter->second)/ss_probes<<"\n";
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
            	cout<<"\n";
            	output<<"\n";
            
        	}
        	output.close();

	}
        expanded.clear();
}

void SSSearch::updateGlobalVariables() {
	cout<<"UPDATING_VARIABLES\n";
	n_probes_global = ss_probes;//the average must always be done with the ss_probes (Marvin's code the ss_probes is sent from the sh)
	
	n_heuristics_global = all_heuristics.size();
	dominio_global = domain_name;
        tarefa_global = problem_name2;
        heuristica_global = heuristic_name2;
	domain_pddl_global = domain_instance_pddl;
}

void SSSearch::updateSSCC() {
        size_t found = tarefa_global.find(".");
        string name = tarefa_global.substr(0, found);
        name+="_F_";
        name+=boost::lexical_cast<std::string>(threshold);
        name += ".csv";
	
	string dirDomain, dirSSCC, outputFile;
	if (is_mov_bound) {
		string nameProbes = "reportss_bounds";
                nameProbes += "_probes_";
                nameProbes += boost::lexical_cast<std::string>(n_probes_global);
                cout<<"nameProbes = "<<nameProbes<<"\n";

        	dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global;
        	dirSSCC = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/bc";
        	outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/bc/"+name;
	} else {
		string nameProbes = "reportss_";
                nameProbes += boost::lexical_cast<std::string>(n_probes_global);
                nameProbes += "_probes_sscc";
                cout<<"nameProbes = "<<nameProbes<<"\n";

		dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global;
        	dirSSCC = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/bc";
        	outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/bc/"+name;
	}

        if (system(dirDomain.c_str())) {
           cout<<"Directory: "<<dirDomain.c_str()<<" created."<<endl;
        }

        if (system(dirSSCC.c_str())) {
           cout<<"Directory: "<<dirSSCC.c_str()<<" created."<<endl;
        }

        ofstream output;
        output.open(outputFile.c_str());

	//print the name of the all_heuristics just to be analyzed later
        for (size_t i = 0; i < all_heuristics.size(); i++) {
            string heur_name = all_heuristics[i]->get_heur_name();
            output<<"h(,"<<i<<"):,"<<heur_name<<"\n";
        }

	n_heuristics_global = all_heuristics.size();
	count_line_sscc = collector.size();

	//get the combintation 1/0/1/0..../1/1
	harray_sscc = new int*[count_line_sscc];
	ccarray_sscc = new double*[count_line_sscc];
	for (int i = 0; i < count_line_sscc; i++) {
		harray_sscc[i] = new int[n_heuristics_global];
		ccarray_sscc[i] = new double[1];
	}

	int counter_line = 0;
        for (map<boost::dynamic_bitset<>, double>::iterator iter = collector.begin(); iter != collector.end(); iter++) {
                boost::dynamic_bitset<> b_node_v = iter->first;
                double cc = iter->second;
                //cout<<"bc(";
                output<<"bc(";
                for (size_t i = 0; i < b_node_v.size(); i++) {
                        //cout<<b_node_v.test(i);
                        output<<b_node_v.test(i);
			harray_sscc[counter_line][i] = b_node_v.test(i);
                        if (i != b_node_v.size() - 1) {
                                //cout<<"/";
                                output<<"/";
                        }
                }
                //cout<<")cc="<<(double)cc/(double)n_probes_global<<"\n";
		double result_cc = (double)cc/(double)n_probes_global; 
                //output<<")cc="<<(double)cc/(double)n_probes_global<<"\n";
                output<<")cc="<<result_cc<<"\n";
		ccarray_sscc[counter_line][0] = result_cc;
		counter_line++;
        }
        output.close();
}

void SSSearch::generateSSCCReport(bool termination) {
	updateSSCC();
	if (termination) {
		//make it work in 30 minutes
		string delimiter = ",";
		//cout<<"heuristic-information\n";
		map<string, double> add_line_map_heuristic;
		map<string, vector<string> > map_info_heur;
		for (int i = 0; i < n_heuristics_global; i++) {
			vector<string> collector;
                	string s =  all_heuristics[i]->get_heur_name();
                	string pot[6];
                	size_t pos = 0;
                	string token;
                	int index = 0;
                	while ((pos = s.find(delimiter)) != std::string::npos) {
                        	token = s.substr(0, pos);
				//cout<<"token="<<token<<"\n";
                        	pot[index] = token;
                        	s.erase(0, pos + delimiter.length());
                        	index++;
                	}
                	pot[index] = s;

           		string heuristic_name_created = pot[0],
	    		number_h = std::to_string(i), //consider this order because SS commands
			mutation_rate,
			mutation_rate_aux,
			size_gapdb,
			size_gapdb_aux,
			wd,
			wd_aux,
			name;

	    		if (heuristic_name_created == "ipdb") {
            			name = "ipdb_" + number_h;
            		} else if (heuristic_name_created == "lmcut") {
            			name = "lmcut_" + number_h;
            		} else if (heuristic_name_created == "merge_and_shrink") {
            			name = "mands_" + number_h;
            		} else {
            			name = "gapdb_" + number_h;
            		}
			//cout<<"name="<<name<<"\n";
			mutation_rate_aux = pot[1];
                	size_t t2 = mutation_rate_aux.find(":");
                	mutation_rate = mutation_rate_aux.substr(t2 + 1, mutation_rate_aux.length());
                	//cout<<"mutation_rate = "<<mutation_rate<<"\n";

			size_gapdb_aux = pot[2];
                	size_t t3 = size_gapdb_aux.find("=");
                	size_gapdb = size_gapdb_aux.substr(t3 + 1);
                	//cout<<"size_gapdb = "<<size_gapdb<<"\n";

                	//without disjoint patterns
                	wd_aux = pot[3];
                	//cout<<"wd_aux = "<<wd_aux<<"\n";
                	size_t t4 = wd_aux.find("out");
                	if (t4 < 100) {
                        	wd = "false";
                	} else {
                        	wd = "true";
                	}
                	//cout<<"wd = "<<wd<<"\n\n";

			collector.push_back(number_h);
                	collector.push_back(mutation_rate);
                	collector.push_back(size_gapdb);
                	collector.push_back(wd);

                	map_info_heur.insert(pair<string, vector<string> >(name, collector));

			//number of nodes expanded by each heuristic
                	double sum_ones = 0;
                	for (int j = 0; j < count_line_sscc; j++) {
                        	if (harray_sscc[j][i] == 1) {
                                	sum_ones += ccarray_sscc[j][0];
                        	}
                	}
                	add_line_map_heuristic.insert(pair<string, double>(name, sum_ones));
        	}
		//cout<<"\n\n";

		//cout<<"printing m:\n";
        	map<string, double>::iterator iter_test;
        	double min_number_expanded =  std::numeric_limits<double>::max();
        	string min_number_heuristic;
		vector<string> number_gapdb_heurs;
        	for (iter_test = add_line_map_heuristic.begin(); iter_test != add_line_map_heuristic.end(); iter_test++) {
                	string s = iter_test->first;
                	double d = iter_test->second;
			number_gapdb_heurs.push_back(s);
                	//cout<<s<<", "<<d<<"\n";
                	if (min_number_expanded > d) {
                        	min_number_expanded = d;
                        	min_number_heuristic = s;
                	}
        	}
        	//cout<<"min_number_expanded = "<<min_number_expanded<<"\n";
        	//cout<<"min_number_heuristic = "<<min_number_heuristic<<"\n";
		//Implement the evaluation time selection for the meta-reasoning


		double min_eval_time = std::numeric_limits<double>::max();
		int index_min_eval_time=0;
		boost::dynamic_bitset<> b_comb(n_heuristics_global + lmcut_heuristic.size());
		for (int i = 0; i < n_heuristics_global; i++) {
			b_comb.reset();
			b_comb.set(i);
			double cost_heur = calculate_time_costs_specific(b_comb);
			cout<<"heuristic="<<i<<",b_comb.size()="<<b_comb.size()<<"\t"<<b_comb<<"\tcost_heur="<<cost_heur<<"\n";
			if (min_eval_time > cost_heur) {
				min_eval_time = cost_heur;
				index_min_eval_time = i;
			}
		}
		string min_eval_time_heur = getHeuristicInfo(index_min_eval_time);
		cout<<"min_eval_time="<<min_eval_time<<",index_min_eval_time="<<index_min_eval_time<<"\theuristic_name="<<min_eval_time_heur<<"\n";
		//end evaluation time for meta-reasoning

        	vector<string> v_gapdb_string;
		string heuristic_good = "gapdb_good";

		int counter_just_ga_heur = 0;
		int total_gapdb_heuristics = getTotalGAHeurs(number_gapdb_heurs);
        	map<string, vector<string> >::iterator iter;
        	for (iter = map_info_heur.begin(); iter != map_info_heur.end(); iter++) {
                	string gapdb_string; //the name of each heuristic, just remember that the fd only support gapdb and do not gapdb_deep or gapdb_good => both need to be changed to gapdb
                	string s = iter->first;
                	vector<string> info = iter->second;

			if (run_min_eval_time_approach) {
				if (s == min_eval_time_heur) {
					gapdb_string = processHeuristicProperties(s, info, heuristic_good);
                                	v_gapdb_string.push_back(gapdb_string);	
                        		//cout<<"gapdb_string = "<<gapdb_string<<"\n";
				}// s == min_eval_time_heur
			} else {
				if (run_min_heuristic) {
                			if (s == min_number_heuristic) {
                				gapdb_string = processHeuristicProperties(s, info, heuristic_good);
                                		v_gapdb_string.push_back(gapdb_string);		
                        			//cout<<"gapdb_string = "<<gapdb_string<<"\n";
                			}// s == min_number_heuristic
				} else { //end run_min_heuristic
					if (isGAPDB(s)) {
						string all_heur_gapdb_good = "gapdb";	
						gapdb_string = processAllHeuristicProperties(s, info, all_heur_gapdb_good, total_gapdb_heuristics, counter_just_ga_heur);
                                		v_gapdb_string.push_back(gapdb_string);	
                        			//cout<<"gapdb_string = "<<gapdb_string<<"\n";
					}//end isGAPDB
				}//end run_min_heuristic
			} //end run_min_eval_time_approach
        	}//end loop map_info_heur
        	//cout<<"v_gapdb_string.size() = "<<v_gapdb_string.size()<<"\n";
        	//end astar_gpdb call the bc from ss

		string PROB_GOOD = "problemas_";
                PROB_GOOD += boost::lexical_cast<std::string>(n_probes_global);
                PROB_GOOD += "_probes_sscc";
                //cout<<"PROB_GOOD = "<<PROB_GOOD<<"\n";
		int deep_F_boundary = threshold;
		string method = "sscc";
		//create directories, running individual heuristic and running all gapdb heuristics using the same heuristic_good
		mkdirAstar(method, heuristic_good, PROB_GOOD);

		if (run_min_eval_time_approach) {
			for (size_t i = 0; i < v_gapdb_string.size(); i++) {
                		//get the real name
                		string real_heur = v_gapdb_string.at(i);
                		string task = real_heur;
                		//cout<<"task = "<<task<<"\n";
                		//size_t found_task_deep = task.find("deep");
                		size_t found_task_good = task.find("good");
                		string final_real_heur, final_number_heur;
                		string delimiter = "_";
                		if (found_task_good > 1000) {
					string t0 = real_heur;
                        		size_t found_t0 = t0.find("_");
                        		string previous_real_heur = t0.substr(0, found_t0);			
					final_real_heur = previous_real_heur;
                        		//cout<<"previous_real_heur = "<<final_real_heur<<"\n";
                        		if (previous_real_heur == "mands()") {
                                		final_real_heur = "merge_and_shrink(shrink_strategy=shrink_bisimulation(max_states=100000,threshold=1,greedy=false),merge_strategy=merge_dfp())";
                                		//final_real_heur = "merge_and_shrink()";
                        		} else if (previous_real_heur == "ipdb()") {
                                		final_real_heur = "ipdb(max_time=200)";
                        		}

                        		//get the heuristic number
                        		string t1 = real_heur;
                        		size_t found_t1 = t1.find("_");
                        		final_number_heur = t1.substr(found_t1 + 1, t1.length());
                		} else {
                        		string s2 = real_heur;
                        		string pot[6];
                        		size_t pos = 0;
                        		string token;
                        		int index = 0;
                        		while ((pos = s2.find(delimiter)) != std::string::npos) {
                                		token = s2.substr(0, pos);
                                		pot[index] = token;
                                		s2.erase(0, pos + delimiter.length());
                                		index++;
                        		}
                        		//cout<<"index = "<<index<<"\n";
                        		pot[index] = s2;
                        		//remove deep from pot[1]
                        		string pot1 = pot[1];
                        		size_t found_pot1 = pot1.find("(");
                        		string new_pot1 = pot1.substr(found_pot1, pot1.length());
                        		//end remove deep from pot[1]

                        		final_real_heur = "gapdb" + new_pot1;
                        		final_number_heur = pot[2];
				}
				cout<<"final_real_heur = "<<final_real_heur<<"\n";
                		cout<<"final_number_heur = "<<final_number_heur<<"\n\n";

				//begin
                		string new_problem_name = tarefa_global; //--problem_name == problema.c_str()
                		string t = new_problem_name;
                		size_t found = t.find(".");
                		string new_problem_name_mod = t.substr(0, found);

				string prob_name_gapdb = new_problem_name_mod + "_gapdb_" + final_number_heur  + ".pddl";

				//end get real name
				//create the directory of the problemas_500_probes_good

				string dirSASPLAN = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan/";
                		if (system(dirSASPLAN.c_str())) {
                			cout<<"create directory "<<dirSASPLAN.c_str()<<"\n";
                		}

				string dirSASPLANDomain = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan/"+dominio_global;
                		if (system(dirSASPLANDomain.c_str())) {
                			cout<<"create directory "<<dirSASPLANDomain.c_str()<<"\n";
                		}

                		//creation of each sh file for the gapdb heuristic
                		string arquivo;
                		arquivo = new_problem_name_mod + "_gapdb_" + final_number_heur + ".sh";
                		arquivo = "/" + arquivo;
                		arquivo = dominio_global + arquivo;
                		arquivo = "astar/"+heuristic_good+"/" + PROB_GOOD  +  "/" + arquivo;
                		arquivo = "marvin/" + arquivo;
                		arquivo = "marvin/"+ arquivo;
                		arquivo =  _HOME_INFO+"/" + arquivo;
				cout<<"arquivo_1695="<<arquivo<<"\n";

				string plan_dir_file = "/plan/"+dominio_global+"/"+tarefa_global;
				executeQsub(arquivo, final_real_heur, heuristic_good, PROB_GOOD, prob_name_gapdb, deep_F_boundary, method, plan_dir_file, false);
			}//v_gapdb_string for loop
		} else {
			if (run_min_heuristic) {
	 			for (size_t i = 0; i < v_gapdb_string.size(); i++) {
                			//get the real name
                			string real_heur = v_gapdb_string.at(i);
                			string task = real_heur;
                			//cout<<"task = "<<task<<"\n";
                			//size_t found_task_deep = task.find("deep");
                			size_t found_task_good = task.find("good");
                			string final_real_heur, final_number_heur;
                			string delimiter = "_";
                			if (found_task_good > 1000) {
						string t0 = real_heur;
                        			size_t found_t0 = t0.find("_");
                        			string previous_real_heur = t0.substr(0, found_t0);			
						final_real_heur = previous_real_heur;
                        			//cout<<"previous_real_heur = "<<final_real_heur<<"\n";
                        			if (previous_real_heur == "mands()") {
                                			final_real_heur = "merge_and_shrink(shrink_strategy=shrink_bisimulation(max_states=100000,threshold=1,greedy=false),merge_strategy=merge_dfp())";
                                			//final_real_heur = "merge_and_shrink()";
                        			} else if (previous_real_heur == "ipdb()") {
                                			final_real_heur = "ipdb(max_time=200)";
                        			}

                        			//get the heuristic number
                        			string t1 = real_heur;
                        			size_t found_t1 = t1.find("_");
                        			final_number_heur = t1.substr(found_t1 + 1, t1.length());
                			} else {
                        			string s2 = real_heur;
                        			string pot[6];
                        			size_t pos = 0;
                        			string token;
                        			int index = 0;
                        			while ((pos = s2.find(delimiter)) != std::string::npos) {
                                			token = s2.substr(0, pos);
                                			pot[index] = token;
                                			s2.erase(0, pos + delimiter.length());
                                			index++;
                        			}
                        			//cout<<"index = "<<index<<"\n";
                        			pot[index] = s2;
                        			//remove deep from pot[1]
                        			string pot1 = pot[1];
                        			size_t found_pot1 = pot1.find("(");
                        			string new_pot1 = pot1.substr(found_pot1, pot1.length());
                        			//end remove deep from pot[1]

                        			final_real_heur = "gapdb" + new_pot1;
                        			final_number_heur = pot[2];
					}
					//cout<<"final_real_heur = "<<final_real_heur<<"\n";
                			//cout<<"final_number_heur = "<<final_number_heur<<"\n\n";

					//begin
                			string new_problem_name = tarefa_global; //--problem_name == problema.c_str()
                			string t = new_problem_name;
                			size_t found = t.find(".");
                			string new_problem_name_mod = t.substr(0, found);

					string prob_name_gapdb = new_problem_name_mod + "_gapdb_" + final_number_heur  + ".pddl";

					//end get real name
					//create the directory of the problemas_500_probes_good

					string dirSASPLAN = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan/";
                			if (system(dirSASPLAN.c_str())) {
                				cout<<"create directory "<<dirSASPLAN.c_str()<<"\n";
                			}

					string dirSASPLANDomain = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan/"+dominio_global;
                			if (system(dirSASPLANDomain.c_str())) {
                				cout<<"create directory "<<dirSASPLANDomain.c_str()<<"\n";
                			}

                			//creation of each sh file for the gapdb heuristic
                			string arquivo;
                			arquivo = new_problem_name_mod + "_gapdb_" + final_number_heur + ".sh";
                			arquivo = "/" + arquivo;
                			arquivo = dominio_global + arquivo;
                			arquivo = "astar/"+heuristic_good+"/" + PROB_GOOD  +  "/" + arquivo;
                			arquivo = "marvin/" + arquivo;
                			arquivo = "marvin/"+ arquivo;
                			arquivo =  _HOME_INFO+"/" + arquivo;
					cout<<"arquivo="<<arquivo<<"\n";

					string plan_dir_file = "/plan/"+dominio_global+"/"+tarefa_global;
					executeQsub(arquivo, final_real_heur, heuristic_good, PROB_GOOD, prob_name_gapdb, deep_F_boundary, method, plan_dir_file, false);
				}//v_gapdb_string for loop
			} else { //end else run_min_heuristic
				string dirSASPLAN = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan_"+heuristic_good+"/";
                		if (system(dirSASPLAN.c_str())) {
                			cout<<"create directory "<<dirSASPLAN.c_str()<<"\n";
                		}

				string dirSASPLANDomain = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan_"+heuristic_good+"/"+dominio_global;
                		if (system(dirSASPLANDomain.c_str())) {
                			cout<<"create directory "<<dirSASPLANDomain.c_str()<<"\n";
                		}		

				string heuristic_generator;
        			for (size_t i = 0; i < v_gapdb_string.size(); i++) {
                			string heur = v_gapdb_string.at(i);
                			heuristic_generator += heur;
        			}
        			//heuristic_generator += ")";

        			cout<<"heuristic_genertor= "<<heuristic_generator<<"\n";
        			//begin
                		string new_problem_name = tarefa_global; //--problem_name == problema.c_str()
        			string t = new_problem_name;
        			size_t found = t.find(".");
        			string new_problem_name_mod = t.substr(0, found);

        			string prob_name_gapdb = new_problem_name_mod + "_gapdb_all.pddl";
        			//cout<<"prob_name_gapdb = "<<prob_name_gapdb<<"\n\n\n";
        			//end

        			//creation of each sh file for the gapdb heuristic
        			string arquivo;

        			arquivo = new_problem_name_mod + "_gapdb_all.sh";
				arquivo = "/" + arquivo;
        			arquivo = dominio_global + arquivo;
        			arquivo = "astar/"+heuristic_good+"/" + PROB_GOOD  +  "/" + arquivo;
        			arquivo = "marvin/" + arquivo;
        			arquivo = "marvin/"+ arquivo;
        			arquivo =  _HOME_INFO+"/" + arquivo;
				cout<<"arquivo="<<arquivo<<"\n";
				string plan_dir_file = "/plan_"+heuristic_good+"/"+dominio_global+"/"+tarefa_global;

				executeQsub(arquivo, heuristic_generator, heuristic_good, PROB_GOOD, prob_name_gapdb, deep_F_boundary, method, plan_dir_file, true);
			}//end run_min_heuristic
		}//end min_eval_time_approach
	}//end termination

	delete harray_sscc; //freed memory
	harray_sscc = NULL; //pointed dangling pt to NULL
	delete ccarray_sscc;
	ccarray_sscc = NULL;
}


string SSSearch::getHeuristicInfo(int index) {
	map<int, string> map_heuristic_name;
	string delimiter = ",";
	for (int i = 0; i < n_heuristics_global; i++) {
		vector<string> collector;
                string s =  all_heuristics[i]->get_heur_name();
                string pot[6];
                size_t pos = 0;
                string token;
                int index = 0;
                while ((pos = s.find(delimiter)) != std::string::npos) {
                	token = s.substr(0, pos);
			//cout<<"token="<<token<<"\n";
                        pot[index] = token;
                        s.erase(0, pos + delimiter.length());
                        index++;
                }
                pot[index] = s;

           	string heuristic_name_created = pot[0],
	    	number_h = std::to_string(i), //consider this order because SS commands	
		mutation_rate,
		mutation_rate_aux,
		size_gapdb,
		size_gapdb_aux,
		wd,
		wd_aux,
		name;

		if (false) {
	    		if (heuristic_name_created == "ipdb") {
            			name = number_h + "_ipdb";
            		} else if (heuristic_name_created == "lmcut") {
            			name = number_h + "_lmcut";
            		} else if (heuristic_name_created == "merge_and_shrink") {
            			name = number_h + "_mands";
            		} else {
            			name = number_h + "_gapdb";
            		}
		} else {
			if (heuristic_name_created == "ipdb") {
            			name = "ipdb_" + number_h;
            		} else if (heuristic_name_created == "lmcut") {
            			name = "lmcut_" + number_h;
            		} else if (heuristic_name_created == "merge_and_shrink") {
            			name = "mands_" + number_h;
            		} else {
            			name = "gapdb_" + number_h;
            		}
		}
		int index_heur = std::atoi(number_h.c_str());
		cout<<"pair:"<<index_heur<<","<<name<<"\n";
		map_heuristic_name.insert(pair<int, string>(index_heur, name));
        }
	string result;
	//get the name of the heuristic
	map<int, string>::iterator iter = map_heuristic_name.find(index);
	if (iter != map_heuristic_name.end()) {
		result = iter->second;
	} else {
		result = "no_name";
	}
	return result;
}

string SSSearch::processHeuristicProperties(string s, vector<string> info, string heuristic_good) {
	string result_string;
	string gapdb_string =  heuristic_good + "(mp=";
        //cout<<"heuristic (s) = "<<s<<"\n";
        //find the number
        string t = s;
        size_t found = t.find("_");
        string t_final = t.substr(found + 1, t.length());
        //cout<<"t_final = "<<t_final<<"\n";

        bool is_blind_heuristic = false;
        for (size_t i = 0; i < info.size(); i++) {
        	string parameter = info.at(i);
                //cout<<"\t"<<parameter<<"\n";
        	if (i == 1) {
        		gapdb_string += parameter;
                } else if (i == 2) {
                	gapdb_string += ",size="+parameter;
                	if (parameter == "") {
                        	is_blind_heuristic = true;
                        }
                } else if (i == 3) {
                	gapdb_string += ",disjoint="+parameter;
                }
	}
        gapdb_string+=")_" + t_final;
        //gapdb_string+=",eps=120,colls=5)";
        //cout<<"\tgapdb_string = "<<gapdb_string<<"\n\nif (is_blind_heuristic) {

	if (is_blind_heuristic) {
        	//Workaround
                string task2 = s;
		size_t found_task2 =  task2.find("_");
                string new_s = task2.substr(0, found_task2);

                string heur_blind = "blind()_" + t_final;
                if (new_s == "ipdb") {
                	heur_blind = "ipdb()_" + t_final;
                } else if (new_s == "lmcut") {
                	heur_blind = "lmcut()_" + t_final;
                } else if (new_s == "mands") {
                	heur_blind = "mands()_" + t_final;
                }
		result_string = heur_blind;
	} else {
		result_string = gapdb_string;
        }
        //cout<<"gapdb_string = "<<gapdb_string<<"\n";
	return result_string;
}


string SSSearch::processAllHeuristicProperties(string s, vector<string> info, string heuristic_good, int total_gapdb_heuristics, int &counter_just_ga_heur) {	
	string result_string;
	string gapdb_string =  heuristic_good +"(mp=";
	//cout<<"heuristic (s) = "<<s<<"\n";
        //find the number
        string t = s;
        size_t found = t.find("_");
        string t_final = t.substr(found + 1, t.length());
        //cout<<"t_final = "<<t_final<<"\n";

        bool is_blind_heuristic = false;
        for (size_t i = 0; i < info.size(); i++) {
        	string parameter = info.at(i);
                //cout<<"\t"<<parameter<<"\n";
                if (i == 1) {
                	gapdb_string += parameter;
                } else if (i == 2) {
                	gapdb_string += ",size="+parameter;
                        if (parameter == "") {
                        	is_blind_heuristic = true;
                        }
		} else if (i == 3) {
                	gapdb_string += ",disjoint="+parameter;
		}
	}
	gapdb_string+=")";//+ t_final;
        if (counter_just_ga_heur != total_gapdb_heuristics - 1) {
        	gapdb_string+=",";//+ t_final;
        }
        cout<<"counter_just_ga_heur = "<<counter_just_ga_heur<<"\n";
        counter_just_ga_heur++;

	if (is_blind_heuristic) {
        	//Workaround
                string task2 = s;

                size_t found_task2 =  task2.find("_");
                string new_s = task2.substr(0, found_task2);

                string heur_blind = "blind()_" + t_final;
                if (new_s == "ipdb") {
                	heur_blind = "ipdb()_" + t_final;
                } else if (new_s == "lmcut") {
                	heur_blind = "lmcut()_" + t_final;
                } else if (new_s == "mands") {
                	heur_blind = "mands()_" + t_final;
                }
		result_string = heur_blind;
	} else {
		result_string = gapdb_string;
	}
	cout<<"result_string="<<result_string<<"\n";
	return result_string;
        //cout<<"gapdb_string = "<<gapdb_string<<"\n";	
}

void SSSearch::mkdirAstar(string  method, string heuristic, string probLogs) {
	string dirProbGood = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic+"/";
        if (system(dirProbGood.c_str())) {
                cout<<"create directory "<<dirProbGood.c_str()<<"\n";
        }

        string dirProblema = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic+"/" + probLogs;
        if (system(dirProblema.c_str())) {
                cout<<"create directory "<<dirProblema.c_str()<<"\n";
        }

	string dirResultado = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic+"/reportastar_"+method;
        if (system(dirResultado.c_str())) {
        	cout<<"create directory "<<dirResultado.c_str()<<"\n";
        }

	string pastaProblema = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic+"/" + probLogs  + "/"+dominio_global;
	if (system(pastaProblema.c_str())) {
		cout<<"create directory "<<pastaProblema.c_str()<<"\n";
	}

	string pastaResultado = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic+"/" + probLogs  +  "/"+dominio_global+"/resultado";
	if (system(pastaResultado.c_str())) {
		cout<<"create directory "<<pastaResultado.c_str()<<"\n";
	}
}

void SSSearch::executeQsub(string arquivo, string final_real_heur, string heuristic_good, string PROB_GOOD, string prob_name_gapdb, int deep_F_boundary, string method, string plan_dir_file, bool apply_max) {
	string ASTAR_GOOD_NAME = "_SS_ASTAR";

	ofstream outfile(arquivo.c_str(), ios::out);
        string parameter =  final_real_heur;

	//Calculate the time to execute the process	
	//Timer good_time = 1800 - g_timer();
	double duration = 1800 - g_timer();
	string good_timer;
	ostringstream convert;
	convert<<duration;
	good_timer = convert.str();
	cout<<"search_time()="<<search_time()<<"\n";
	cout<<"g_timer()="<<g_timer()<<"\n";
	cout<<"duration="<<good_timer<<"\n";
	
	//Begin construction of the sh file
	outfile<<"#!/bin/bash\n\n";
        outfile<<"#PBS -N "<<ASTAR_GOOD_NAME<<"\n\n#PBS -m a\n\n#PBS -l walltime="<<good_timer<<"\n\n";
	outfile<<"#PBS -M marvin.zarate@ufv.br\n\n";
	outfile<<"source /usr/share/modules/init/bash\n\n";
	outfile<<"module load python\nmodule load mercurial\n\n";

	//cout<<"pasta = "<<dominio_global<<"\n\n;
	outfile<<"FD_ROOT="<<_HOME_INFO<<"/marvin"<<_FD_INFO<<"\n\n";
        outfile<<"TEMP="<<_HOME_INFO<<"/marvin"<<_FD_INFO<<"/temp\n\n";
        outfile<<"DIR=$(mktemp  --tmpdir=${TEMP})\n\n";	
        outfile<<"RESULTS="<<_HOME_INFO<<"/marvin/marvin/astar/"<<heuristic_good<<"/" + PROB_GOOD  +  "/"<<dominio_global<<"/resultado\n\n";
	outfile<<"cd ${DIR}\n\n";
        outfile<<"python3 ${FD_ROOT}/src/translate/translate.py ${FD_ROOT}/benchmarks/"<<dominio_global<<"/"<<domain_pddl_global<<" ${FD_ROOT}/benchmarks/"<<dominio_global<<"/"<<tarefa_global<<"\n\n";

        outfile<<"${FD_ROOT}/src/preprocess/preprocess < output.sas"<<"\n\n"; 
	//Santiago's code to find the F_boundary on the fly


	if (apply_max) {
		outfile<<"${FD_ROOT}/src/search/downward-release --use_saved_pdbs --domain_name "<<dominio_global<<" --problem_name "<<tarefa_global<<" --heuristic_name "<<heuristic_good<<" --problem_name_gapdb "<<prob_name_gapdb<<" --deep_F_boundary "<<deep_F_boundary<<"  --dir_creation "<<method<<"  --search \"astar_original(max(["<<parameter<<"]))\" <  output > ${RESULTS}/"<<prob_name_gapdb<<"\n\n";
	} else {
	 	outfile<<"${FD_ROOT}/src/search/downward-release --use_saved_pdbs --domain_name "<<dominio_global<<" --problem_name "<<tarefa_global<<" --heuristic_name "<<heuristic_good<<" --problem_name_gapdb "<<prob_name_gapdb<<" --deep_F_boundary "<<deep_F_boundary<<"  --dir_creation "<<method<<"  --search \"astar_original("<<parameter<<")\" <  output > ${RESULTS}/"<<prob_name_gapdb<<"\n\n";
	}

	outfile<<"\n\nrm ${DIR}\n\n";
        outfile<<"\n\nmv sas_plan ${FD_ROOT}"<<plan_dir_file<<"\n\n";

        outfile.close();

	string executeFile;
	bool is_in_cluster = false;

        if (is_in_cluster) {
        	executeFile = "qsub -l select=1:ncpus=1:mem=6GB "+arquivo;
        	cout<<executeFile<<"\n\n";
		if(system(executeFile.c_str())) {
			cout<<"running in the cluster...\n";
		}
	} else {
        	string allow;
                allow = "chmod +x "+arquivo;
                cout<<allow<<"\n";
                if(system(allow.c_str())) {
			cout<<"adding permition...\n";
		}

		executeFile = "timeout "+ good_timer +" sh "+arquivo; //setting the limit time	
                cout<<executeFile<<"\n\n";
                if(system(executeFile.c_str())) {
			cout<<"running in the local...\n";
		}
	}
}

double SSSearch::getProbingResult() {
 //static map<boost::dynamic_bitset<>,double> prev_collector;
        double expansions = 0;
        
        for (size_t i = 0; i < vweight.size(); i++) {
             double n = vweight.at(i);
             expansions += n;
        }
	cout << endl;
        cout<<"expansions = "<<expansions<<endl;
	boost::dynamic_bitset<> max_comb(heuristics.size()+lmcut_heuristic.size());max_comb.set();
	cout <<std::scientific<<",max_comb_nodes:"<<collector[max_comb]<<endl;
	  /*if(collector[max_comb]>pow(10.0,100.0)){
	  cout<<"Sampling went too far, selecting heuristics now with data from prev next_f_boundary!"<<endl;
	  //collector=prev_collector;
	  select_best_heuristics_greedy();
	  exit(0);
	}*/
	/*  else{
	  prev_collector=collector;
	}*/
        return expansions;
}

double SSSearch::getMeanHeurResult() {
	double mean = 0;
	for (size_t i = 0; i < vmeanheur.size(); i++) {
		SSQueue n = vmeanheur.at(i);
		SSNode node = n.getNode();
		Type t = n.getT();
		int heur_value = t.getH();
		double w = node.getWeight();
		mean += w*heur_value;	
	}
	return mean;
}

void SSSearch::printQueue() {
        cout<<"\nPrintQueue\n";
	for (map<Type, SSNode>::iterator iter = queue.begin(); iter !=  queue.end(); iter++) {
            SSNode t2  = iter->second;
            cout<<"\t\t  g = "<<t2.getGreal()<<"  w = "<<t2.getWeight()<<"\n"; 
        }
        cout<<"\n";
        cout<<"\nEnd PrintQueue\n";
}

bool SSSearch::isGAPDB(string heur) {
        //insert all the heuristics in a vector
        string ga_pdb_s = heur;
        size_t ga_pdb_t = ga_pdb_s.find("_");
        string ga_pdb_new_s = ga_pdb_s.substr(0, ga_pdb_t);
        //end insert all the heuristics in a vector
        bool is_gapdb = false;
        if ("gapdb" == ga_pdb_new_s) {
                is_gapdb = true;
        }
        return is_gapdb;
}

int SSSearch::getTotalGAHeurs(vector<string> v) {
        int total_ga_heur = 0;
        //insert all the heuristics in a vector
        for (size_t i = 0; i < v.size(); i++) {
                string ga_pdb_s = v.at(i);
                size_t ga_pdb_t = ga_pdb_s.find("_");
                string ga_pdb_new_s = ga_pdb_s.substr(0, ga_pdb_t);
                //end insert all the heuristics in a vector
                if ("gapdb" == ga_pdb_new_s) {
                        total_ga_heur++;
                }
        }
        return total_ga_heur;
}

void SSSearch::select_random_greedy(bool termination) {	
	if (termination) {	
		//cout<<"call heuristicCombinator:\n";
		vector<pair<string, double> > Z_subset;
		vector<pair<string, double> > Z_full_vector_combiner;
		vector<pair<string, double> > Z_cut_vector;	
		bool call_first_time = true;
		while ((int)Z_subset.size() < N_default) {
			if (N_default < (int)all_heuristics.size()) {
				
				map<string, double> Z_full_map = heuristicCombinator(call_first_time, Z_subset, Z_full_vector_combiner);
				call_first_time = false;
				
				//order the map
				vector<pair<string, double> > Z_full_vector(Z_full_map.begin(), Z_full_map.end());
				sort(Z_full_vector.begin(), Z_full_vector.end(), greater_second<string, double>());
				cout<<"Imprimiendo el Z_full_vector\n";
				printVectorPair(Z_full_vector);
				cout<<"End Imprimiendo el Z_full_vector\n";

				int index_max = Z_full_vector.size() - N_default;
	
				cout<<"Z_full_vector.size()="<<Z_full_vector.size()<<"\n";	
				cout<<"N_default: "<<N_default<<"\n";
				cout<<"index_max=Z_full_vector.size()-N_default="<<index_max<<"\n";
				if (index_max < 0) {
					std::vector<pair<string, double> > Z_cut_vector_aux(Z_full_vector.begin(), Z_full_vector.end());
					Z_cut_vector = Z_cut_vector_aux;
				} else {
					std::vector<pair<string, double> > Z_cut_vector_aux(Z_full_vector.begin(), Z_full_vector.end() - index_max);
					Z_cut_vector = Z_cut_vector_aux;
				}

				int size_Z_cut_vector = Z_cut_vector.size();
				cout<<"printing: Z_cut:"<<size_Z_cut_vector<<"\n";
				printVectorPair(Z_cut_vector);
				cout<<"end printing: Z_cut\n";
				int random_index;
				if (size_Z_cut_vector > 0) {	
					random_index =  RanGen2->IRandom(0, size_Z_cut_vector - 1);//g_rng() * (index_max-1);
				}  else {
					//Then Z_full is empty
					cout<<"The size of the Z_cut is zero.\n";
					break;
				}
				cout<<"random_index="<<random_index<<"\n";

				pair<string, double> Z_choosed = Z_cut_vector.at(random_index);
				string s_choosed = Z_choosed.first;
				double d_choosed = Z_choosed.second;
				cout<<"\t\t\tstring elegido para remover del Z_full= "<<s_choosed<<"\n";
				cout<<"\t\t\tdouble elegido para remover del Z_full= "<<d_choosed<<"\n";
				Z_subset.push_back(pair<string, double>(s_choosed, d_choosed));
				//remove the Z_choosed from Z_full_map and update the Z_full_vector
				map<string, double>::iterator zIter = Z_full_map.find(Z_choosed.first);
				if (zIter != Z_full_map.end()) {
					//Is found, then remove it.
					Z_full_map.erase(Z_choosed.first);
				} else {
					//Is not found.
				}
				//clear the Z_full_vector
				//Z_full_vector.clear();

				Z_full_vector.erase(std::remove(Z_full_vector.begin(), Z_full_vector.end(), Z_choosed), Z_full_vector.end());
				
				//Remove from the Z_full the heuristics that when combined with the heuristics you already have in Z_subset do not yeald in any improvements to the objective function.	
				
				cout<<"Begin sum_Heur_by_Z_subset\n";
				double sumHeur_by_Z_subset = getSumSubset(Z_subset);
				cout<<"End sum_Heur_by_Z_subset\n";
				typedef std::vector<std::pair<std::string, double> > vector_type;
        			for (vector_type::const_iterator pos = Z_full_vector.begin();
        			pos != Z_full_vector.end(); ++pos)
        			{
					
					vector<pair<string, double> > Z_subset_aux(Z_subset.begin(), Z_subset.end());
					Z_subset_aux.push_back(*pos);
        				if (sumHeur_by_Z_subset == getSumSubset(Z_subset_aux)) {
						//remove pos
						cout<<"Removing pos\n";
						Z_full_vector.erase(std::remove(Z_full_vector.begin(), Z_full_vector.end(), *pos), Z_full_vector.end());
					} else {
						cout<<"not removing\n";
					}
					Z_subset_aux.clear();
        			}
	
				Z_full_vector_combiner = Z_full_vector;
				//Z_full_vector_combiner must contains at this point only the heuristics without combine.. The combination only happens in heuristicCombiner. Here only must have the following:
				//h1, h2, h3, h5 - The Z_subet contains h4 then in the heuristicCombiner happens the magic of combination h1_h4, h2_h4, h3_h4, h5_h4
			}
		}
		cout<<"OUT OF WHILE\n";
		cout<<"PRINTING THE SUBSET OF HEURISTICS\n";
		printVectorPair(Z_subset);
		cout<<"ENDING PRINTING THE SUBSET OF HEURISTICS\n";
		string STORE_PLAN = "GRHS";
		string heuristic_good = "gapdb_good";
		//create new variable called deep_F_boundary
		int deep_F_boundary = threshold;
      		//cout<<"deep_F_boundary = "<<deep_F_boundary<<"\n";
		string PROB_GOOD = "problemas_";
                PROB_GOOD += boost::lexical_cast<std::string>(n_probes_global);
                PROB_GOOD += "_probes_grhs";
                //cout<<"PROB_GOOD = "<<PROB_GOOD<<"\n";
        	//begin
		vector<string> v_gapdb_string;//store the heuristics with properties

		string method = "grhs";	
		mkdirAstar(method, heuristic_good, PROB_GOOD);
	
		string dirSASPLAN = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan_"+STORE_PLAN+"/";
                if (system(dirSASPLAN.c_str())) {
                	cout<<"create directory "<<dirSASPLAN.c_str()<<"\n";
                }

		string dirSASPLANDomain = "mkdir "+_HOME_INFO+"/marvin"+_FD_INFO+"/plan_"+STORE_PLAN+"/"+dominio_global;
                if (system(dirSASPLANDomain.c_str())) {
                	cout<<"create directory "<<dirSASPLANDomain.c_str()<<"\n";
                }	

        	for (size_t i = 0; i < Z_subset.size(); i++) {
			pair<string, double> p = Z_subset.at(i);
			string s_subset = p.first;

			map<string, vector<string> >::iterator iter;
			for (iter = map_heur_properties.begin(); iter != map_heur_properties.end(); iter++) {
				string gapdb_string; //the name of each heuristic, just remember that the fd only support gapdb and do not gapdb_deep or gapdb_good => both need to be changed to gapdb
                		string s_prop = iter->first;
                		vector<string> info = iter->second;
				if (s_subset == s_prop) {
					gapdb_string = processHeuristicProperties(s_subset, info, heuristic_good);
                                	v_gapdb_string.push_back(gapdb_string);
					cout<<"s_subset="<<s_subset<<"\n";	
				}
			}
        	}

		string heuristic_generator;
		for (size_t i = 0; i < v_gapdb_string.size(); i++) {
			//get the real name
                	string real_heur = v_gapdb_string.at(i);
                	string task = real_heur;
                	//cout<<"task = "<<task<<"\n";
                	//size_t found_task_deep = task.find("deep");
                	size_t found_task_good = task.find("good");
                	string final_real_heur, final_number_heur;
                	string delimiter = "_";
                	if (found_task_good > 1000) {
				string t0 = real_heur;
                        	size_t found_t0 = t0.find("_");
                        	string previous_real_heur = t0.substr(0, found_t0);			
				final_real_heur = previous_real_heur;
                        	//cout<<"previous_real_heur = "<<final_real_heur<<"\n";
                        	if (previous_real_heur == "mands()") {
                                	final_real_heur = "merge_and_shrink(shrink_strategy=shrink_bisimulation(max_states=100000,threshold=1,greedy=false),merge_strategy=merge_dfp())";
                                	//final_real_heur = "merge_and_shrink()";
                        	} else if (previous_real_heur == "ipdb()") {
                                	final_real_heur = "ipdb(max_time=200)";
                        	}

                        	//get the heuristic number
                        	string t1 = real_heur;
                        	size_t found_t1 = t1.find("_");
                        	final_number_heur = t1.substr(found_t1 + 1, t1.length());
                	} else {
                        	string s2 = real_heur;
                        	string pot[6];
                        	size_t pos = 0;
                        	string token;
                        	int index = 0;
                        	while ((pos = s2.find(delimiter)) != std::string::npos) {
                                	token = s2.substr(0, pos);
                                	pot[index] = token;
                                	s2.erase(0, pos + delimiter.length());
                                	index++;
                        	}
                        	//cout<<"index = "<<index<<"\n";
                        	pot[index] = s2;
                        	//remove deep from pot[1]
                        	string pot1 = pot[1];
                        	size_t found_pot1 = pot1.find("(");
                        	string new_pot1 = pot1.substr(found_pot1, pot1.length());
                        	//end remove deep from pot[1]

                        	final_real_heur = "gapdb" + new_pot1;
                        	final_number_heur = pot[2];
			}

			if (i != v_gapdb_string.size() - 1) {
				heuristic_generator+=final_real_heur+",";
			} else {
				heuristic_generator+=final_real_heur;
			}
			cout<<"final_real_heur = "<<final_real_heur<<"\n";
                	cout<<"final_number_heur = "<<final_number_heur<<"\n\n";
		}

		cout<<"heuristic_generator="<<heuristic_generator<<"\n";
		
                string new_problem_name = tarefa_global;//--problem_name == problema.c_str()
		string t = new_problem_name;
        	size_t found = t.find(".");
        	string new_problem_name_mod = t.substr(0, found);
        	string prob_name_gapdb = new_problem_name_mod + "_gapdb_grhs.pddl";
        	string arquivo;
        	arquivo = new_problem_name_mod + "_gapdb_grhs.sh";
		arquivo = "/" + arquivo;
        	arquivo = dominio_global + arquivo;
        	arquivo = "astar/"+heuristic_good+"/" + PROB_GOOD  +  "/" + arquivo;
        	arquivo = "marvin/" + arquivo;
        	arquivo = "marvin/"+ arquivo;
        	arquivo =  _HOME_INFO+"/" + arquivo;
        	ofstream outfile(arquivo.c_str(), ios::out);

		string plan_dir_file = "/plan_"+STORE_PLAN+"/"+dominio_global+"/"+tarefa_global;
		
		executeQsub(arquivo, heuristic_generator, heuristic_good, PROB_GOOD, prob_name_gapdb, deep_F_boundary, method, plan_dir_file, true);	
	}//end termination
}

void SSSearch::updateGRHS() {
        size_t found = tarefa_global.find(".");
        string name = tarefa_global.substr(0, found);
        name+="_F_";
        name+=boost::lexical_cast<std::string>(threshold);
        name += ".csv";
	
	string dirDomain_greedy, dirDomain, dirSSCC, outputFile;
	string nameProbes = "reportss_";
        nameProbes += boost::lexical_cast<std::string>(n_probes_global);
       	nameProbes += "_probes_grhs";
        cout<<"nameProbes = "<<nameProbes<<"\n";

	dirDomain_greedy = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+nameProbes;
	dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global;
        dirSSCC = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/bc";
        outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica_global+"/"+ nameProbes +"/"+dominio_global+"/bc/"+name;

	if (system(dirDomain_greedy.c_str())) {
        	cout<<"Directory: "<<dirDomain_greedy.c_str()<<" created."<<endl;
        }

        if (system(dirDomain.c_str())) {
        	cout<<"Directory: "<<dirDomain.c_str()<<" created."<<endl;
        }

        if (system(dirSSCC.c_str())) {
           	cout<<"Directory: SSCC created."<<endl;
        }

	ofstream output;
        output.open(outputFile.c_str());
	
	//print the name of the all_heuristics just to be analyzed later
        for (size_t i = 0; i < all_heuristics.size(); i++) {
            	string heur_name = all_heuristics[i]->get_heur_name();
            	output<<"h(,"<<i<<"):,"<<heur_name<<"\n";
        }

	count_line_grhs = collector_heur.size();
	cout<<"n_heuristics_global="<<n_heuristics_global<<"\n";
	cout<<"count_line_grhs="<<count_line_grhs<<"\n";

	//get the combintation 1/0/1/0..../1/1
        harray_grhs = new int*[count_line_grhs];
        ccarray_grhs = new double*[count_line_grhs];
        for (int i = 0; i < count_line_grhs; i++) {
                harray_grhs[i] = new int[n_heuristics_global];
               	ccarray_grhs[i] = new double[1];
        }

        int counter_line = 0;
	for (map<std::vector<int>, double>::iterator iter = collector_heur.begin(); iter != collector_heur.end(); iter++) {
		vector<int> h_node_v = iter->first;
		output<<"bc(";
		double cc = iter->second;
		for (size_t i = 0; i < h_node_v.size(); i++) {
			int h = h_node_v.at(i);
			output<<h;
			//cout<<h;	
			harray_grhs[counter_line][i] = h_node_v.at(i);
			if (i != h_node_v.size() - 1) {
				output<<"/";
				//cout<<"/";
			}
		}
                double result_cc = (double)cc/(double)n_probes_global;
                output<<")cc="<<result_cc<<"\n";
                ccarray_grhs[counter_line][0] = result_cc;
                counter_line++;
	}
	output.close();
}

//double SSSearch::calculate_time_costs_specific(boost::dynamic_bitset<> h_comb, vector<Heuristic *> heuristics){
double SSSearch::calculate_time_costs_specific(boost::dynamic_bitset<> h_comb){
  //double time_cost=1.1/pow(10,6);//node creation and expanding average cost
  double time_cost=node_gen_and_exp_cost;//node creation and expanding meassured avg cost for current problem
 // cout<<"\t node_gen_and_exp_cost:"<<time_cost<<",";
  for (size_t i = 0; i < h_comb.size(); i++){
    if(h_comb.test(i)){//add heuristic time cost to total if heuristic is active
      time_cost+=(heuristics[i]->get_measured_TPN()/gen_to_eval_ratio);
      //cout<<",h["<<i<<"]:"<<heuristics[i]->get_measured_TPN();
    }
  }
  //cout<<",TPN:"<<time_cost<<endl;
  return time_cost;
}

map<string, double> SSSearch::heuristicCombinator(bool call_first_time, vector<pair<string, double> > Z_subset, vector<pair<string, double> > Z_full) {
	cout<<"combine heuristics..."<<call_first_time<<"\n";	
	map<string, double> Z_full_map;
	if (call_first_time) {
		updateGRHS();	
		//make it work in 30 minutes
		string delimiter = ",";
		//cout<<"heuristic-information\n";
		for (int i = 0; i < n_heuristics_global; i++) {
			vector<string> collector;
                	string s =  all_heuristics[i]->get_heur_name();
                	string pot[6];
                	size_t pos = 0;
                	string token;
                	int index = 0;
                	while ((pos = s.find(delimiter)) != std::string::npos) {
                		token = s.substr(0, pos);
				//cout<<"token="<<token<<"\n";
                        	pot[index] = token;
                        	s.erase(0, pos + delimiter.length());
                        	index++;
                	}
                	pot[index] = s;

           		string heuristic_name_created = pot[0],
	    		number_h = std::to_string(i), //consider this order because SS commands	
			mutation_rate,
			mutation_rate_aux,
			size_gapdb,
			size_gapdb_aux,
			wd,
			wd_aux,
			name;


			if (false) {
	    			if (heuristic_name_created == "ipdb") {
            				name = number_h + "_ipdb";
            			} else if (heuristic_name_created == "lmcut") {
            				name = number_h + "_lmcut";
            			} else if (heuristic_name_created == "merge_and_shrink") {
            				name = number_h + "_mands";
            			} else {
            				name = number_h + "_gapdb";
            			}
			} else {
				if (heuristic_name_created == "ipdb") {
            				name = "ipdb_" + number_h;
            			} else if (heuristic_name_created == "lmcut") {
            				name = "lmcut_" + number_h;
            			} else if (heuristic_name_created == "merge_and_shrink") {
            				name = "mands_" + number_h;
            			} else {
            				name = "gapdb_" + number_h;
            			}
			}

			//cout<<"name="<<name<<"\n";
			mutation_rate_aux = pot[1];
                	size_t t2 = mutation_rate_aux.find(":");
                	mutation_rate = mutation_rate_aux.substr(t2 + 1, mutation_rate_aux.length());
                	//cout<<"mutation_rate = "<<mutation_rate<<"\n";

			size_gapdb_aux = pot[2];
                	size_t t3 = size_gapdb_aux.find("=");
                	size_gapdb = size_gapdb_aux.substr(t3 + 1);
                	//cout<<"size_gapdb = "<<size_gapdb<<"\n";

                	//without disjoint patterns
                	wd_aux = pot[3];
                	//cout<<"wd_aux = "<<wd_aux<<"\n";
                	size_t t4 = wd_aux.find("out");
                	if (t4 < 100) {
                        	wd = "false";
                	} else {
                        	wd = "true";
                	}
                	//cout<<"wd = "<<wd<<"\n\n";

			collector.push_back(number_h);
                	collector.push_back(mutation_rate);
                	collector.push_back(size_gapdb);
                	collector.push_back(wd);

			map_heur_properties.insert(pair<string, vector<string> >(name, collector));
			//number of heuristics values in the search tree.
			double sum_heur_values = 0;
			for (int j = 0; j < count_line_grhs; j++) {
				sum_heur_values += harray_grhs[j][i]*ccarray_grhs[j][0];
			}
                	Z_full_map.insert(pair<string, double>(name, sum_heur_values));
        	}
	}  else {//end_true call_first_time
		
		cout<<"Dentro del heuristicCombinator Imprimiendo Z_subset y Z_full respectivamente\n";
		printVectorPair(Z_subset);
		printVectorPair(Z_full);
		cout<<"End Dentro del heuristicCombinator Imprimiendo Z_subset y Z_full respectivamente\n";
		bool use_sum_of_the_max = true; //use the sum of the max instead of the sum of the sum
	
		cout<<"-----------------------BEGIN PRINT EACH COMBINATION------------------------\n";
		typedef std::vector<std::pair<std::string, double> > vector_type;		
		for (vector_type::const_iterator pos = Z_full.begin();
                pos != Z_full.end(); ++pos)
                {
                	string s = pos->first;
			string new_s = s;
			size_t f_s = new_s.find("_");
			//string first_key_string = new_s.substr(0, f_s);
			string first_key_string = new_s.substr(f_s + 1, new_s.length());
			//cout<<"first_key_string="<<first_key_string<<"\n";
			int first_key = std::atoi(first_key_string.c_str());
			double sum_combination_heur = 0, sum_max_combination_heur = 0;
			string key_name = s;
			map<int, vector<int> > map_max;
			map<int, double> map_cc;
			typedef std::vector<std::pair<std::string, double> > vector_type_inner;
			for (vector_type_inner::const_iterator posinner = Z_subset.begin();
			posinner != Z_subset.end(); ++posinner) {
				string s_inner = posinner->first;
				string new_s_inner = s_inner;
				size_t f_s_inner = new_s_inner.find("_");
				//string second_key_string = new_s_inner.substr(0, f_s_inner);
				string second_key_string = new_s_inner.substr(f_s_inner + 1, new_s_inner.length());
				//cout<<"second_key_string="<<second_key_string<<"\n";
				int second_key = std::atoi(second_key_string.c_str());
				string delimiter = ",";
				//cout<<"first_key = "<<first_key<<", second_key = "<<second_key<<"\n";
				//cout<<"heuristic-information\n";
				for (int i = 0; i < n_heuristics_global; i++) {
					if (first_key == i || second_key == i) {
                				string s =  all_heuristics[i]->get_heur_name();
                				string pot[6];//six is the max string allowed for heuristic string
                				size_t pos = 0;
                				string token;
                				int index = 0;
                				while ((pos = s.find(delimiter)) != std::string::npos) {
                					token = s.substr(0, pos);
							//cout<<"token="<<token<<"\n";
                        				pot[index] = token;
                        				s.erase(0, pos + delimiter.length());
                        				index++;
                				}
                				pot[index] = s;

           					string heuristic_name_created = pot[0],
	    					number_h = std::to_string(i), //consider this order because SS commands	
						name;

						if (false) {
	    						if (heuristic_name_created == "ipdb") {
            							name = number_h + "_ipdb";
            						} else if (heuristic_name_created == "lmcut") {
            							name = number_h + "_lmcut";
            						} else if (heuristic_name_created == "merge_and_shrink") {
            							name = number_h + "_mands";
            						} else {
            							name = number_h + "_gapdb";
            						}
						} else {
							if (heuristic_name_created == "ipdb") {
            							name = "ipdb_" + number_h;
            						} else if (heuristic_name_created == "lmcut") {
            							name = "lmcut_" + number_h;
            						} else if (heuristic_name_created == "merge_and_shrink") {
            							name = "mands_" + number_h;
            						} else {
            							name = "gapdb_" + number_h;
            						}
						}
						//cout<<"name="<<name<<"\n";
						//number of heuristics values in the search tree.
						double sum_heur_values = 0;
						if (use_sum_of_the_max) {
							for (int j = 0; j < count_line_grhs; j++) {
								map<int, vector<int> >::iterator itmap = map_max.find(j);
								if (itmap != map_max.end()) {
									//int bring_c = itmap->first;
									vector<int> bring_v = itmap->second;
									bring_v.push_back(harray_grhs[j][i]);
									itmap->second = bring_v;
									//cout<<"count_line_repeated\n";
								} else {
									//cout<<"count_line="<<j<<"\n";
									vector<int> add_max_heuristic_values;
									add_max_heuristic_values.push_back(harray_grhs[j][i]);
									map_max.insert(pair<int, vector<int> >(j, add_max_heuristic_values));	
								}

								map<int, double>::iterator itmap2 = map_cc.find(j);
								if (itmap2 != map_cc.end()) {
									//DO NOTHING
								} else {
									map_cc.insert(pair<int, double>(j, ccarray_grhs[j][0]));
								}
							}
						} else { 
							for (int j = 0; j < count_line_grhs; j++) {
								sum_heur_values += harray_grhs[j][i]*ccarray_grhs[j][0];
							}
							cout<<"i="<<i<<", sum_heur_values="<<sum_heur_values<<"\n";
							sum_combination_heur += sum_heur_values;
						}
						if (first_key == i) {
							first_key = -1;//set -1 in order to avoid enter again
						}
					} // key comparator
        			}//end n_heuristics loop
				//cout<<"count_line="<<count_line<<"\n";
				//cout<<"map_max.size()="<<map_max.size()<<"\n";
				if (use_sum_of_the_max) {
					for (map<int, vector<int> >::iterator it = map_max.begin(); it != map_max.end(); it++) {
						int row = it->first;
						vector<int> v = it->second;
						int max_heuristic = getMaxHeur(v);
						//cout<<"max_heuristic="<<max_heuristic<<"\n";
						map<int, double>::iterator it_inner = map_cc.find(row);
						double cc = it_inner->second;
						//cout<<"cc="<<cc<<"\n";
						sum_max_combination_heur += max_heuristic*cc;
						//cout<<"max_heuristic="<<max_heuristic<<", cc="<<cc<<"\n";
					}	
					//cout<<"sum_max_combination_heur="<<sum_max_combination_heur<<"\n";
					map_max.clear();
					map_cc.clear();
				}
			}//end Z_subset

			if (use_sum_of_the_max) {
				cout<<"key_name ending="<<key_name<<", sum_max="<<sum_max_combination_heur<<"\n";
				Z_full_map.insert(pair<string, double>(key_name, sum_max_combination_heur));
			} else {
				cout<<"key_name ending="<<key_name<<", sum_combination="<<sum_combination_heur<<"\n";
				Z_full_map.insert(pair<string, double>(key_name, sum_combination_heur));
			}
		}//end Z_full
		cout<<"-------------------------END PRINT EACH COMBINATION------------------------\n";	
	}	
	cout<<"Z_full_map empty\n";
	return Z_full_map;
}

double SSSearch::getSumSubset(vector<pair<string, double> > Z_subset) {
	bool use_sum_of_the_max = true;
	double sum_combination_heur = 0, sum_max_combination_heur = 0;
	map<int, vector<int> > map_max;
	map<int, double> map_cc;
	typedef std::vector<std::pair<std::string, double> > vector_type;
	for (vector_type::const_iterator posinner = Z_subset.begin();
	posinner != Z_subset.end(); ++posinner) {
		string s_inner = posinner->first;
		string new_s_inner = s_inner;
		size_t f_s_inner = new_s_inner.find("_");
		string key_string = new_s_inner.substr(f_s_inner + 1, new_s_inner.length());
		int key = std::atoi(key_string.c_str());
		string delimiter = ",";
		cout<<"key = "<<key<<"\n";
		//cout<<"heuristic-information\n";
		for (int i = 0; i < n_heuristics_global; i++) {
			if (key == i) {
                		string s =  all_heuristics[i]->get_heur_name();
                		string pot[6];//six is the max string allowed for heuristic string
                		size_t pos = 0;
                		string token;
                		int index = 0;
                		while ((pos = s.find(delimiter)) != std::string::npos) {
                			token = s.substr(0, pos);
					//cout<<"token="<<token<<"\n";
                        		pot[index] = token;
                        		s.erase(0, pos + delimiter.length());
                        		index++;
                		}
                		pot[index] = s;

           			string heuristic_name_created = pot[0],
	    			number_h = std::to_string(i), //consider this order because SS commands	
				name;

				if (false) {
	    				if (heuristic_name_created == "ipdb") {
            					name = number_h + "_ipdb";
            				} else if (heuristic_name_created == "lmcut") {
            					name = number_h + "_lmcut";
            				} else if (heuristic_name_created == "merge_and_shrink") {
            					name = number_h + "_mands";
            				} else {
            					name = number_h + "_gapdb";
            				}
				} else {
					if (heuristic_name_created == "ipdb") {
            					name = "ipdb_" + number_h;
            				} else if (heuristic_name_created == "lmcut") {
            					name = "lmcut_" + number_h;
            				} else if (heuristic_name_created == "merge_and_shrink") {
            					name = "mands_" + number_h;
            				} else {
            					name = "gapdb_" + number_h;
            				}
				}

				//cout<<"name="<<name<<"\n";
				//number of heuristics values in the search tree.
				double sum_heur_values = 0;
				if (use_sum_of_the_max) {
					for (int j = 0; j < count_line_grhs; j++) {
						map<int, vector<int> >::iterator itmap = map_max.find(j);
						if (itmap != map_max.end()) {
							//int bring_c = itmap->first;
							vector<int> bring_v = itmap->second;
							bring_v.push_back(harray_grhs[j][i]);
							itmap->second = bring_v;
							//cout<<"count_line_repeated\n";
						} else {
							//cout<<"count_line="<<j<<"\n";
							vector<int> add_max_heuristic_values;
							add_max_heuristic_values.push_back(harray_grhs[j][i]);
							map_max.insert(pair<int, vector<int> >(j, add_max_heuristic_values));	
						}

						map<int, double>::iterator itmap2 = map_cc.find(j);
						if (itmap2 != map_cc.end()) {
							//DO NOTHING
						} else {
							map_cc.insert(pair<int, double>(j, ccarray_grhs[j][0]));
						}
					}
				} else { 
					for (int j = 0; j < count_line_grhs; j++) {
						sum_heur_values += harray_grhs[j][i]*ccarray_grhs[j][0];
					}
					cout<<"i="<<i<<", sum_heur_values="<<sum_heur_values<<"\n";
					sum_combination_heur += sum_heur_values;
				}	
			} // key comparator
        	}//end n_heuristics loop
		//cout<<"count_line="<<count_line<<"\n";
		//cout<<"map_max.size()="<<map_max.size()<<"\n";
		if (use_sum_of_the_max) {
			for (map<int, vector<int> >::iterator it = map_max.begin(); it != map_max.end(); it++) {
				int row = it->first;
				vector<int> v = it->second;
				int max_heuristic = getMaxHeur(v);
				//cout<<"max_heuristic="<<max_heuristic<<"\n";
				map<int, double>::iterator it_inner = map_cc.find(row);
				double cc = it_inner->second;
				//cout<<"cc="<<cc<<"\n";
				sum_max_combination_heur += max_heuristic*cc;
				//cout<<"max_heuristic="<<max_heuristic<<", cc="<<cc<<"\n";
			}	
			//cout<<"sum_max_combination_heur="<<sum_max_combination_heur<<"\n";
		}
	}//end Z_subset

	if (use_sum_of_the_max) {
		cout<<"returning the max.\n";
		return sum_max_combination_heur;
	} else {
		cout<<"returning no max.\n";
		return sum_combination_heur;
	}
}

void SSSearch::printVectorPair(vector<pair<string, double> > vpair) {	
	cout<<"\n\t-------------beging print vector pair---------\n";
	typedef std::vector<std::pair<std::string, double> > vector_type;
        for (vector_type::const_iterator pos = vpair.begin();
        pos != vpair.end(); ++pos)
        {
        	string s = pos->first;
                double d = pos->second;
                cout<<"\t\t\t"<<s<<"  -  "<<d<<"\n";
        }
	cout<<"\t-------------end print vector pair-----------\n";
}

void SSSearch::select_best_heuristics_greedy(){
  //h_comb_current_time.clear();
  vector<pair<int,double> > selected_heuristics;
  cout<<"Starting select_best_h_combs_RIDA_Sampling_Greedy,bcs:"<<collector.size()<<endl;
  /*  double probe_percentage=(double)(last_probe+1)/(double)ss_probes;
  if(full_collector.size()>0&&probe_percentage<=0.75){
    cout<<"using last fully expanded collector"<<endl;
    collector=full_collector;
  }
  cout<<"probe_percentage:"<<probe_percentage<<endl;*/
  
  //double best_current_time=0;
  //long best_initial_nodes=0;
  boost::dynamic_bitset<> current_h_comb(heuristics.size()+lmcut_heuristic.size());
  boost::dynamic_bitset<> best_h_comb(heuristics.size()+lmcut_heuristic.size());
  //double last_best_h_comb_time=0;
  
  //First Greedy selection
  Timer progressive_greedy_timer;
  best_h_comb.reset();


  //Timer calculate_heur_timer;
  vector<int> active_heurs;
  vector<int> candidate_heurs;
  static set<int> skip_heurs;
  bool delete_counter;
    
  cout<<"finding greedily best h_comb, starting degree:"<<best_h_comb.count()<<endl;
  double best_comb_nodes=DBL_MAX;
  cout<<"starting best_comb_nodes="<<best_comb_nodes<<endl;
  cout<<"lmcut_heuristic.size():"<<lmcut_heuristic.size()<<",heuristics.size:"<<heuristics.size()<<",collector bc size:"<<(collector.begin())->first.size()<<endl;
  
	  

  vector<double> initial_heur_sizes(heuristics.size()+lmcut_heuristic.size());
  if(best_h_comb.count()==0){
    for (map<boost::dynamic_bitset<>,double>::const_iterator it_h_values= collector.begin(); it_h_values != collector.end();it_h_values++){
	//cout<<"it_h_values:"<<it_h_values->first<<endl;
      for(size_t heur=0;heur<heuristics.size()+lmcut_heuristic.size();heur++){
	if(it_h_values->first.test(heur)){//if indiv heuristic is active and culling the counter then this counter does not apply to the h_comb
	  initial_heur_sizes[heur]+=it_h_values->second;
	  //cout<<"added "<<it_h_values->second<<" to initial_heur_"<<heur<<endl;
	}
	else{
	  //cout<<"skipped added "<<it_h_values->second<<" to initial_heur_sizes of "<<heur<<endl;
	}
      }
    }
  //Now calculate which one was best
    std::cout << "scientific:\n" << std::scientific;
    for(size_t heur=0;heur<initial_heur_sizes.size();heur++){
      cout<<"initial_heur_sizes("<<heur<<"):"<<initial_heur_sizes.at(heur)<<endl;
#ifdef _LMCUT_EARLY_TERM
	  if((heur==heuristics.size())&&lmcut_heuristic.size()>0){//so we are dealing with lmcut
	    cout<<"skipping lmcut for now, will check later once best ipdb/gapdbs comb is found"<<endl;
	    continue;
	  }
#endif
      if(initial_heur_sizes.at(heur)<best_comb_nodes){
	best_comb_nodes=initial_heur_sizes.at(heur);
	best_h_comb.resize(initial_heur_sizes.size());best_h_comb.reset();best_h_comb.set(heur);
	cout<<"temp best_comb_nodes:"<<best_comb_nodes<<",temp best heur:";print_h_comb(best_h_comb);cout<<endl;
      }
    }
  }
    
  for(unsigned i=0;i<best_h_comb.size();i++){
    if(best_h_comb.test(i)){
      selected_heuristics.push_back(make_pair(i,best_comb_nodes));
    }
  }
  

  while(true){
    cout<<"prev best_h_comb:";print_h_comb(best_h_comb);cout<<",nodes:"<<best_comb_nodes<<endl;

    active_heurs.clear();
    candidate_heurs.clear();
    for(unsigned i=0;i<best_h_comb.size();i++){
      if(best_h_comb.test(i)){
	active_heurs.push_back(i);
      }
      else if(skip_heurs.find(i)==skip_heurs.end()){
	candidate_heurs.push_back(i);
      }
    }
    cout<<"active_heurs"<<active_heurs<<",remaining candidates:"<<candidate_heurs.size()<<endl;
    //cout<<"candidate_heurs:"<<candidate_heurs<<endl;

    vector<double> h_comb_prunning_nodes(heuristics.size()+lmcut_heuristic.size(),0);
    for (map<boost::dynamic_bitset<>,double>::iterator it_h_values= collector.begin(); it_h_values != collector.end();){
      //cout<<"collector.size:"<<it_h_values->first.size()<<",bc:"<<it_h_values->first<<endl;
      delete_counter=false;
      for(size_t heur=0;heur<active_heurs.size();heur++){
	//cout<<"testing :"<<active_heurs[heur]<<endl;
	if(!it_h_values->first.test(active_heurs[heur])){//if indiv heuristic is active and culling the counter then this counter does not apply to the h_comb
	  delete_counter=true;
	  //cout<<"Erasing "<<it_h_values->first<<endl;
	  collector.erase(it_h_values++);//removing CC, this CC is already culled by previous best candidate set
	  break;
	}
      }
      if(delete_counter){
	continue;
      }
      //So this counter is not being culled by any of the preselected heuristics, lets keep track of the heuristics which would prune it


      for(size_t heur=0;heur<candidate_heurs.size();heur++){
	if(!it_h_values->first.test(candidate_heurs[heur])){//if indiv heuristic is active and culling the counter then this counter shows heur actually prunning extra nodes
	  //cout<<"\t working on counter:"<<it_h_values->first<<",prunning heur:"<<heur<<",adding:"<<it_h_values->second<<endl;
	  h_comb_prunning_nodes[candidate_heurs[heur]]+=it_h_values->second;
	}
      }
      ++it_h_values;
    }
    //cout<<"Finished checking collectors"<<endl;fflush(stdout);
    //cout<<"h_comb_prunning_nodes.size:"<<h_comb_prunning_nodes.size()<<endl;
      //now check which heur is best to add to h_comb and wether any of the remaining h_combs can reduce the size of the prev best combinaiton
      double highest_prunner=0;
      int best_candidate=0;
      for(size_t i=0;i<candidate_heurs.size();i++){
	if(h_comb_prunning_nodes.at(candidate_heurs[i])==0){
	  //cout<<"skipping h("<<i<<"), it can not add further prunning"<<endl;fflush(stdout);
	  skip_heurs.insert(candidate_heurs[i]);
	}
	else if(h_comb_prunning_nodes.at(candidate_heurs[i])>highest_prunner){
	  highest_prunner=h_comb_prunning_nodes.at(candidate_heurs[i]);
	  best_candidate=candidate_heurs[i];
	}
      }
      if(highest_prunner==0){
	cout<<"no possible improvements node-wise"<<endl;
	  break;
      }
      else{
	best_h_comb.set(best_candidate);
	best_comb_nodes-=highest_prunner;
	selected_heuristics.push_back(make_pair(best_candidate,best_comb_nodes));
	cout<<"best_h_comb:";print_h_comb(best_h_comb);cout<<",nodes:"<<best_comb_nodes<<",extra pruned:"<<highest_prunner<<endl;
	
      }
  }
  if(best_h_comb.count()==0){
    cout<<"Best h comb is empty, setting first heursitic as default, cabnt send empty heuristic"<<endl;
    best_h_comb.set(0);
    selected_heuristics.push_back(make_pair(0,best_comb_nodes));
  }
  /*if(prev_F_boundaries.size()>1){
    double F_step=prev_F_boundaries.back()-prev_F_boundaries.at(prev_F_boundaries.size()-2);
    double current_F_boundary;
    double best_h_comb_HBF=double(best_comb_nodes)/double(prev_nodes);
    cout<<"prev_nodes:"<<prev_nodes<<"F_step:"<<F_step<<",current_F_boundary:"<<prev_F_boundaries.back()<<",best_h_comb HBF:"<<best_h_comb_HBF;
      
    //double next_nodes=best_comb_nodes;
    ofstream Predictions_file;
    string Predictions_filename="Preds_";
    Predictions_filename+=g_plan_filename;
    Predictions_filename+=".txt";
    Predictions_file.open (Predictions_filename.c_str());
    for(int i=0;i<100;i++){
      current_F_boundary=prev_F_boundaries.back()+i*F_step;
      cout<<"F_boundary:,"<<current_F_boundary<<",size:,"<<best_comb_nodes*pow(best_h_comb_HBF,i)<<endl;
      Predictions_file<<std::fixed;
      Predictions_file<<std::setprecision(0);
      Predictions_file<<"Predictions,F_boundary:,"<<current_F_boundary<<",size:,"<<best_comb_nodes*pow(best_h_comb_HBF,i)<<endl;
    }
    Predictions_file.close();
  }
  else{
    cout<<"only one common F_boundary, so no HBF predictions for future F_boundaries"<<endl;
  }*/

  cout<<"progressive_greedy_timer:,"<<progressive_greedy_timer<<",level:,"<<best_h_comb.count()<<",Final Overall best h_comb node-wise(no lmcut considered yet):";print_h_comb(best_h_comb);cout<<"line-by-line selected heurs:"<<endl;
  for(size_t j=0;j<best_h_comb.size();j++){
    if(best_h_comb.test(j)){
      if(j<heuristics.size()){
	cout<<heuristics[j]->get_heur_call_name();
      }
      else{
	cout<<lmcut_heuristic[0]->get_heur_call_name()<<endl;
      }
    }
  }

  ofstream output;
  string output_file="temp_scripts/selected_heurs_";
  output_file+=tarefa_global;
  output_file+=".sh";
  output.open(output_file.c_str());
  output<<"./downward-release --use_saved_pdbs --domain_name "<<domain_name<<" --problem_name "<<tarefa_global<<" --search \"astar(";
  size_t populated_counter=0;
  output<<"min([";
  for(size_t j=0;j<best_h_comb.size();j++){
    if(best_h_comb.test(j)){
      if(j<heuristics.size()){
	output<<heuristics[j]->get_heur_call_name();
      }
      else{
	output<<lmcut_heuristic[0]->get_heur_call_name();
      }
      if(populated_counter<(best_h_comb.count()-1)){
	output<<",";
	//cout<<"j:"<<j<<"best_h_comb.count:"<<best_h_comb.count()<<", printing comma"<<endl;
      }
      else{
	//cout<<"j:"<<j<<",best_h_comb.size:"<<best_h_comb.count()<<", not printing comma"<<endl;
      }
      populated_counter++;
    }
  }
  output<<"]))\" ";
  output.close();

  /*
  best_current_time=double(best_current_nodes)*calculate_time_costs_specific(best_h_comb,selectable_heuristics);
  cout<<",nodes:,"<<best_current_nodes<<",best_current_time:,"<<best_current_time<<endl;*/

 cout<<"Memory before deleting states:";
 get_peak_memory_in_kb(true);
 delete g_state_registry;
 g_state_registry = new StateRegistry;
 cout<<"Memory after deleting states:";
 get_peak_memory_in_kb(true);
 exit(0);
}
void SSSearch::select_best_heuristics_greedy_stocastic(int cardinality){
  //h_comb_current_time.clear();
  vector<pair<int,double> > selected_heuristics;
  cout<<"Starting select_best_h_combs_RIDA_Sampling_Greedy,bcs:"<<collector.size()<<endl;
  /*  double probe_percentage=(double)(last_probe+1)/(double)ss_probes;
  if(full_collector.size()>0&&probe_percentage<=0.75){
    cout<<"using last fully expanded collector"<<endl;
    collector=full_collector;
  }
  cout<<"probe_percentage:"<<probe_percentage<<endl;*/
  
  //double best_current_time=0;
  //long best_initial_nodes=0;
  boost::dynamic_bitset<> current_h_comb(heuristics.size()+lmcut_heuristic.size());
  boost::dynamic_bitset<> best_h_comb(heuristics.size()+lmcut_heuristic.size());
  //double last_best_h_comb_time=0;
  
  //First Greedy selection
  Timer progressive_greedy_timer;
  best_h_comb.reset();


  //Timer calculate_heur_timer;
  vector<int> active_heurs;
  vector<int> candidate_heurs;
  static set<int> skip_heurs;
  bool delete_counter;
    
  cout<<"finding greedily best h_comb, starting degree:"<<best_h_comb.count()<<endl;
  double best_comb_nodes=DBL_MAX;
  cout<<"starting best_comb_nodes="<<best_comb_nodes<<endl;
  cout<<"lmcut_heuristic.size():"<<lmcut_heuristic.size()<<",heuristics.size:"<<heuristics.size()<<",collector bc size:"<<(collector.begin())->first.size()<<endl;
  
	  

  vector<double> initial_heur_sizes(heuristics.size()+lmcut_heuristic.size());
  vector<double> initial_heur_times(heuristics.size()+lmcut_heuristic.size());
  if(best_h_comb.count()==0){
    for (map<boost::dynamic_bitset<>,double>::const_iterator it_h_values= collector.begin(); it_h_values != collector.end();it_h_values++){
	//cout<<"it_h_values:"<<it_h_values->first<<endl;
      for(size_t heur=0;heur<heuristics.size()+lmcut_heuristic.size();heur++){
	if(it_h_values->first.test(heur)){//if indiv heuristic is active and culling the counter then this counter does not apply to the h_comb
	  initial_heur_sizes[heur]+=it_h_values->second;
	  boost::dynamic_bitset<> current_h_comb(heuristics.size());current_h_comb.set(heur);
	  //initial_heur_times[heur]=double(initial_heur_sizes.at(heur))*heuristics[heur]->get_measured_TPN();
	 // cout<<"initial_heur_times("<<heur<<"):"<<initial_heur_times.at(heur)<<endl;
	  //cout<<"added "<<it_h_values->second<<" to initial_heur_"<<heur<<endl;
	}
	else{
	  //cout<<"skipped added "<<it_h_values->second<<" to initial_heur_sizes of "<<heur<<endl;
	}
      }
    }
  //Now calculate which one was best
    std::cout << "scientific:\n" << std::scientific;
    for(size_t heur=0;heur<initial_heur_sizes.size();heur++){
      cout<<"initial_heur_sizes("<<heur<<"):"<<initial_heur_sizes.at(heur)<<endl;
#ifdef _LMCUT_EARLY_TERM
	  if((heur==heuristics.size())&&lmcut_heuristic.size()>0){//so we are dealing with lmcut
	    cout<<"skipping lmcut for now, will check later once best ipdb/gapdbs comb is found"<<endl;
	    continue;
	  }
#endif
      if(initial_heur_sizes.at(heur)<best_comb_nodes){
	best_comb_nodes=initial_heur_sizes.at(heur);
	best_h_comb.resize(initial_heur_sizes.size());best_h_comb.reset();best_h_comb.set(heur);
	cout<<"temp best_comb_nodes:"<<best_comb_nodes<<",temp best heur:";print_h_comb(best_h_comb);cout<<endl;
      }
    }
  }
    
  for(unsigned i=0;i<best_h_comb.size();i++){
    if(best_h_comb.test(i)){
      selected_heuristics.push_back(make_pair(i,best_comb_nodes));
    }
  }
  
  for(int comb_degree=0;comb_degree<cardinality;comb_degree++){
    cout<<"prev best_h_comb:";print_h_comb(best_h_comb);cout<<",nodes:"<<best_comb_nodes<<endl;

    active_heurs.clear();
    candidate_heurs.clear();
    for(unsigned i=0;i<best_h_comb.size();i++){
      if(best_h_comb.test(i)){
	active_heurs.push_back(i);
      }
      else if(skip_heurs.find(i)==skip_heurs.end()){
	candidate_heurs.push_back(i);
      }
    }
    cout<<"active_heurs"<<active_heurs<<",remaining candidates:"<<candidate_heurs.size()<<endl;
    //cout<<"candidate_heurs:"<<candidate_heurs<<endl;

    vector<double> h_comb_prunning_nodes(heuristics.size()+lmcut_heuristic.size(),0);
    for (map<boost::dynamic_bitset<>,double>::iterator it_h_values= collector.begin(); it_h_values != collector.end();){
      //cout<<"collector.size:"<<it_h_values->first.size()<<",bc:"<<it_h_values->first<<endl;
      delete_counter=false;
      for(size_t heur=0;heur<active_heurs.size();heur++){
	//cout<<"testing :"<<active_heurs[heur]<<endl;
	if(!it_h_values->first.test(active_heurs[heur])){//if indiv heuristic is active and culling the counter then this counter does not apply to the h_comb
	  delete_counter=true;
	  //cout<<"Erasing "<<it_h_values->first<<endl;
	  collector.erase(it_h_values++);//removing CC, this CC is already culled by previous best candidate set
	  break;
	}
      }
      if(delete_counter){
	continue;
      }
      //So this counter is not being culled by any of the preselected heuristics, lets keep track of the heuristics which would prune it


      for(size_t heur=0;heur<candidate_heurs.size();heur++){
	if(!it_h_values->first.test(candidate_heurs[heur])){//if indiv heuristic is active and culling the counter then this counter shows heur actually prunning extra nodes
	  //cout<<"\t working on counter:"<<it_h_values->first<<",prunning heur:"<<heur<<",adding:"<<it_h_values->second<<endl;
	  h_comb_prunning_nodes[candidate_heurs[heur]]+=it_h_values->second;
	}
      }
      ++it_h_values;
    }
    //cout<<"Finished checking collectors"<<endl;fflush(stdout);
    //cout<<"h_comb_prunning_nodes.size:"<<h_comb_prunning_nodes.size()<<endl;
      //now check which heur is best to add to h_comb and wether any of the remaining h_combs can reduce the size of the prev best combinaiton
      double highest_prunner=0;
      int best_candidate=0;
      for(size_t i=0;i<candidate_heurs.size();i++){
	if(h_comb_prunning_nodes.at(candidate_heurs[i])==0){
	  //cout<<"skipping h("<<i<<"), it can not add further prunning"<<endl;fflush(stdout);
	  skip_heurs.insert(candidate_heurs[i]);
	}
	else if(h_comb_prunning_nodes.at(candidate_heurs[i])>highest_prunner){
	  highest_prunner=h_comb_prunning_nodes.at(candidate_heurs[i]);
	  best_candidate=candidate_heurs[i];
	}
      }
      if(highest_prunner==0){
	cout<<"no possible improvements node-wise"<<endl;
	  break;
      }
      else{
	best_h_comb.set(best_candidate);
	best_comb_nodes-=highest_prunner;
	selected_heuristics.push_back(make_pair(best_candidate,best_comb_nodes));
	cout<<"best_h_comb:";print_h_comb(best_h_comb);cout<<",nodes:"<<best_comb_nodes<<",extra pruned:"<<highest_prunner<<endl;
	
      }
  }
  if(best_h_comb.count()==0){
    cout<<"Best h comb is empty, setting first heursitic as default, cabnt send empty heuristic"<<endl;
    best_h_comb.set(0);
    selected_heuristics.push_back(make_pair(0,best_comb_nodes));
  }
  /*if(prev_F_boundaries.size()>1){
    double F_step=prev_F_boundaries.back()-prev_F_boundaries.at(prev_F_boundaries.size()-2);
    double current_F_boundary;
    double best_h_comb_HBF=double(best_comb_nodes)/double(prev_nodes);
    cout<<"prev_nodes:"<<prev_nodes<<"F_step:"<<F_step<<",current_F_boundary:"<<prev_F_boundaries.back()<<",best_h_comb HBF:"<<best_h_comb_HBF;
      
    //double next_nodes=best_comb_nodes;
    ofstream Predictions_file;
    string Predictions_filename="Preds_";
    Predictions_filename+=g_plan_filename;
    Predictions_filename+=".txt";
    Predictions_file.open (Predictions_filename.c_str());
    for(int i=0;i<100;i++){
      current_F_boundary=prev_F_boundaries.back()+i*F_step;
      cout<<"F_boundary:,"<<current_F_boundary<<",size:,"<<best_comb_nodes*pow(best_h_comb_HBF,i)<<endl;
      Predictions_file<<std::fixed;
      Predictions_file<<std::setprecision(0);
      Predictions_file<<"Predictions,F_boundary:,"<<current_F_boundary<<",size:,"<<best_comb_nodes*pow(best_h_comb_HBF,i)<<endl;
    }
    Predictions_file.close();
  }
  else{
    cout<<"only one common F_boundary, so no HBF predictions for future F_boundaries"<<endl;
  }*/

  cout<<"progressive_greedy_timer:,"<<progressive_greedy_timer<<",level:,"<<best_h_comb.count()<<",Final Overall best h_comb node-wise(no lmcut considered yet):";print_h_comb(best_h_comb);cout<<"line-by-line selected heurs:"<<endl;
  for(size_t j=0;j<best_h_comb.size();j++){
    if(best_h_comb.test(j)){
      if(j<heuristics.size()){
	cout<<heuristics[j]->get_heur_call_name();
      }
      else{
	cout<<lmcut_heuristic[0]->get_heur_call_name()<<endl;
      }
    }
  } 

  ofstream output;
  string output_file="temp_scripts/selected_heurs_";
  output_file+=tarefa_global;
  output_file+=".sh";
  output.open(output_file.c_str());
  output<<"./downward-release --use_saved_pdbs --domain_name "<<domain_name<<" --problem_name "<<tarefa_global<<" --search \"astar(";
  size_t populated_counter=0;
  output<<"min([";
  for(size_t j=0;j<best_h_comb.size();j++){
    if(best_h_comb.test(j)){
      if(j<heuristics.size()){
	output<<heuristics[j]->get_heur_call_name();
      }
      else{
	output<<lmcut_heuristic[0]->get_heur_call_name();
      }
      if(populated_counter<(best_h_comb.count()-1)){
	output<<",";
	//cout<<"j:"<<j<<"best_h_comb.count:"<<best_h_comb.count()<<", printing comma"<<endl;
      }
      else{
	//cout<<"j:"<<j<<",best_h_comb.size:"<<best_h_comb.count()<<", not printing comma"<<endl;
      }
      populated_counter++;
    }
  }
  output<<"]))\" ";
  output.close();

  /*
  best_current_time=double(best_current_nodes)*calculate_time_costs_specific(best_h_comb,selectable_heuristics);
  cout<<",nodes:,"<<best_current_nodes<<",best_current_time:,"<<best_current_time<<endl;*/

 cout<<"Memory before deleting states:";
 get_peak_memory_in_kb(true);
 delete g_state_registry;
 g_state_registry = new StateRegistry;
 cout<<"Memory after deleting states:";
 get_peak_memory_in_kb(true);
 exit(0);
}


void SSSearch::initialize() {
	cout << "SSSearch ..." << endl;
	search_time.reset();
	level_time.reset();

	if (gen_to_eval_ratio == 0) {
        	gen_to_eval_ratio=1;
		cout<<"gen_to_eval_ratio="<<gen_to_eval_ratio<<"\n";
	}

        //ss+cc santiago code
	if (use_saved_pdbs) {
		stored_GA_patterns.clear();
		cout<<"cleared store_GA_patterns."<<endl;
	}
}
void SSSearch::BFS(SSNode root, Type type) {
	//std::queue<SSNode> D;
	std::deque<SSQueue> D;
	SSQueue s1;
	s1.setNode(root);
	s1.setT(type);
        D.push_back(s1);
	
	check.insert(root);
	
	//double weight = root.getWeight();
	//int counter = 0;
#ifdef _SS_DEBUG
                
	SSNode BFS_RootNode = D.front().getNode();
        cout<<"Starting BFS,D.size:  "<<D.size()<<",state_id="<<BFS_RootNode.get_id()<<endl;fflush(stdout);
	GlobalState global_state_2 = g_state_registry->lookup_state(BFS_RootNode.get_id());
        global_state_2.dump_pddl();
#endif
        while (!D.empty()) {
		SSQueue sq = D.front();
                SSNode nodecp = sq.getNode();
		Type t = sq.getT();

                double g_real = nodecp.getGreal();
                StateID state_id = nodecp.get_id();
		double w = nodecp.getWeight();
                double level = t.getLevel();
#ifdef _SS_DEBUG
                cout<<"\n\t\t\tBFS: Node expanded: h = "<<t.getH()<<", g_real = "<<nodecp.getGreal()<<", f = "<<t.getH() + nodecp.getGreal()<<", level = "<<level<<", w = "<<w<<", stateID,:"<<state_id<<"\n";
fflush(NULL);
#endif

                //counter++;
		D.pop_front();
#ifdef _SS_DEBUG
		cout<<"D.size:"<<D.size()<<endl;
#endif

                std::vector<const GlobalOperator *> applicable_ops;
                //Recover the global_state
                GlobalState global_state = g_state_registry->lookup_state(state_id);
                g_successor_generator->generate_applicable_ops(global_state, applicable_ops);
#ifdef _SS_DEBUG
                cout<<"\t\t\tBFS:,retrieved state_id:"<<state_id<<","<<"applicable_ops.size() = "<<applicable_ops.size()<<endl;fflush(NULL);
                cout<<"\t\t\t--------------childs-------------\n";fflush(NULL);
#endif
		for (size_t i = 0; i < applicable_ops.size(); i++) {

                        const GlobalOperator *op = applicable_ops[i];
#ifdef _SS_DEBUG
			cout<<"op["<<i<<"]:";op->dump();
#endif
                        GlobalState child =  g_state_registry->get_successor_state(global_state, *op);

                        int hmin_value = INT_MAX/2;
			int cost_op = get_adjusted_cost(*op);
			SSNode succ_node(child.get_id(), w, g_real + cost_op);

			std::pair<std::set<SSNode, classcomp2>::iterator, bool> r1;
#ifdef _SS_DEBUG
			cout<<"\t\t\tBefore check insert\n";fflush(NULL);
			cout<<"\t\t\tsucc_node id:";cout<<child.get_id()<<endl;fflush(NULL);
			GlobalState global_state_3 = g_state_registry->lookup_state(succ_node.get_id());
			global_state_3.dump_pddl();
#endif
			r1 = check.insert(succ_node);
#ifdef _SS_DEBUG
			if(r1.second){
			  cout<<"\t\t\tsucc_node id:";cout<<child.get_id()<<"is new in check"<<endl;fflush(NULL);
			}
#endif
			if (r1.second) {
				for (size_t i = 0; i < heuristics.size(); i++) {
                                	heuristics[i]->evaluate(child);
					if (!heuristics[i]->is_dead_end()) {
						hmin_value = min(hmin_value, heuristics[i]->get_heuristic());
					} else {
						hmin_value = INT_MAX/2;
						break;
					}
                        	}

                        	double succ_h = hmin_value;
				Type object = sampler->getType(child.get_id(), succ_h,  1);
                      		object.setLevel( level + 1 );
				
				SSQueue s3;
				s3.setNode(succ_node);
				s3.setT(object);

				//cout<<"\t\t\tChild_"<<(i+1)<<" : h = "<<succ_h<<", g_real = "<<succ_node.getGreal()<<", f = "<<succ_h + succ_node.getGreal()<<", level = "<<object.getLevel()<<", w = "<<w<<", stateID,:"<<child.get_id()<<"\n";
				if (succ_h + succ_node.getGreal() <= threshold) {
					if (cost_op == 0) {
#ifdef _SS_DEBUG
						cout<<"\t\t\tcost = 0, then is inserted to the D.\n";
#endif
						D.push_back(s3);			
					} else {
#ifdef _SS_DEBUG
						cout<<"\t\t\tcost != 0, then is inserted to the L.\n";
#endif
						L.insert(s3);
					}
				}
				else{
#ifdef _SS_DEBUG
				  int prev_next_f_bound=next_f_bound;
#endif
				  next_f_bound=min(next_f_bound,int(succ_h + succ_node.getGreal()));
#ifdef _SS_DEBUG
				  if(next_f_bound!=prev_next_f_bound){
				    cout<<"BFS updates next_f_bound to:"<<next_f_bound<<endl;
				  }
#endif
				}
			} else {
#ifdef _SS_DEBUG
				cout<<"\t\t\tSSQueue with id = "<<child.get_id()<<" already exists.\n";
#endif
			}
                }//end for
                //cout<<"\t\t\t-------------End childs------------\n";
	//fflush(NULL);
        }
#ifdef _SS_DEBUG
	cout<<"BFS finished, L size:"<<L.size()<<endl;
#endif
	/*double result = counter*weight;
	root.setWeight(result);
	vweight.push_back(root);
	//Insert each node.
	double new_f_value = type.getH() + root.getGreal();
        Node2 node2(new_f_value, type.getLevel());
        
        //count nodes expanded
        if (new_f_value <= threshold) {
		std::pair<std::map<Node2, double>::iterator, bool> ret0;
                std::map<Node2, double>::iterator it0;

                ret0 = expanded.insert(pair<Node2, double>(node2, weight));
                it0 = ret0.first;

                if (ret0.second) {
                    //cout<<"new node expanded is added."<<endl;
                } else {
                    //cout<<"node expanded is being updated."<<endl;
                    it0->second += weight;
                    //cout<<"it0->second = "<<it0->second<<endl;
                }
                //end count node
        }*/
	//cout<<"\t\t\tbefore memory error."<<endl;
        //cout<<"\t\t\tis std::queue empty? D.empty() == "<<D.empty()<<endl;
	//cout<<"\t\t\t return L.size() = "<<L.size()<<endl;
}

static SearchEngine *_parse(OptionParser &parser) {
        
	parser.add_option<ScalarEvaluator *>("eval");
	
	SearchEngine::add_options_to_parser(parser);
	Options opts = parser.parse();
	SSSearch *engine = 0;
	if (!parser.dry_run()) {
		engine = new SSSearch(opts);
	}
	return engine;
}


static Plugin<SearchEngine> _plugin("ss", _parse);
