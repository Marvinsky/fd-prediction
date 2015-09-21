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

//ss+culprits
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

//Root and fd information
string _HOME_INFO = "/home";
string _FD_INFO = "/fd";

SSSearch::SSSearch(const Options &opts) : SearchEngine(opts), current_state(g_initial_state()) {

	ScalarEvaluator * evaluator = opts.get<ScalarEvaluator *>("eval");
	std::set<Heuristic *> hset;
	evaluator->get_involved_heuristics(hset);
	int heuristics_active=0;
	for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); it++) {
	  //Eliminate any heuristics which were not generated because we ran out of time
	  //currently this is hacked to return not using heuristics
	  if((*it)->is_using()){
	    heuristics.push_back(*it);
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
	  else{ ////Only using 10 GA + ipdb + lmcut + merge_and_shrink, the ss_probes must b     e passed from the sh
	     cout<<"initialize:ss_probes="<<ss_probes<<"\n";    
	     //ss_probes=500;
	  }
	  cout<<"Not doing domination_check, setting probes to :"<<ss_probes<<endl;
	}
}

SSSearch::~SSSearch() {
}


SearchStatus SSSearch::step() {
  bool at_least_one_dominated=false;
  vector<int> demotted_heurs;
  int original_threshold=f_boundary;
  while(search_time()<sampling_time_limit){
	this->RanGen2 = new CRandomMersenne(1);
        	
	//clear all data from previous iteration
	queue.clear();expanded.clear();generated.clear();
	if(!domination_check){
	  collector.clear();
	  cout<<"Cleared F_culprits&b_culprits"<<endl;F_culprits.clear();
	}
        predict(ss_probes);
	if(next_f_bound==INT_MAX/2){
	  cout<<"next_f_bound was not updated!, check code!"<<endl;exit(1);
	}
	threshold=next_f_bound;
	f_boundary=next_f_bound;
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
	  else{//Running 10 GA + ipdb + lmcut + merge_and_shrink
            cout<<"into_step:"<<ss_probes<<"\n";
	    //ss_probes=500;
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
  cout<<"selecting best heuristc after,"<<search_time()<<", seconds"<<endl;
  //select_best_heuristics_greedy();
  generateSSCCReport(threshold, true);
        return SOLVED;
}

void SSSearch::predict(int probes) {
        totalPrediction = 0;
        cout<<"#probes : "<<probes<<",g_timer:"<<g_timer()<<endl;
	cout<<"input heuristics:"<<heuristics.size()<<endl;
	int current_n_probes = 0;
        for (int i = 0; i < probes; i++) {
	  cout<<"#probe:"<<i<<",g_timer:"<<g_timer()<<",search_time:"<<search_time()<<endl;
	  if(search_time()>sampling_time_limit || g_timer()>overall_time_limit){ //300 default
	    cout<<"Search_timer past maximum sampling_time"<<endl;
	    cout<<"selecting best heuristic after search_time: "<<search_time()<<", seconds,g_timer:"<<g_timer()<<endl;
	    //select_best_heuristics_greedy();	
	    generateSSCCReport(current_n_probes, true);
	    exit(0);
	  }
	  else if(search_time()>5.0&&domination_check){
	    cout<<"Search_time:"<<search_time()<<",returning for domination_check"<<endl;
	    return;
	  }
            vweight.clear();
            probe();
	    current_n_probes++;
            double p = getProbingResult();
            totalPrediction = totalPrediction + (p - totalPrediction)/(i + 1);
            cout<<"**********"<<endl;
            cout<<"p = "<<p<<endl;
            cout<<"prePre_"<<(i+1)<<" = "<<totalPrediction<<endl;
            cout<<"**********"<<endl;
        }
	cout<<"next_f_bound:"<<next_f_bound<<endl;
        ss_timer_value = ss_timer();	
        cout<<"\ntotalPrediction : "<<totalPrediction<<"\n";
        cout<<"ss_timer: "<<ss_timer_value<<"\n";
	cout<<"probes: "<<probes<<"\n"; 
	cout<<"threshold : "<<threshold<<"\n";
	generateSSCCReport(probes, false);
	//generateGeneratedReport();
        //generateExpandedReport();
	
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
	boost::dynamic_bitset<> b_initial_v(heuristics.size()); 
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


        int max_h_initial_value = 0;
        for (size_t i = 0; i < h_initial_v.size(); i++) {
            int a = h_initial_v.at(i);
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

	//Initializing the operators to determinate how many childs the root node will have
	std::vector<const GlobalOperator*> applicable_ops0;
        GlobalState global_state0 = g_state_registry->lookup_state(g_initial_state().get_id());
        g_successor_generator->generate_applicable_ops(global_state0, applicable_ops0);
	//count nodes generated
	double amount_initial = (double)applicable_ops0.size();

	collector.insert(std::pair<boost::dynamic_bitset<>, double>(b_initial_v, 1 + amount_initial));

	//end counting nodes generated by the root node
	
        SSNode node;
        const GlobalState &initial_state = g_initial_state();
        StateID initial_state_id = initial_state.get_id();
        node.set_id(initial_state_id);
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
 
	long queue_counter = 0; 
	while( !queue.empty() )
	{
		queue_counter++;
		if (queue_counter%1000==0) {
			if(search_time()>sampling_time_limit||g_timer()>overall_time_limit){
				cout<<"Search_timer past maximum sampling_time"<<endl;
	      			cout<<"selecting best heuristic after search_time: "<<search_time()<<", seconds,g_timer:"<<g_timer()<<endl;
				generateSSCCReport(threshold, true);
	      			//select_best_heuristics_greedy();
	      			exit(0);
			}
		}
#ifdef _SS_DEBUG
	  cout<<"queue.size:"<<queue.size()<<endl;
#endif
		Type out = queue.begin()->first;
		SSNode s = queue.begin()->second;

               	int g_real =  s.getGreal();
                int level = out.getLevel();
		double w = s.getWeight(); 
                //std::vector<int> h_global_v = s.getHC();
		int min_h=s.getH();


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
		vweight.push_back(s);
		
	        
                //Insert each node.
                Node2 node(min_h + g_real, level);
                double new_f_value = min_h + g_real;
                //count nodes expanded
                if (new_f_value <= threshold) {
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
                	//end count node
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

			//vector<int> h_child_v;
                  	boost::dynamic_bitset<> b_child_v(heuristics.size());b_child_v.set();
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
				  //h_child_v.push_back(INT_MAX/2);
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
				  if(domination_check){
				    F_culprit.push_back(heuristics[i]->get_heuristic() + g_real + get_adjusted_cost(*op) );
				  }
				  //h_child_v.push_back(new_heur);
				  //string heur_name = heuristics[i]->get_heur_name();
				  //heur_name_v.push_back(heur_name);
				}
			}

                  	//cout<<"dump beging"<<endl;
                  	//child.dump_inline();
                  	//cout<<"dump end"<<endl;
                  	//heur_name_v.clear();
			//  cout<<", g_real = "<<g_real + get_adjusted_cost(*op)<<" f_min = "<< h + g_real + get_adjusted_cost(*op)<<",f_max:"<<max_h + g_real + get_adjusted_cost(*op)<<",b_child_v.count:"<<b_child_v.count()<<endl;
#ifdef _SS_DEBUG
			  cout<<", g_real = "<<g_real + get_adjusted_cost(*op)<<" f_min = "<< h + g_real + get_adjusted_cost(*op)<<",b_child_v.count:"<<b_child_v.count()<<endl;
			  cout<<"\tget_adjusted_cost(*op) = "<<get_adjusted_cost(*op)<<"\n";
			  cout<<"\tChild_"<<(i+1)<<" : h = "<<h<<",b_child_v:"<<b_child_v<<endl; 
                        /*  for (size_t i = 0; i < h_child_v.size(); i++) {
                            int h_value = h_child_v.at(i);
			      cout<<h_value + g_real + get_adjusted_cost(*op);
 		            if (i != h_child_v.size() -1) {
                         	cout<<"/";
                      	    }
                        }*/
			  cout<<h + g_real + get_adjusted_cost(*op)<<endl;
			  cout<<", level = "<<(level + 1);
			  cout<<", w = "<<w<<"\n";
#endif

	     		std::vector<const GlobalOperator *> applicable_ops_2;
             

             		GlobalState global_state_2 = g_state_registry->lookup_state(child.get_id());
			//cout<<"S:"<<endl;global_state_2.dump_inline();
             		g_successor_generator->generate_applicable_ops(global_state_2, applicable_ops_2);
             
             		int amount = applicable_ops_2.size();
                
             		std::pair<std::map<boost::dynamic_bitset<>, double>::iterator, bool> ret2;
          
             		std::map<boost::dynamic_bitset<>, double>::iterator it2; 
			if(domination_check){
			  //cout<<"inserting F_culprit:"<<F_culprit<<endl;
			  F_culprits.insert(F_culprit);
			}
      
                        //add to the collector
             		if (b_child_v.count() > 0) {
 	     			ret2 = collector.insert(std::pair<boost::dynamic_bitset<>, double>(b_child_v, amount*w));
             			it2 = ret2.first;

             			if (ret2.second) {
					//cout<<"raiz bc new is added"<<endl;
             			} else {
                			//cout<<"raiz bc old is being updated"<<endl; 
                			it2->second += amount*w;
                			//cout<<", newcc : "<<it2->second<<"\n"; 
             			}

                        //Make pruning
			   Type object = sampler->getType(child.get_id(), h, 1);
			   
                           object.setLevel( level + 1 );
                           
                           SSNode child_node;
                           StateID child_state_id = child.get_id();
                           child_node.set_id(child_state_id);
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

void SSSearch::generateGeneratedReport() {
       double count_nodes = 0;
       for (map<Node2, double>::iterator iter = generated.begin(); iter != generated.end(); iter++) {
           double n = iter->second;
           count_nodes += n;
       }
       cout<<"count nodes generates : "<<count_nodes/(double)ss_probes<<endl;
       generated.clear();
}

void SSSearch::generateExpandedReport() {
        
        string dominio = domain_name;
        string tarefa = problem_name2;
        string heuristica = heuristic_name2;
        double total_nodes = 0.0;
        for (map<Node2, double>::iterator iter = expanded.begin(); iter != expanded.end(); iter++) {
            double n = iter->second;
            total_nodes += n;
        }
	double predictionExpanded = (double)total_nodes/(double)ss_probes;
        cout<<"count nodes expanded : "<<predictionExpanded <<endl;

        cout<<"dominio = "<<dominio<<endl;
        cout<<"tarefa = "<<tarefa<<endl;
        cout<<"heuristica = "<<heuristica<<endl;

	string dirDomain, dirfDist, outputFile;

	if (is_mov_bound) {
		/*if (ss_probes == 1000) {
                        dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/reportss_bounds/"+dominio;
                        dirfDist = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/reportss_bounds/"+dominio+"/fdist";
                        outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/reportss_bounds/"+dominio+"/fdist/"+tarefa;
                } else {*/
                        string nameProbes = "reportss_bounds";
                        nameProbes += "_probes_";
                        nameProbes += boost::lexical_cast<std::string>(ss_probes);
                        cout<<"nameProbes = "<<nameProbes<<"\n";
                        dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes  +"/"+dominio;
                        dirfDist = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes  +"/"+dominio+"/fdist";
                        outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes  +"/"+dominio+"/fdist/"+tarefa;
                //}
	} else {
		string nameProbes = "reportss_";
                nameProbes += boost::lexical_cast<std::string>(ss_probes);
                nameProbes += "_probes";
                cout<<"nameProbes = "<<nameProbes<<"\n";

        	dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio;
        	dirfDist = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio+"/fdist";
        	outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio+"/fdist/"+tarefa;
	}

        ofstream output;

        output.open(outputFile.c_str());
        output<<"\t"<<outputFile.c_str()<<"\n";
        //output<<"predictionSS: "<<totalPrediction<<"\n";
        output<<"predictionSS: "<<predictionExpanded<<"\n";
        output<<"ss_timer: "<<ss_timer_value<<"\n";
        

        if (system(dirDomain.c_str())) {
           cout<<"Directory: "<<heuristica<<" created."<<endl;
        }

        if (system(dirfDist.c_str())) {
           cout<<"Directory: fdist created."<<endl;
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
                    q.push_back(((double)iter->second)/(double)ss_probes);

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
        expanded.clear();
}


double SSSearch::getProbingResult() {
        double expansions = 0;
        
        for (size_t i = 0; i < vweight.size(); i++) {
             SSNode n = vweight.at(i);
             expansions += n.getWeight();
        }
	//cout << endl;
        //cout<<"expansions = "<<expansions<<endl;
        return expansions;
}

void SSSearch::generateSSCCReport(int n_probes, bool termination) {

        string dominio = domain_name;
        string tarefa = problem_name2;
        string heuristica = heuristic_name2;
	string domain_pddl = domain_instance_pddl;

        cout<<"dominio = "<<dominio<<endl;
        cout<<"tarefa_pddl = "<<tarefa<<endl;
	cout<<"domain_pddl = "<<domain_pddl<<"\n";
        cout<<"heuristica = "<<heuristica<<endl;
        //int length = tarefa.size();
        //cout<<"length = "<<length<<endl;
        size_t found = tarefa.find(".");
        //cout<<"found = "<<found<<endl;
        string name = tarefa.substr(0, found);
        name+="_F_";
        name+=boost::lexical_cast<std::string>(threshold);
        name += ".csv";
        //cout<<"name = "<<name<<endl;
	
	string dirDomain, dirSSCC, outputFile;
	if (is_mov_bound) {
		string nameProbes = "reportss_bounds";
                nameProbes += "_probes_";
                nameProbes += boost::lexical_cast<std::string>(ss_probes);
                cout<<"nameProbes = "<<nameProbes<<"\n";

        	dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio;
        	dirSSCC = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio+"/bc";
        	outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio+"/bc/"+name;
	} else {
		string nameProbes = "reportss_";
                nameProbes += boost::lexical_cast<std::string>(ss_probes);
                nameProbes += "_probes";
                cout<<"nameProbes = "<<nameProbes<<"\n";

		dirDomain = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio;
        	dirSSCC = "mkdir "+_HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio+"/bc";
        	outputFile = _HOME_INFO+"/marvin/marvin/testss/"+heuristica+"/"+ nameProbes +"/"+dominio+"/bc/"+name;
	}

        if (system(dirDomain.c_str())) {
           cout<<"Directory: "<<heuristica<<" created."<<endl;
        }

        if (system(dirSSCC.c_str())) {
           cout<<"Directory: SSCC created."<<endl;
        }

        ofstream output;
        output.open(outputFile.c_str());

	//print the name of the heuristics just to be analyzed later
        for (size_t i = 0; i < heuristics.size(); i++) {
            string heur_name = heuristics[i]->get_heur_name();
            output<<"h(,"<<i<<"):,"<<heur_name<<"\n";
        }

	int n_heuristics = heuristics.size();
	int count_line = collector.size();
	//cout<<"n_heuristics="<<n_heuristics<<"\n";
	//cout<<"count_line="<<count_line<<"\n";

	//get the combintation 1/0/1/0..../1/1
	int** harray = new int*[count_line];
	double** ccarray = new double*[count_line];
	for (int i = 0; i < count_line; i++) {
		harray[i] = new int[n_heuristics];
		ccarray[i] = new double[1];
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
			harray[counter_line][i] = b_node_v.test(i);
                        if (i != b_node_v.size() - 1) {
                                //cout<<"/";
                                output<<"/";
                        }
                }
                //cout<<")cc="<<(double)cc/(double)n_probes<<"\n";
		double result_cc = (double)cc/(double)n_probes; 
                //output<<")cc="<<(double)cc/(double)n_probes<<"\n";
                output<<")cc="<<result_cc<<"\n";
		ccarray[counter_line][0] = result_cc;
		counter_line++;
        }
        output.close();

	if (termination) {
		//make it work in 30 minutes
		string delimiter = ",";
		//cout<<"heuristic-information\n";
		map<string, double> add_line_map_heuristic;
		map<string, vector<string> > map_info_heur;
		for (size_t i = 0; i < heuristics.size(); i++) {
			vector<string> collector;
                	string s =  heuristics[i]->get_heur_name();
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
                	for (int j = 0; j < count_line; j++) {
                        	if (harray[j][i] == 1) {
                                	sum_ones += ccarray[j][0];
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
        	//cout<<"ending m:\n";


        	vector<string> v_gapdb_string;
		string heuristic_good = "gapdb_good";

		bool run_min_heuristic = true;
		int counter_just_ga_heur = 0;
		int total_gapdb_heuristics = getTotalGAHeurs(number_gapdb_heurs);
        	map<string, vector<string> >::iterator iter;
        	for (iter = map_info_heur.begin(); iter != map_info_heur.end(); iter++) {
                	string gapdb_string; //the name of each heuristic, just remember that the fd only support gapdb and do not gapdb_deep or gapdb_good => both need to be changed to gapdb
                	string s = iter->first;
                	vector<string> info = iter->second;

			if (run_min_heuristic) {
                		if (s == min_number_heuristic) {
                			gapdb_string =  heuristic_good + "(mp=";
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
                        		//cout<<"\tgapdb_string = "<<gapdb_string<<"\n\n";

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
                                		v_gapdb_string.push_back(heur_blind);
                        		} else {
                                		v_gapdb_string.push_back(gapdb_string);
                        		}
                        		//cout<<"gapdb_string = "<<gapdb_string<<"\n";
                		}// s == min_number_heuristic
			} else { //end run_min_heuristic
				if (isGAPDB(s)) {
                			gapdb_string = "gapdb(mp=";
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
                        		//cout<<"counter_just_ga_heur = "<<counter_just_ga_heur<<"\n";
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
                                		v_gapdb_string.push_back(heur_blind);
                        		} else {
                                		v_gapdb_string.push_back(gapdb_string);
                        		}
                        	//cout<<"gapdb_string = "<<gapdb_string<<"\n";
				}//end isGAPDB
			}//end run_min_heuristic
        	}
        	//cout<<"v_gapdb_string.size() = "<<v_gapdb_string.size()<<"\n";
        	//end astar_gpdb call the bc from ss

		string PROB_GOOD = "problemas_";
                PROB_GOOD += boost::lexical_cast<std::string>(ss_probes);
                PROB_GOOD += "_probes";
                //cout<<"PROB_GOOD = "<<PROB_GOOD<<"\n";
		string ASTAR_GOOD_NAME = "_SS_ASTAR";
		int deep_F_boundary = threshold;

		//create directories, running individual heuristic and running all gapdb heuristics using the same heuristic_good

		string dirProbGood = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic_good+"/";
                if (system(dirProbGood.c_str())) {
                	cout<<"create directory "<<dirProbGood.c_str()<<"\n";
                }

                string dirProblema = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic_good+"/" + PROB_GOOD;
                if (system(dirProblema.c_str())) {
                	cout<<"create directory "<<dirProblema.c_str()<<"\n";
                }

		string dirResultado = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic_good+"/reportastar";
                if (system(dirResultado.c_str())) {
                	cout<<"create directory "<<dirResultado.c_str()<<"\n";
                }

		string pastaProblema = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic_good+"/" + PROB_GOOD  + "/"+dominio;
		if (system(pastaProblema.c_str())) {
			cout<<"create directory "<<pastaProblema.c_str()<<"\n";
		}

		string pastaResultado = "mkdir "+_HOME_INFO+"/marvin/marvin/astar/"+heuristic_good+"/" + PROB_GOOD  +  "/"+dominio+"/resultado";
		if (system(pastaResultado.c_str())) {
			cout<<"create directory "<<pastaResultado.c_str()<<"\n";
		}

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

                        		//cout<<"pot[0] = "<<pot[0]<<"\n";
                        		//cout<<"pot[1] = "<<pot[1]<<"\n";
                        		//cout<<"pot[2] = "<<pot[2]<<"\n";
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
                		string new_problem_name = problem_name2; //--problem_name == problema.c_str()
                		string t = new_problem_name;
                		size_t found = t.find(".");
                		string new_problem_name_mod = t.substr(0, found);

				string prob_name_gapdb = new_problem_name_mod + "_gapdb_" + final_number_heur  + ".pddl";

				//end get real name
				//create the directory of the problemas_500_probes_good

				string dirSASPLAN = "mkdir "+_HOME_INFO+"/marvin/fd/plan/";
                		if (system(dirSASPLAN.c_str())) {
                			cout<<"create directory "<<dirSASPLAN.c_str()<<"\n";
                		}

				string dirSASPLANDomain = "mkdir "+_HOME_INFO+"/marvin/fd/plan/"+dominio;
                		if (system(dirSASPLANDomain.c_str())) {
                			cout<<"create directory "<<dirSASPLANDomain.c_str()<<"\n";
                		}

                		//creation of each sh file for the gapdb heuristic
                		string arquivo;
                		string sas;
                		stringstream Resultado;

                		Resultado<<i+1;

                		arquivo = new_problem_name_mod + "_gapdb_" + final_number_heur + ".sh";
                		arquivo = "/" + arquivo;
                		arquivo = dominio + arquivo;
                		arquivo = "astar/"+heuristic_good+"/" + PROB_GOOD  +  "/" + arquivo;
                		arquivo = "marvin/" + arquivo;
                		arquivo = "marvin/"+ arquivo;
                		arquivo =  _HOME_INFO+"/" + arquivo;
				cout<<"arquivo="<<arquivo<<"\n";
                		ofstream outfile(arquivo.c_str(), ios::out);

				string newdominio = dominio + "_" + final_number_heur + "_" + new_problem_name_mod;

                		sas = "Inside_Astar";
                		sas += newdominio;

				sas += Resultado.str();
                		//End creation of each sh file for the gapdb heuristic

                		string parameter =  final_real_heur;//v_gapdb_string.at(i);
                		cout<<"parameter_"<<i<<" = "<<parameter<<"\n";

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

				//cout<<"pasta = "<<dominio<<"\n\n;
				outfile<<"FD_ROOT="<<_HOME_INFO<<"/marvin/fd\n\n";
        			outfile<<"TEMP="<<_HOME_INFO<<"/marvin/fd/temp\n\n";
        			outfile<<"DIR=$(mktemp  --tmpdir=${TEMP})\n\n";	
                		outfile<<"RESULTS="<<_HOME_INFO<<"/marvin/marvin/astar/"<<heuristic_good<<"/" + PROB_GOOD  +  "/"<<dominio<<"/resultado\n\n";
				outfile<<"cd ${DIR}\n\n";
                		outfile<<"python3 ${FD_ROOT}/src/translate/translate.py ${FD_ROOT}/benchmarks/"<<dominio<<"/"<<domain_pddl<<" ${FD_ROOT}/benchmarks/"<<dominio<<"/"<<tarefa<<"\n\n";

                		outfile<<"${FD_ROOT}/src/preprocess/preprocess < output.sas"<<"\n\n"; 
				//Santiago's code to find the F_boundary on the fly
                		outfile<<"${FD_ROOT}/src/search/downward-release --use_saved_pdbs --domain_name "<<dominio<<" --problem_name "<<tarefa<<" --heuristic_name "<<heuristic_good<<" --problem_name_gapdb "<<prob_name_gapdb<<" --deep_F_boundary "<<deep_F_boundary<<"  --search \"astar_original("<<parameter<<")\" <  output > ${RESULTS}/"<<prob_name_gapdb<<"\n\n";
                	

				outfile<<"\n\nrm ${DIR}\n\n";
                		outfile<<"\n\nmv sas_plan ${FD_ROOT}/plan/"<<dominio<<"/"<<tarefa<<"\n\n";

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
			}//v_gapdb_string for loop
		} else { //end else run_min_heuristic
			string dirSASPLAN = "mkdir "+_HOME_INFO+"/marvin/fd/plan_"+heuristic_good+"/";
                	if (system(dirSASPLAN.c_str())) {
                		cout<<"create directory "<<dirSASPLAN.c_str()<<"\n";
                	}

			string dirSASPLANDomain = "mkdir "+_HOME_INFO+"/marvin/fd/plan_"+heuristic_good+"/"+dominio;
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

			//create new variable called deep_F_boundary
        		//cout<<"deep_F_boundary = "<<deep_F_boundary<<"\n";

        		//begin
			
                	string new_problem_name = problem_name2; //--problem_name == problema.c_str()
        		string t = new_problem_name;
        		size_t found = t.find(".");
        		string new_problem_name_mod = t.substr(0, found);
        		//cout<<"new_problem_name_mod = "<<new_problem_name_mod<<"\n";
        		//stringstream number;
        		//number<<i; //this should contains the real number
        		//name that will be used in the backend
        		//string prob_name_gapdb = new_problem_name_mod + "_gapdb_" + number.str() + ".pddl";
        		string prob_name_gapdb = new_problem_name_mod + "_gapdb_all.pddl";
        		//cout<<"prob_name_gapdb = "<<prob_name_gapdb<<"\n\n\n";
        		//end

        		//end get real name
        		//creation of each sh file for the gapdb heuristic
        		string arquivo;

        		arquivo = new_problem_name_mod + "_gapdb_all.sh";
			arquivo = "/" + arquivo;
        		arquivo = dominio + arquivo;
        		arquivo = "astar/"+heuristic_good+"/" + PROB_GOOD  +  "/" + arquivo;
        		arquivo = "marvin/" + arquivo;
        		arquivo = "marvin/"+ arquivo;
        		arquivo =  _HOME_INFO+"/" + arquivo;
        		ofstream outfile(arquivo.c_str(), ios::out);

        		string parameter =  heuristic_generator;

        		//Begin construction of the sh file
        		outfile<<"#!/bin/bash\n\n";
        		outfile<<"#PBS -N "<<ASTAR_GOOD_NAME<<"\n\n#PBS -m a\n\n#PBS -l walltime=00:30:00\n\n";
        		outfile<<"#PBS -M marvin.zarate@ufv.br\n\n";
        		outfile<<"cd $PBS_O_WORKDIR\n\n";
        		outfile<<"source /usr/share/modules/init/bash\n\n";
        		outfile<<"module load python\nmodule load mercurial\n\n";

        		outfile<<"FD_ROOT="<<_HOME_INFO<<"/marvin/fd\n\n";
        		outfile<<"TEMP="<<_HOME_INFO<<"/marvin/fd/temp\n\n";
        		outfile<<"DIR=$(mktemp  --tmpdir=${TEMP})\n\n";
        		//cout<<"pasta = "<<pasta.c_str()<<"\n\n";

        		outfile<<"RESULTS="<<_HOME_INFO<<"/marvin/marvin/astar/"<<heuristic_good<<"/" + PROB_GOOD  +  "/"<<dominio<<"/resultado"<<"\n\n";
        		//outfile<<"cd "<<_HOME_INFO<<"/marvin/fd\n\n";
        		outfile<<"cd ${DIR}\n\n";
                	outfile<<"python3 ${FD_ROOT}/src/translate/translate.py ${FD_ROOT}/benchmarks/"<<dominio<<"/"<<domain_pddl<<" ${FD_ROOT}/benchmarks/"<<dominio<<"/"<<tarefa<<"\n\n";

        		outfile<<"${FD_ROOT}/src/preprocess/preprocess < output.sas"<<"\n\n";

        		outfile<<"${FD_ROOT}/src/search/downward-release --use_saved_pdbs --domain_name "<<dominio<<" --problem_name "<<tarefa<<" --heuristic_name "<<heuristic_good<<" --problem_name_gapdb "<<prob_name_gapdb<<" --deep_F_boundary "<<deep_F_boundary<<"  --search \"astar_original(max(["<<parameter<<"]))\" <  output > ${RESULTS}/"<<prob_name_gapdb<<"\n\n";

        		outfile<<"\n\nrm ${DIR}\n\n";
			outfile<<"\n\nmv sas_plan ${FD_ROOT}/plan_"<<heuristic_good<<"/"<<dominio<<"/"<<tarefa<<"\n\n";

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
                		executeFile = "timeout 1800 sh "+arquivo; //setting the limit time
				if(system(executeFile.c_str())) {
					cout<<"running in the local...\n";
				}
        		}
		}//end run_min_heuristic
	}//end termination
}

void SSSearch::printQueue() {
        cout<<"\nPrintQueue\n";
	for (map<Type, SSNode>::iterator iter = queue.begin(); iter !=  queue.end(); iter++) {
            Type t = iter->first;
            SSNode t2  = iter->second;
            cout<<"\t\t h = "<<t.getH()<<" g = "<<t.getLevel()<<" f = "<<t.getH() + t2.getGreal()<<" w = "<<t2.getWeight()<<"\n"; 
        }
        cout<<"\n";
        cout<<"\nEnd PrintQueue\n";
}

void SSSearch::select_best_heuristics_greedy(){
  //h_comb_current_time.clear();
  vector<Heuristic *> selected_heuristics;
  cout<<"Starting select_best_h_combs_RIDA_Sampling_Greedy,bcs:"<<collector.size()<<endl;
  
  //double best_current_time=0;
  //long best_initial_nodes=0;
  boost::dynamic_bitset<> current_h_comb(heuristics.size());
  boost::dynamic_bitset<> best_h_comb(heuristics.size());
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
  

  if(best_h_comb.count()==0){
    vector<double> initial_heur_sizes(heuristics.size(),0);
    for (map<boost::dynamic_bitset<>,double>::const_iterator it_h_values= collector.begin(); it_h_values != collector.end();it_h_values++){
      for(size_t heur=0;heur<heuristics.size();heur++){
	//cout<<"it_h_values:"<<it_h_values->first<<endl;
	if(it_h_values->first.test(heur)){//if indiv heuristic is active and culling the counter then this counter does not apply to the h_comb
	  initial_heur_sizes[heur]+=it_h_values->second;
	  //cout<<"added "<<it_h_values->second<<" to initial_heur_"<<heur<<endl;
	}
	else{
	  //cout<<"skipped added "<<it_h_values->second<<" to initial_heur_sizeds of "<<heur<<endl;

	}
      }
    }
    //Now calculate which one was best
    std::cout << "scientific:\n" << std::scientific;
    for(size_t heur=0;heur<heuristics.size();heur++){
      cout<<"initial_heur_sizes("<<heur<<"):"<<initial_heur_sizes.at(heur)<<endl;
      if(initial_heur_sizes.at(heur)<best_comb_nodes){
	best_comb_nodes=initial_heur_sizes.at(heur);
	best_h_comb.resize(heuristics.size());best_h_comb.reset();best_h_comb.set(heur);
	cout<<"temp best_comb_nodes:"<<best_comb_nodes<<",temp best heur:";print_h_comb(best_h_comb);cout<<endl;
    }
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
    cout<<"active_heurs"<<active_heurs<<endl;
    cout<<"candidate_heurs:"<<candidate_heurs<<endl;

    vector<double> h_comb_prunning_nodes(heuristics.size(),0);
    for (map<boost::dynamic_bitset<>,double>::iterator it_h_values= collector.begin(); it_h_values != collector.end();){
      delete_counter=false;
      for(size_t heur=0;heur<active_heurs.size();heur++){
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
	cout<<"best_h_comb:";print_h_comb(best_h_comb);cout<<",nodes:"<<best_comb_nodes<<",extra pruned:"<<highest_prunner<<endl;
      }
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

  cout<<"progressive_greedy_timer:,"<<progressive_greedy_timer<<",level:,"<<best_h_comb.count()<<",Final Overall best h_comb(translated):";print_h_comb(best_h_comb);
  for(size_t j=0;j<best_h_comb.size();j++){
    if(best_h_comb.test(j)){
      cout<<heuristics[j]->get_heur_call_name()<<endl;
      //selected_heuristics.push_back(selectable_heuristics.at(j));
    }
  }
        
  ofstream output;
  string output_file="temp_scripts/selected_heurs_";
  output_file+=problem_name2;
  output_file+=".sh";
  output.open(output_file.c_str());
  output<<"./downward-release --use_saved_pdbs --domain_name "<<domain_name<<" --problem_name "<<problem_name2<<" --search \"astar(";
  size_t populated_counter=0;
  output<<"min([";
  for(size_t j=0;j<best_h_comb.size();j++){
    if(best_h_comb.test(j)){
      output<<heuristics[j]->get_heur_call_name();
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

void SSSearch::initialize() {	
	cout << "SSSearch ..."<<endl;
	search_time.reset();
	level_time.reset();

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
