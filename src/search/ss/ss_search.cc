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

SSSearch::SSSearch(const Options &opts) : SearchEngine(opts), current_state(g_initial_state()) {
	//rg = opts.get<double>("rg");
	//rl = opts.get<double>("rl");
	//lookahead = opts.get<int>("lookahead");
	//beamsize = opts.get<int>("beamsize");
	//maxlevels = opts.get<int>("maxlevels");
	//timelimit = opts.get<int>("timelimit");
	
	ScalarEvaluator * evaluator = opts.get<ScalarEvaluator *>("eval");
	std::set<Heuristic *> hset;
	evaluator->get_involved_heuristics(hset);
	for (set<Heuristic *>::iterator it = hset.begin(); it != hset.end(); it++) {
		heuristics.push_back(*it);
	}
	assert(heuristics.size() == 1);
	heuristic = heuristics[0];

	sampler = new TypeSystem(heuristic);
	this->RanGen2 = new CRandomMersenne((unsigned)time(NULL));        
}

SSSearch::~SSSearch() {
}


SearchStatus SSSearch::step() {
        predict(ss_probes);
        return SOLVED;
}

void SSSearch::predict(int probes) {
        totalPrediction = 0;

        for (int i = 0; i < probes; i++) {
            vweight.clear();
            probe();
            double p = getProbingResult();
            totalPrediction = totalPrediction + (p - totalPrediction)/(i + 1);
            cout<<"**********"<<endl;
            cout<<"p = "<<p<<endl;
            cout<<"prePre = "<<totalPrediction<<endl;
            cout<<"**********"<<endl;
        }
        cout<<"\ntotalPrediction : "<<totalPrediction<<"\n";
        generateGeneratedReport();
        generateExpandedReport();
}

void SSSearch::probe()
{
	/*
	 * Probing is done based on the types of the children of the root state
	 */

        queue.clear();
	// evaluating the initial state
	heuristic->evaluate(g_initial_state());
	if (heuristic->is_dead_end())
	{
		assert(heuristic->dead_ends_are_reliable());
		cout << "Initial state is a dead end." << endl;
		exit(0);
	}
	initial_value = heuristic->get_value();

        //for the open domains the heuristic is set to six
        threshold = 2*initial_value;

        SSNode node;
        const GlobalState &initial_state = g_initial_state();
        StateID initial_state_id = initial_state.get_id();
        node.setId(initial_state_id);
        node.setWeight(1.0);
	/*
	 * Seeding the prediction with the children of the start state
	 *
	 */
	Type type = sampler->getType(node.getId(), initial_value, 1);
          
 
	type.setLevel( 0 );

	queue.insert( pair<Type, SSNode>( type, node ) );

        int nraiz = 1;
  
	while( !queue.empty() )
	{
		Type out = queue.begin()->first;
		SSNode s = queue.begin()->second;
               	int g = (int)out.getLevel();

                //printQueue();

                std::map<Type, SSNode>::iterator rt;
                rt = queue.find(out);


		queue.erase( rt );
                //cout<<nraiz<<": Raiz: h = "<<out.getH()<<" g = "<<out.getLevel()<<" f = "<<out.getH() + out.getLevel()<<" w  = "<<s.getWeight()<<endl;   
                nraiz++;                
                
		vweight.push_back(s);
		
	        
                //Insert each node.
                Node2 node2(out.getH() + g, g);

                //count nodes expanded

                std::pair<std::map<Node2, double>::iterator, bool> ret0;

                std::map<Node2, double>::iterator it0;

                ret0 = expanded.insert(pair<Node2, double>(node2, s.getWeight()));
                it0 = ret0.first;

                if (ret0.second) {
                    cout<<"new node expanded is added."<<endl;
                } else {
                    cout<<"node expanded is being updated."<<endl;
                    it0->second += s.getWeight();
                    cout<<"it0->second = "<<it0->second<<endl;
                }

                //end count node
 
		int h = -1;
		double w = (double)s.getWeight();
		//cout<<"w = "<<w<<endl;
 
                std::vector<const GlobalOperator*> applicable_ops;
                GlobalState global_state = g_state_registry->lookup_state(s.getId()); 
		g_successor_generator->generate_applicable_ops(global_state, applicable_ops);

                //count nodes generated
                double amount = (double)applicable_ops.size();

                std::pair<std::map<Node2, double>::iterator, bool> ret;
                std::map<Node2, double>::iterator it;

                ret = generated.insert(std::pair<Node2, double>(node2, amount*w));
                it = ret.first;


                if (ret.second) {
                   cout<<"new node is added."<<endl;
                } else {
                   cout<<"old is being updated."<<endl;
                   it->second += amount*w;
                   cout<<"new = "<<it->second<<endl;
                }
                 //end count nodes generated
		for (size_t i = 0; i < applicable_ops.size(); ++i)
		{
                        const GlobalOperator *op = applicable_ops[i];
                        GlobalState child = g_state_registry->get_successor_state(global_state, *op);

                        heuristic->evaluate(child);	

			if(!heuristic->is_dead_end())
			{
				h = heuristic->get_heuristic();

			}	
                        
			//cout<<"\tChild: h = "<< h <<" g = "<< g + 1 <<" f = "<< h + g + 1 <<" w = "<<w<<endl; 

                        if (h + g + 1 <= threshold) {
			   Type object = sampler->getType(child.get_id(), h,  1);
			    
                           object.setLevel( g + 1 );

                           SSNode child_node;
                           StateID child_state_id = child.get_id();
                           child_node.setId(child_state_id);
                           child_node.setWeight(w);

                           //cout<<"\tChild: h = "<<object.getH()<<" g = "<<object.getLevel()<<" f = "<<object.getH() + object.getLevel()<<"\n"; 
			   

                           map<Type, SSNode>::iterator queueIt = queue.find( object );
			   if( queueIt != queue.end() )
			   {

                                //cout<<"\tis duplicate: h = "<<queueIt->first.getH()<<" g = "<<queueIt->first.getLevel()<<" f = "<< queueIt->first.getH()   +  queueIt->first.getLevel()<<"\n";
                                SSNode snode = queueIt->second;
				double wa = (double)snode.getWeight();
				//snode.setWeight( wa + w);
                                queueIt->second.setWeight(wa + w); 
                                //std::pair<std::map<Type, SSNode>::iterator, bool> ret0;

                                //ret0 = queue.insert(pair<Type, SSNode>(object, snode));
                                //cout<<"\tsnode.getWeight() = "<<snode.getWeight()<<endl;
                                //queueIt->second.setWeight(snode.getWeight());
 
 
				double prob = ( double )w / (double)( wa + w );
				int rand_100 =  RanGen2->IRandom(0, 99);  //(int)g_rng.next(100);
                          	 
                                double a = (( double )rand_100) / 100;
                                //cout<<"a = "<<a<<" prob = "<<prob<<endl;
                                
                                
                                
                                
				if (a < prob) 
				{
                                        //cout<<"\t\tAdded even though is duplicate.\n";
                                        
				        child_node.setWeight( wa + w);
                                        
                                        std::pair<std::map<Type, SSNode>::iterator, bool> ret;
                                     	queue.erase(object); 

                                        ret = queue.insert( pair<Type, SSNode>( object, child_node ));      

                                        queueIt = ret.first;
                                        queueIt->second.setWeight(child_node.getWeight());
                                        
                                      	
				} else {
                                        //cout<<"\t\tNot added.\n";
                                }
			   } 
			   else
			   {
                                //cout<<"\t\tNew node added\n";
				queue.insert( pair<Type, SSNode>( object, child_node ) );
                                //cout<<"\t\tchild_node.getWeight() = "<<child_node.getWeight()<<"\n";
                                //S.insert( pair<Type, double>( object,12.44 ) );
                                
                                //cout<<"\t\tChild: h = "<< h <<" g = "<< g + 1 <<" f = "<< h + g + 1 << " threshold: " << threshold <<" w = "<<child_node.getWeight()<<endl;
                           }
                        }
			else 
			{
				//cout << "\tNode was pruned!" << endl;
				//cout<<"\tChild: h = "<< h <<" g = "<< g + 1 <<" f = "<< h + g + 1 << " threshold: " << threshold <<"\n";
			}
		}
	}
        
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
        cout<<"count nodes expanded : "<<(double)total_nodes/(double)ss_probes<<endl;

        cout<<"dominio = "<<dominio<<endl;
        cout<<"tarefa = "<<tarefa<<endl;
        cout<<"heuristica = "<<heuristica<<endl;

        string dirDomain = "mkdir /home/marvin/marvin/testss/"+heuristica+"/reportss/"+dominio;
        string dirfDist = "mkdir /home/marvin/marvin/testss/"+heuristica+"/reportss/"+dominio+"/fdist";
       
        string outputFile = "/home/marvin/marvin/testss/"+heuristica+"/reportss/"+dominio+"/fdist/"+tarefa;

        ofstream output;

        output.open(outputFile.c_str());
        output<<"\t"<<outputFile.c_str()<<"\n" ;
        output<<"predictionSS: "<<totalPrediction<<"\n";
        output<<"threshold: "<<threshold<<"\n";
        

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
	cout << endl;
        cout<<"expansions = "<<expansions<<endl;
        return expansions;
}

void SSSearch::printQueue() {
        cout<<"\nPrintQueue\n";
	for (map<Type, SSNode>::iterator iter = queue.begin(); iter !=  queue.end(); iter++) {
            Type t = iter->first;
            SSNode t2  = iter->second;
            cout<<"\t\t h = "<<t.getH()<<" g = "<<t.getLevel()<<" f = "<<t.getH() + t.getLevel()<<" w = "<<t2.getWeight()<<"\n"; 
        }
        cout<<"\n";
        cout<<"\nEnd PrintQueue\n";
}


void SSSearch::report_progress()
{
	cout << "h_min: " << total_min << " depth: " << depth << " #states: " << queue.size() << " time: " << search_time << endl;
}

void SSSearch::initialize() {
	cout << "SSSearch ..." << endl;
	search_time.reset();
	level_time.reset();
        
	//queue.clear();
	// evaluating the initial state
	heuristic->evaluate(g_initial_state());
	if (heuristic->is_dead_end())
	{
		assert(heuristic->dead_ends_are_reliable());
		cout << "Initial state is a dead end." << endl;
		exit(0);
	}
	initial_value = heuristic->get_value();
	total_min = initial_value;
       
	progress = true;
	cout << "Initial heuristic value: ";
	cout << initial_value << endl;
	depth = 0;
	report_progress();
	depth ++;
}

static SearchEngine *_parse(OptionParser &parser) {
        
	parser.add_option<ScalarEvaluator *>("eval");
	//parser.add_option<double>("rg", DEFAULT_SS_RG, "the global restarting rate");
	//parser.add_option<double>("rl", DEFAULT_SS_RL, "the local restarting rate");
	//parser.add_option<int>("lookahead", DEFAULT_SS_LOOKAHEAD, "lookahead that defines the type system being used");
	//parser.add_option<int>("beamsize", DEFAULT_SS_BEAMSIZE, "maximum number of nodes expanded by level");
	//parser.add_option<int>("maxlevels", DEFAULT_SS_MAXLEVELS, "maximum number of nodes expanded by level");
	//parser.add_option<int>("timelimit", DEFAULT_SS_TIMELIMIT, "time limit in seconds to solve the problem");
	SearchEngine::add_options_to_parser(parser);
	Options opts = parser.parse();
	SSSearch *engine = 0;
	if (!parser.dry_run()) {
		engine = new SSSearch(opts);
	}
	return engine;
}


static Plugin<SearchEngine> _plugin("ss", _parse);
