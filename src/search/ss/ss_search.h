#ifndef STRATIFIED_SAMPLING_H_
#define STRATIFIED_SAMPLING_H_


#include "../global_state.h"
#include "../heuristic.h"
#include "../global_operator.h"
#include "../search_engine.h"
#include "type.h"
#include "type_system.h"
#include "node2.h"

#include <algorithm>

#include <map>
#include <queue>
#include <set>
#include <deque>

#include "../randomc/randomc.h"
#include "../randomc/mersenne.cpp"
#include "../state_id.h"
#include "../timer.h"

//ss+culprits
#include "../ext/boost/dynamic_bitset.hpp"

class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;


using namespace std;

template <typename T1, typename T2>
struct less_second {
    typedef pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second < b.second;
    }
};

template <typename T1, typename T2>
struct greater_second {
    typedef pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second > b.second;
    }
};

class SSNode{
private:
	StateID id;
	double weight;
        int g_real;
	int h;
	bool lmcut_expanding;//for early term
public:
        SSNode(): id(StateID::no_state), weight(0.0), g_real(0),lmcut_expanding(false) {}
        SSNode(StateID identifier, double w, int g) : id(identifier), weight(w), g_real(g){}
        StateID get_id() const {return this->id;}
	void setId(StateID identifier) {this->id = identifier;}
        double getWeight()  {return this->weight;}
        void setWeight(double w) {this->weight = w;}
        int getGreal() const {return this->g_real;}
        void setGreal(int g) {this->g_real = g;}
	int getH() {return this->h;}
	void setH(int H) {this->h = H;}
	int get_lmcut_expanding() {return this->lmcut_expanding;}
	void set_lmcut_expanding(bool status=true) {lmcut_expanding=status;}
};

class SSQueue {
private:
	SSNode node;
	Type type;
public:
	SSNode getNode() const {return this->node;}
	void setNode(SSNode n) {this->node = n;}
	Type getT() const {return this->type;}
	void setT(Type t) {this->type = t;}
};

struct classcomp {
        bool operator() (const SSQueue& lhs, const SSQueue& rhs) const {
		return lhs.getNode().get_id() < rhs.getNode().get_id(); 
        }
};

struct classcomp2 {
        bool operator() (const SSNode& lhs, const SSNode& rhs) const {
	  //cout<<"\t\t\tCalling classcomp2"<<endl;fflush(stdout);
	  //cout<<"\t\t\tlhs.get_id():"<<lhs.get_id();fflush(stdout);
	  //cout<<"\t\t\trhs.get_id():"<<rhs.get_id();fflush(stdout);
		return lhs.get_id() < rhs.get_id(); 
        }
};

#endif /*MRW_H_*/

class SSSearch : public SearchEngine { 
private:


	map<Type, SSNode> queue;
 
        vector<double> vweight;
	vector<SSQueue> vmeanheur;
        std::map<Node2, double> expanded;

        std::map<Node2, double> generated;
        double totalPrediction;
	double totalPredictionMean;

        std::vector<Heuristic*> heuristics;
	std::vector<Heuristic*> lmcut_heuristic; 
	std::vector<Heuristic*> all_heuristics; 
	Heuristic* heuristic;
	
	
	int initial_value;
        
        GlobalState current_state;
	Timer search_time;
	Timer level_time; //time required to expand an entire level
	Timer ss_timer;
	double ss_timer_value;

	TypeSystem * sampler;

        CRandomMersenne* RanGen2;

	//IDA* - BFS
	std::set<SSQueue, classcomp> L;
	std::set<SSNode, classcomp2> check;

        //ss+culprits
        int threshold;
        map<boost::dynamic_bitset<>,  double> collector;
        map<boost::dynamic_bitset<>,  double> full_collector;
        map<std::vector<int>,  double> collector_heur;
	int last_probe;
	double last_n_expanded;

	//parameters GRHS
	int** harray_grhs;
	double** ccarray_grhs;	
	int count_line_grhs;                

	//parameters SSCC
	int** harray_sscc;
	double** ccarray_sscc;
	int n_heuristics_global2;
	int count_line_sscc;

	//global parameters
	string dominio_global;
	string tarefa_global;
	string heuristica_global;
	string domain_pddl_global;
	int n_heuristics_global;
	int n_probes_global;

protected:

	virtual SearchStatus step();
	virtual void initialize();

public:
	
	SSSearch(const Options &opts);
	virtual ~SSSearch();
        void printQueue(); 
	void printVectorPair(vector<pair<string, double> > vpair);
        void generateExpandedReport();
        void generateGeneratedReport();
	void generateSSCCReport(bool termination);
        double getProbingResult();
	double getMeanHeurResult();
        void probe();
        void predict(int probes);
	int getMinHeur(vector<int> v);
	int getMaxHeur(vector<int> v);
	void select_best_heuristics_greedy();
	void BFS(SSNode root, Type type);
	bool isGAPDB(string heuristic);
        int getTotalGAHeurs(vector<string> v);
	void select_best_heuristics_greedy_stocastic(int cardinality);
	void select_random_greedy(bool termination);
	map<string, double>  heuristicCombinator(bool call_first_time, vector<pair<string, double> > Z_subset, vector<pair<string, double> > Z_full);
	double getSumSubset(vector<pair<string, double> > Z_subset);
	void updateGRHS();
	void updateSSCC();
	void updateGlobalVariables();
};
