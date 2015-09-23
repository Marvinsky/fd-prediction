#ifndef STRATIFIED_SAMPLING_H_
#define STRATIFIED_SAMPLING_H_


#include "../global_state.h"
#include "../heuristic.h"
#include "../global_operator.h"
#include "../search_engine.h"
#include "type.h"
#include "type_system.h"
#include "node2.h"

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

class classcomp3 {
	public:
		bool operator()(const vector<int> x, const vector<int> y) {
			bool two_keys_are_different = false;
			if (x.size() == y.size()) {
				for (size_t i = 0; i < x.size(); i++) {
					int hx = x.at(i), hy = y.at(i);
					if (hx != hy) {
						two_keys_are_different = true;
					}
				}
			}  else {
				two_keys_are_different = true;
			}
			return two_keys_are_different;	
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
        map<std::vector<int>,  double, classcomp3> collector_heur;
	int last_probe;
	double last_n_expanded;
                

protected:

	virtual SearchStatus step();
	virtual void initialize();

public:
	
	SSSearch(const Options &opts);
	virtual ~SSSearch();
        void printQueue(); 
        void generateExpandedReport();
        void generateGeneratedReport();
	void generateSSCCReport(bool termination);
        double getProbingResult();
	double getMeanHeurResult();
        void probe();
        void predict(int probes);
	int getMinHeur(vector<int> v);
	void select_best_heuristics_greedy();
	void BFS(SSNode root, Type type);
	bool isGAPDB(string heuristic);
        int getTotalGAHeurs(vector<string> v);
	void select_best_heuristics_greedy_stocastic(int cardinality);
	void select_random_greedy();
};
