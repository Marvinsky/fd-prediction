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
	bool iPDB_expanding;//for early term
	bool lmcut_expanding;//for early term
public:
        SSNode(): id(StateID::no_state), weight(0.0), g_real(0),iPDB_expanding(false),lmcut_expanding(false) {}
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
        std::map<Node2, double> expanded;

        std::map<Node2, double> generated;
        double totalPrediction;         

	std::vector<Heuristic*> heuristics; 
	std::vector<Heuristic*> ipdb_heuristics; 
	std::vector<Heuristic*> lmcut_heuristic; 
	Heuristic* heuristic;
	
	
	int initial_value;
        
        GlobalState current_state;
	Timer search_time;
	Timer level_time; //time required to expand an entire level

	TypeSystem * sampler;

        CRandomMersenne* RanGen2;

	//IDA* - BFS
	std::set<SSQueue, classcomp> L;
	std::set<SSNode, classcomp2> check;

        //ss+culprits
        int threshold;
        map<boost::dynamic_bitset<>,  double> collector;
        map<boost::dynamic_bitset<>,  double> full_collector;
	int last_probe;
                

protected:

	virtual SearchStatus step();
	virtual void initialize();

public:
	
	SSSearch(const Options &opts);
	virtual ~SSSearch();
        void printQueue(); 
        void generateExpandedReport();
        void generateGeneratedReport();
	void generateSSCCReport();
        double getProbingResult();
        void probe();
        void predict(int probes);
	int getMinHeur(vector<int> v);
	void select_best_heuristics_greedy();
	void BFS(SSNode root, Type type);
	void select_best_heuristics_greedy_stocastic(int cardinality);
};
