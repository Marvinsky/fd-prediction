#ifndef STRATIFIED_SAMPLING_MEAN_H_
#define STRATIFIED_SAMPLING_MEAN_H_

#define INT_MAX   2147483647

#include "../global_state.h"
#include "../heuristic.h"
#include "../global_operator.h"
#include "../search_engine.h"
#include "../ss/type.h"
#include "../ss/type_system.h"
#include "../ss/node2.h"

#include <map>
#include <queue>
#include <set>
#include <deque>

//#include "../randomc/randomc.h"
//#include "../randomc/mersenne.cpp"
#include "../state_id.h"
#include "../timer.h"
#include "../rng.h"

class GlobalOperator;
class Heuristic;
class Options;
class ScalarEvaluator;

using namespace std;

class SSNode{
private:
	StateID id;
	double weight;
        double g_real;
public:
        SSNode(): id(StateID::no_state), weight(0.0), g_real(0){}
        SSNode(StateID identifier, double w, double g) : id(identifier), weight(w), g_real(g){}
        StateID get_id() const {return this->id;}
        void set_id(StateID identifier) {this->id = identifier;}
        double getWeight()  {return this->weight;}
        void setWeight(double w) {this->weight = w;}
        double getGreal() const {return this->g_real;}
        void setGreal(double g) {this->g_real = g;}
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
		return lhs.get_id() < rhs.get_id(); 
        }
};

class SSMean : public SearchEngine { 
private:

	map<Type, SSNode> queue; 
        vector<SSNode> vweight;
	vector<SSQueue> vmeanheur;
        std::map<Node2, double> expanded;
        std::map<Node2, double> generated;
        double totalPrediction;
	double totalPredictionMean;      

	std::vector<Heuristic*> heuristics; 
	Heuristic* heuristic;

	double initial_value;
        double threshold;
        GlobalState current_state;

        //A* prediction using ss:
	TypeSystem * sampler;
        //CRandomMersenne* RanGen3;
        Timer ss_timer;
        double ss_timer_value;
	
	//IDA* - BFS
	std::set<SSQueue, classcomp> L;
	std::set<SSNode, classcomp2> check;

protected:

	virtual SearchStatus step();
	virtual void initialize();

public:
	SSMean(const Options &opts);
	virtual ~SSMean();
        void printQueue(); 
        void generateExpandedReport();
        void generateGeneratedReport();
        double getProbingResult();
	double getMeanHeurResult();
        void probe();
        void predict(int probes);
	//BFS
	void BFS(SSNode root, Type type);
};

#endif /*MRW_H_*/
