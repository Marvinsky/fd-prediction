#ifndef STRATIFIED_SAMPLING_H_
#define STRATIFIED_SAMPLING_H_

#define INT_MAX   2147483647

#include "../global_state.h"
#include "../heuristic.h"
#include "../global_operator.h"
#include "../search_engine.h"
#include "type.h"
#include "type_system.h"
#include "node2.h"

#include "map"

#include "../randomc/randomc.h"
#include "../randomc/mersenne.cpp"
#include "../state_id.h"
#include "../timer.h"


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
public:
        SSNode(): id(StateID::no_state), weight(0.0), g_real(0){}
        SSNode(StateID identifier, double w, int g) : id(identifier), weight(w), g_real(g){}
        StateID getId() const {return this->id;}
        void setId(StateID identifier) {this->id = identifier;}
        double getWeight()  {return this->weight;}
        void setWeight(double w) {this->weight = w;}
        int getGreal() const {return this->g_real;}
        void setGreal(int g) {this->g_real = g;}
};

class SSSearch : public SearchEngine {
public: 
private:
	
	std::map<int, SSNode> open;
	map<Type, SSNode> queue; 
        vector<SSNode> vweight;
        std::map<Node2, double> expanded;

        std::map<Node2, double> generated;
        double totalPrediction;         

	std::vector<Heuristic*> heuristics; 
	Heuristic* heuristic;

	int initial_value;
        int threshold;
        GlobalState current_state;

        //A* prediction using ss:
	TypeSystem * sampler;
        CRandomMersenne* RanGen2;
        Timer ss_timer;
        double ss_timer_value;	

protected:

	virtual SearchStatus step();
	virtual void initialize();

public:
	SSSearch(const Options &opts);
	virtual ~SSSearch();
        void printQueue(); 
        void generateExpandedReport();
        void generateGeneratedReport();
        double getProbingResult();
        void probe();
        void predict(int probes);
};

#endif /*MRW_H_*/
