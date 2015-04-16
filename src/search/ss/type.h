#ifndef TYPE_H_
#define TYPE_H_

#include <map>

#include "type_children.h"

using namespace std;

class Type {

private:
	TypeChildren children;
	long p;
	unsigned long long h;
	unsigned long long level;
	int best_h;
	double distance;
	int label;
	int random;

public:

	Type();
	Type(long parent, unsigned long long heuristic);

	void addAddtionalInfo(int);
	const TypeChildren& getChildren() const {return children;}
	TypeChildren getConstChildren() const {return children;}
	void setChildren(TypeChildren c){this->children = c;}
	void print() const;

	friend bool operator< (const Type&, const Type&);


        Type &operator=(const Type &rhs);
        bool operator==(const Type &rhs) const;
        //bool operator<(const Type &rhs) const;        


	unsigned long long getH() const {return h;}
	void setH(unsigned long long i) {h = i;}
	int getP() const {return p;}
	void setP(int i) {p = i;}
	unsigned long long getLevel() const {return level;}
	void setLevel(unsigned long long i) {level = i;}
	double getDistance() const {return distance;}
	void setDistance(double i) {distance = i;}
	int getLabel() const {return label;}
	void setLabel(int i) {label = i;}
	int getRandom() const {return random;}
	void setRandom(int r) {this->random = r;}
	int getBestH() const {return best_h;}
	void setBestH(int r) {this->best_h = r;}

	static int lookahead;
};


#endif
