#include "node2.h"

Node2::Node2() {
	this->f_value = -1;
	this->level = -1;
}

Node2::Node2(int f, int l) {
	this->f_value = f;
	this->level = l;
}

int Node2::getF()  {
	return this->f_value;
}

void Node2::setF(int f) {
	this->f_value = f;
}

int Node2::getL() {
	return this->level;
}

void Node2::setL(int l) {
	this->level = l;
}

int Node2::getQ() {
	return this->q_value;
}

void Node2::setQ(int q) {
	this->q_value = q;
}

bool operator< (const Node2 &n1, const Node2 &n2) {
	if (n1.f_value != n2.f_value) {
           return n1.f_value < n2.f_value;
        }

	if (n1.level != n2.level) {
	   return n1.level < n2.level;
	}
	return false;
}
