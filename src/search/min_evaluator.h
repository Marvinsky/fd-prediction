#ifndef MIN_EVALUATOR_H
#define MIN_EVALUATOR_H

#include "combining_evaluator.h"

#include <vector>

class Options;

class MinEvaluator : public CombiningEvaluator {
protected:
    virtual int combine_values(const std::vector<int> &values);
public:
    MinEvaluator(const Options &opts);
    MinEvaluator(const std::vector<ScalarEvaluator *> &subevaluators);
    ~MinEvaluator();
};

#endif
