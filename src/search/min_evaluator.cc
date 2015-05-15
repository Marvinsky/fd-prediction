#include "min_evaluator.h"

#include "option_parser.h"
#include "plugin.h"

#include <cassert>

using namespace std;


MinEvaluator::MinEvaluator(const Options &opts)
    : CombiningEvaluator(opts.get_list<ScalarEvaluator *>("evals")) {
}
MinEvaluator::MinEvaluator(const vector<ScalarEvaluator *> &evals)
    : CombiningEvaluator(evals) {
}

MinEvaluator::~MinEvaluator() {
}

int MinEvaluator::combine_values(const vector<int> &values) {
    int result = 0;
    for (size_t i = 0; i < values.size(); ++i) {
      cout<<"\tMining_h["<<i<<"]:"<<values[i]<<endl;
        assert(values[i] >= 0);
        result = min(result, values[i]);
    }
    //exit(0);
    return result;
}

static ScalarEvaluator *_parse(OptionParser &parser) {
  cout<<"calling ScalarEvaluator Min parser"<<endl;
    parser.add_list_option<ScalarEvaluator *>("evals");
    Options opts = parser.parse();

    opts.verify_list_non_empty<ScalarEvaluator *>("evals");

    if (parser.dry_run())
        return 0;
    else{
      cout<<"returning MinEvaluator"<<endl;fflush(stdout);
        return new MinEvaluator(opts);
    }
}

static Plugin<ScalarEvaluator> _plugin("min", _parse);

/*commented out to silence compiler warning while this is unused.

static ScalarEvaluator *create(const vector<string> &config,
                               int start, int &end, bool dry_run) {
    if (config[start + 1] != "(")
        throw ParseError(start + 1);

    // create evaluators
    vector<ScalarEvaluator *> evals;
    OptionParser::instance()->parse_scalar_evaluator_list(
        config, start + 2, end, false, evals, dry_run);

    if (evals.empty())
        throw ParseError(end);
    // need at least one evaluator

    end++;
    if (config[end] != ")")
        throw ParseError(end);

    if (dry_run)
        return 0;
    else
        return new MinEvaluator(evals);
}*/


// TODO: Comment this in once the max/sum evaluator stuff is fixed.
//       For now, it's commented out to use the IPC implementation of
//       max again, see issue181.
// static Plugin<ScalarEvaluator> plugin("max", create);
