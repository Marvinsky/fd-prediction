#include "dijkstra_search_heuristic.h"
#include "globals.h"
#include "global_operator.h"
#include "global_state.h"

#include "option_parser.h"
#include "plugin.h"

#include "limits"
#include <utility>


using namespace std;

DijkstraSearchHeuristic::DijkstraSearchHeuristic(const Options &opts)
	: Heuristic(opts) {
 	min_operator_cost = 0;//numeric_limits<int>::max();
 	for (size_t i = 0; i < g_operators.size(); ++i) {
 		min_operator_cost = min(min_operator_cost, get_adjusted_cost(g_operators[i]));
  	}
}

DijkstraSearchHeuristic::~DijkstraSearchHeuristic() {
}

void DijkstraSearchHeuristic::initialize() {
  cout << "Initializing DijkstraSearchHeuristic ..." << endl;
}

int DijkstraSearchHeuristic::compute_heuristic(const GlobalState &state) {
  if (test_goal(state)) {
    return 0;
  } else {
    return min_operator_cost;
  }
}

static ScalarEvaluator *_parse(OptionParser &parser) {
  Heuristic::add_options_to_parser(parser);
  Options opts = parser.parse();
  if (parser.dry_run())
    return 0;
  else
    return new DijkstraSearchHeuristic(opts);
}


static Plugin<ScalarEvaluator> _plugin("dijkstra", _parse);
