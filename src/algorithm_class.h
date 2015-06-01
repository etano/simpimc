#ifndef SIMPIMC_ALGORITHM_CLASS_H_
#define SIMPIMC_ALGORITHM_CLASS_H_

#include "loop_class.h"
#include "write_class.h"
#include "actions/kinetic_class.h"
#include "actions/free_nodal_class.h"
#include "actions/optimized_free_nodal_class.h"
#include "actions/optimized_sho_nodal_class.h"
#include "actions/bare_pair_action_class.h"
#include "actions/david_pair_action_class.h"
#include "actions/ilkka_pair_action_class.h"
#include "actions/trap_class.h"
#include "moves/bisect_class.h"
#include "moves/displace_particle_class.h"
#include "moves/perm_bisect_table_class.h"
#include "moves/perm_bisect_iterative_class.h"
#include "moves/shift_ref_slice_class.h"
#include "moves/vary_action_class.h"
#include "observables/contact_density_class.h"
#include "observables/energy_class.h"
#include "observables/importance_weight_class.h"
#include "observables/pair_correlation_class.h"
#include "observables/path_dump_class.h"
#include "observables/permutation_class.h"
#include "observables/record_optimized_action_class.h"
#include "observables/sign_class.h"
#include "observables/structure_factor_class.h"
#include "observables/time_class.h"

/// Class the actually holds all the events and data objects. This includes the path object which is the main container of all relevant information to the path.
class Algorithm
{
private:
  Loop main_loop; ///< Main loop
  std::vector<std::shared_ptr<Event>> events; /// < Vector of events
  std::vector<std::shared_ptr<Action>> actions; ///< Vector of actions
  Path path; ///< Main datastructure

  /// Initializes all events and objects
  void Init(Input &in, IO &out, RNG &rng, const uint32_t proc_i);
public:
  /// Constructor calls Init.
  Algorithm(Input &in, IO &out, RNG &rng, const uint32_t proc_i)
    : path(proc_i)
  {
    Init(in, out, rng, proc_i);
  }

  /// Runs the algorithm by calling the DoEvent of the main_loop
  void Run() { main_loop.DoEvent(); };
};

#endif // SIMPIMC_ALGORITHM_CLASS_H_
