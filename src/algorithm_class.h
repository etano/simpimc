#ifndef SIMPIMC_ALGORITHM_CLASS_H_
#define SIMPIMC_ALGORITHM_CLASS_H_

#include "loop_class.h"
#include "write_class.h"
#include "actions/kinetic_class.h"
#include "actions/free_nodal_class.h"
#include "actions/optimized_nodal_class.h"
#include "actions/bare_pair_action_class.h"
#include "actions/david_pair_action_class.h"
#include "actions/ilkka_pair_action_class.h"
#include "actions/importance_pair_action_class.h"
#include "actions/trap_class.h"
#include "moves/bisect_class.h"
#include "moves/displace_particle_class.h"
#include "moves/perm_bisect_class.h"
#include "moves/perm_bisect_iterative_class.h"
#include "moves/shift_ref_slice_class.h"
#include "moves/vary_optimized_nodal_class.h"
#include "observables/contact_density_class.h"
#include "observables/energy_class.h"
#include "observables/importance_weight_class.h"
#include "observables/pair_correlation_class.h"
#include "observables/path_dump_class.h"
#include "observables/permutation_class.h"
#include "observables/record_optimized_nodal_class.h"
#include "observables/sign_class.h"
#include "observables/structure_factor_class.h"
#include "observables/time_class.h"

class Algorithm
{
public:
  // Constructor
  Algorithm(Communicator &world_comm, Communicator &inter_comm, Communicator &intra_comm, Input &in, IO &out, RNG &rng)
   : path(world_comm,inter_comm,intra_comm)
  {
    Init(in, out, rng, inter_comm);
  }

  void Init(Input &in, IO &out, RNG &rng, Communicator &inter_comm);
  void Run() { main_loop.DoEvent(); };

  // Events
  std::vector<std::shared_ptr<Event>> events;

  // Main loop
  Loop main_loop;

  // Actions
  std::vector<std::shared_ptr<Action>> actions;

  // Datastructure
  Path path;
};

#endif // SIMPIMC_ALGORITHM_CLASS_H_
