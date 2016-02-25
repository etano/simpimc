#ifndef SIMPIMC_MOVES_H_
#define SIMPIMC_MOVES_H_

#include "single_species_move/displace_particle_class.h"
#include "single_species_move/open_class.h"
#include "single_species_move/shift_ref_slice_class.h"
#include "single_species_move/bisect/bisect_class.h"
#include "single_species_move/bisect/perm_bisect/perm_bisect_class.h"
#include "single_species_move/bisect/perm_bisect/perm_bisect_iterative_class.h"
#include "single_species_move/bisect/perm_bisect/perm_bisect_table_class.h"
#include "vary_action_class.h"

std::shared_ptr<Move> MoveFactory(Input &in, IO &out, Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &actions) {
    std::string type = in.GetAttribute<std::string>("type");
    if (type == "Bisect")
        return std::make_shared<Bisect>(path, rng, actions, in, out);
    else if (type == "DisplaceParticle")
        return std::make_shared<DisplaceParticle>(path, rng, actions, in, out);
    else if (type == "Open")
        return std::make_shared<Open>(path, rng, actions, in, out);
    else if (type == "PermBisectIterative")
        return std::make_shared<PermBisectIterative>(path, rng, actions, in, out);
    else if (type == "PermBisectTable")
        return std::make_shared<PermBisectTable>(path, rng, actions, in, out);
    else if (type == "ShiftRefSlice")
        return std::make_shared<ShiftRefSlice>(path, rng, actions, in, out);
    else if (type == "VaryAction")
        return std::make_shared<VaryAction>(path, rng, actions, in, out);
    else {
        std::cerr << "ERROR: Unrecognized Move, " << type << std::endl;
        exit(1);
    }
}

#endif
