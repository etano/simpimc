#ifndef SIMPIMC_ACTIONS_H_
#define SIMPIMC_ACTIONS_H_

#include "pair_action/bare_pair_action_class.h"
#include "pair_action/david_pair_action_class.h"
#include "pair_action/ilkka_pair_action_class.h"
#include "single_action/kinetic_class.h"
#include "single_action/nodal/free_nodal_class.h"
#include "single_action/nodal/optimized_free_nodal_class.h"
#include "single_action/nodal/optimized_sho_nodal_class.h"
#include "single_action/trap_class.h"

std::shared_ptr<Action> ActionFactory(Input &in, IO &out, Path &path) {
    std::string type = in.GetAttribute<std::string>("type");
    if (type == "Kinetic")
        return std::make_shared<Kinetic>(path, in, out);
    else if (type == "HarmonicTrap")
        return std::make_shared<Trap>(path, in, out);
    else if (type == "FreeNodal")
        return std::make_shared<FreeNodal>(path, in, out);
    else if (type == "OptimizedFreeNodal")
        return std::make_shared<OptimizedFreeNodal>(path, in, out);
    else if (type == "OptimizedSHONodal")
        return std::make_shared<OptimizedSHONodal>(path, in, out);
    else if (type == "BarePairAction")
        return std::make_shared<BarePairAction>(path, in, out);
    else if (type == "DavidPairAction")
        return std::make_shared<DavidPairAction>(path, in, out);
    else if (type == "IlkkaPairAction")
        return std::make_shared<IlkkaPairAction>(path, in, out);
    else {
        std::cerr << "ERROR: Unrecognized Action, " << type << std::endl;
        exit(1);
    }
}

#endif
