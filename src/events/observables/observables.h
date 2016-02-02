#ifndef SIMPIMC_OBSERVABLES_H_
#define SIMPIMC_OBSERVABLES_H_

#include "contact_density_class.h"
#include "contact_density_naive_class.h"
#include "energy_class.h"
#include "importance_weight_class.h"
#include "pair_correlation_class.h"
#include "path_dump_class.h"
#include "permutation_class.h"
#include "record_optimized_action_class.h"
#include "sign_class.h"
#include "structure_factor_class.h"
#include "time_class.h"

std::shared_ptr<Observable> ObservableFactory(Input &in, IO &out, Path& path, std::vector<std::shared_ptr<Action>>& actions, std::vector<std::shared_ptr<Event>>& events)
{
  std::string type = in.GetAttribute<std::string>("type");
  if (type == "ContactDensity")
    return std::make_shared<ContactDensity>(path,actions,in,out);
  else if (type == "ContactDensityNaive")
    return std::make_shared<ContactDensityNaive>(path,actions,in,out);
  else if (type == "Energy")
    return std::make_shared<Energy>(path,actions,in,out);
  else if (type == "ImportanceWeight")
    return std::make_shared<ImportanceWeight>(path,actions,in,out);
  else if (type == "PairCorrelation")
    return std::make_shared<PairCorrelation>(path,in,out);
  else if (type == "PathDump")
    return std::make_shared<PathDump>(path,in,out);
  else if (type == "Permutation")
    return std::make_shared<Permutation>(path,in,out);
  else if (type == "RecordOptimizedAction")
    return std::make_shared<RecordOptimizedAction>(path,actions,in,out);
  else if (type == "Sign")
    return std::make_shared<Sign>(path,in,out);
  else if (type == "StructureFactor")
    return std::make_shared<StructureFactor>(path,in,out);
  else if (type == "Time")
    return std::make_shared<Time>(path,events,in,out);
  else {
    std::cerr << "ERROR: Unrecognized observable, " << type << std::endl;
    exit(1);
  }
}

#endif
