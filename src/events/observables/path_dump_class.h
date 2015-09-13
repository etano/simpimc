#ifndef SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_
#define SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_

#include "observable_class.h"

/// Dump out all path information
class PathDump : public Observable
{
private:
  uint32_t n_dump; ///< Number of times dumped
  uint32_t n_write_calls; ///< Number of times write has been called

  /// Accumulate the observable
  virtual void Accumulate() {};

  /// Reset the observable's counters
  virtual void Reset()
  {
    n_write_calls = 0;
  }
public:
  /// Constructor calls Init
  PathDump(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    n_dump = 0;
    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (!(n_write_calls % skip)) {
      n_dump += 1;

      // Loop through species
      for (const auto& species : path.GetSpecies()) {
        // Get positions
        cube<double> path_positions(path.GetND(),species->GetNBead(),species->GetNPart());
        for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i)
          for (uint32_t b_i=0; b_i<species->GetNBead(); ++b_i)
            for (uint32_t d_i=0; d_i<path.GetND(); ++d_i)
              path_positions(d_i,b_i,p_i) = species->GetBead(p_i,b_i)->GetR()(d_i);

        // Get permutation
        mat<double> path_permutation(2,species->GetNPart());
        for (uint32_t p_i=0; p_i<species->GetNPart(); ++p_i) {
          path_permutation(0,p_i) = species->GetBead(p_i,0)->GetPrevBead(1)->GetP();
          path_permutation(1,p_i) = species->GetBead(p_i,species->GetNBead()-1)->GetNextBead(1)->GetP();
        }

        // Write to file
        if (first_time) {
          out.CreateGroup(prefix+species->GetName()+"/");
          out.Write(prefix+species->GetName()+"/n_dump", n_dump);
          out.CreateExtendableDataSet(prefix+species->GetName()+"/", "positions", path_positions);
          out.CreateExtendableDataSet(prefix+species->GetName()+"/", "permutation", path_permutation);
        } else {
          out.Rewrite(prefix+species->GetName()+"/n_dump", n_dump);
          out.AppendDataSet(prefix+species->GetName()+"/", "positions", path_positions);
          out.AppendDataSet(prefix+species->GetName()+"/", "permutation", path_permutation);
        }

      }

      if (first_time)
        first_time = 0;

      Reset();
    }

    n_write_calls += 1;
  }
};

#endif // SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_
