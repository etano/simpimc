#include "path_dump_class.h"

void PathDump::Init(Input &in)
{
  n_dump = 0;
  Reset();
}

void PathDump::Reset()
{
  n_write_calls = 0;
}

void PathDump::Accumulate() {}

void PathDump::Write()
{
  if (!(n_write_calls % skip)) {
    n_dump += 1;

    // Loop through species
    for (uint32_t s_i=0; s_i<path.n_species; ++s_i) {
       // Get species info
      path.GetSpeciesInfo(path.species_list[s_i]->name,s_i);
      uint32_t n_part = path.species_list[s_i]->n_part;

      // Get positions
      cube<double> path_positions(path.n_d,path.n_bead,n_part);
      for (uint32_t p_i=0; p_i<n_part; ++p_i)
        for (uint32_t b_i=0; b_i<path.n_bead; ++b_i)
          for (uint32_t d_i=0; d_i<path.n_d; ++d_i)
            path_positions(d_i,b_i,p_i) = path(s_i,p_i,b_i)->r(d_i);

      // Get permutation
      mat<double> path_permutation(2,n_part);
      for (uint32_t p_i=0; p_i<n_part; ++p_i) {
        path_permutation(0,p_i) = path(s_i,p_i,0)->prev->p;
        path_permutation(1,p_i) = path(s_i,p_i,path.n_bead-1)->next->p;
      }

      // Write to file
      if (first_time) {
        out.CreateGroup(prefix+path.species_list[s_i]->name+"/");
        out.Write(prefix+path.species_list[s_i]->name+"/n_dump", n_dump);
        out.CreateExtendableDataSet(prefix+path.species_list[s_i]->name+"/", "positions", path_positions);
        out.CreateExtendableDataSet(prefix+path.species_list[s_i]->name+"/", "permutation", path_permutation);
      } else {
        out.Rewrite(prefix+path.species_list[s_i]->name+"/n_dump", n_dump);
        out.AppendDataSet(prefix+path.species_list[s_i]->name+"/", "positions", path_positions);
        out.AppendDataSet(prefix+path.species_list[s_i]->name+"/", "permutation", path_permutation);
      }

    }

    if (first_time)
      first_time = 0;

    Reset();
  }

  n_write_calls += 1;
}
