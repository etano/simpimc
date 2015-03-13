#include "structure_factor_class.h"

void StructureFactor::Init(Input &in)
{
  // Read in species info
  species_a = in.GetAttribute<std::string>("species_a");
  species_b = in.GetAttribute<std::string>("species_b");
  path.GetSpeciesInfo(species_a, species_a_i);
  path.GetSpeciesInfo(species_b, species_b_i);
  k_cut = in.GetAttribute<double>("k_cut", path.k_c);

  // Resize
  path.SetupKs(k_cut);
  sk.zeros(path.ks.size());

  // Write things to file
  out.Write(prefix+"/species_a", species_a);
  out.Write(prefix+"/species_b", species_b);
  out.Write(prefix+"/k_cut", k_cut);
  out.CreateExtendableDataSet("/"+prefix+"/", "x", path.mag_ks[0]);
  for (uint32_t k_i=1; k_i<path.mag_ks.size(); ++k_i)
    out.AppendDataSet("/"+prefix+"/", "x", path.mag_ks[k_i]);

  Reset();
}

void StructureFactor::Reset()
{
  n_measure = 0;
  sk.zeros();
}

void StructureFactor::Accumulate()
{
  path.SetMode(1);
  for (uint32_t k_i=0; k_i<path.k_indices.size(); k_i++) {
    if (path.mag_ks[k_i] < k_cut) {
      for (uint32_t b_i=0; b_i<path.n_bead; ++b_i) {
        sk(k_i) += path.sign*path.importance_weight*CMag2(path.rho_k(b_i,species_a_i)(k_i),path.rho_k(b_i,species_b_i)(k_i));
      }
    }
  }

  //if (species_b_i != species_a_i)
  //  sk = 2.*sk;

  n_measure += 1;
}

void StructureFactor::Write()
{
  if (n_measure > 0) {
    // Normalize histograms
    uint32_t N_a = path.species_list[species_a_i]->n_part;
    uint32_t N_b = path.species_list[species_b_i]->n_part;
    double norm = n_measure*path.n_bead*N_a*N_b;
    sk = sk/norm;

    // Write to file
    if (first_time) {
      first_time = 0;
      out.CreateExtendableDataSet("/"+prefix+"/", "y", sk);
    } else {
      out.AppendDataSet("/"+prefix+"/", "y", sk);
    }

    Reset();
  }
}
