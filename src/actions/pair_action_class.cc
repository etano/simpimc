#include "pair_action_class.h"

void PairAction::Init(Input &in)
{
  // Read in things
  n_images = in.GetAttribute<int>("n_images");
  n_order = in.GetAttribute<uint>("n_order",0);
  species_a = in.GetAttribute<std::string>("species_a");
  species_list.push_back(species_a);
  species_b = in.GetAttribute<std::string>("species_b");
  species_list.push_back(species_b);
  std::cout << "Setting up pair action between " << species_a << " and " << species_b << "..." << std::endl;
  max_level = in.GetAttribute<uint>("max_level",0);
  use_long_range = in.GetAttribute<bool>("use_long_range",0);
  if (use_long_range) {
    k_cut = in.GetAttribute<double>("k_cut",path.k_c);
    path.SetupKs(k_cut);
  }
  path.GetSpeciesInfo(species_a,species_a_i);
  path.GetSpeciesInfo(species_b,species_b_i);
  if (species_a_i >= 0 && species_b_i >= 0) {
    is_constant = ((species_a_i == species_b_i)
                 && (path.species_list[species_a_i]->n_part == 1
                    || path.species_list[species_a_i]->lambda == 0.));
    is_first_time = true;

    std::string file_name = in.GetAttribute<std::string>("file");
    ReadFile(file_name);

    // Write things to file
    out.Write("Actions/"+name+"/file", file_name);
    out.Write("Actions/"+name+"/n_images", n_images);
    out.Write("Actions/"+name+"/n_order", n_order);
    out.Write("Actions/"+name+"/species_a", species_a);
    out.Write("Actions/"+name+"/species_b", species_b);
    out.Write("Actions/"+name+"/max_level", max_level);
    out.Write("Actions/"+name+"/use_long_range", use_long_range);
    if (use_long_range)
      out.Write("Actions/"+name+"/k_cut", k_cut);
  } else {
    is_constant = true;
    is_first_time = false;
    dUdB_constant = 0.;
    potential_constant = 0.;
  }

}

double PairAction::Potential()
{
  if (is_constant && !is_first_time)
    return potential_constant;
  else {
    double tot = 0.;
    if (species_a_i == species_b_i) {
      for (uint p_i=0; p_i<path.species_list[species_a_i]->n_part-1; ++p_i) {
        for (uint p_j=p_i+1; p_j<path.species_list[species_a_i]->n_part; ++p_j) {
          for (uint b_i=0; b_i<path.n_bead; b_i+=1) {
            uint b_j = b_i + 1;
            vec<double> dr(path.Dr(path(species_a_i,p_i,b_i),path(species_a_i,p_j,b_i)));
            double rMag = mag(dr);
            dr = path.Dr(path(species_a_i,p_i,b_j),path(species_a_i,p_j,b_j));
            double r_p_mag = mag(dr);
            tot += CalcV(rMag,r_p_mag,0);
          }
        }
      }
    } else {
      for (uint p_i=0; p_i<path.species_list[species_a_i]->n_part; ++p_i) {
        for (uint p_j=0; p_j<path.species_list[species_b_i]->n_part; ++p_j) {
          for (uint b_i=0; b_i<path.n_bead; b_i+=1) {
            uint b_j = b_i + 1;
            vec<double> dr(path.Dr(path(species_a_i,p_i,b_i),path(species_b_i,p_j,b_i)));
            double rMag = mag(dr);
            dr = path.Dr(path(species_a_i,p_i,b_j),path(species_b_i,p_j,b_j));
            double r_p_mag = mag(dr);
            tot += CalcV(rMag,r_p_mag,0);
          }
        }
      }
    }

    if (use_long_range)
      tot += CalcVLong();

    if (is_first_time) {
      is_first_time = false;
      potential_constant = tot;
    }

    return tot;
  }
}

double PairAction::DActionDBeta()
{
  if (is_constant && !is_first_time)
    return dUdB_constant;
  else {
    double tot = 0.;
    if (species_a_i == species_b_i) {
      #pragma omp parallel
      {
        for (uint p_i=0; p_i<path.species_list[species_a_i]->n_part-1; ++p_i) {
          #pragma omp for collapse(2) reduction(+:tot) schedule(dynamic)
          for (uint p_j=p_i+1; p_j<path.species_list[species_a_i]->n_part; ++p_j) {
            for (uint b_i=0; b_i<path.n_bead; ++b_i) {
              double rMag, r_p_mag, r_r_p_mag;
              path.DrDrpDrrp(b_i,b_i+1,species_a_i,species_a_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
              double dUdB = CalcdUdBeta(rMag,r_p_mag,r_r_p_mag,0);
              tot += dUdB;
            }
          }
        }
      }
    } else {
      #pragma omp parallel for collapse(3) reduction(+:tot) schedule(dynamic)
      for (uint p_i=0; p_i<path.species_list[species_a_i]->n_part; ++p_i) {
        for (uint p_j=0; p_j<path.species_list[species_b_i]->n_part; ++p_j) {
          for (uint b_i=0; b_i<path.n_bead; ++b_i) {
            double rMag, r_p_mag, r_r_p_mag;
            path.DrDrpDrrp(b_i,b_i+1,species_a_i,species_b_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
            double dUdB = CalcdUdBeta(rMag,r_p_mag,r_r_p_mag,0);
            tot += dUdB;
          }
        }
      }
    }

    if (use_long_range)
      tot += CalcdUdBetaLong();

    if (is_first_time) {
      is_first_time = false;
      dUdB_constant = tot;
    }

    return tot;
  }
}

void PairAction::GenerateParticlePairs(const std::vector<std::pair<uint,uint>> &particles, std::vector<uint> &particles_a, std::vector<uint> &particles_b, std::vector<std::pair<uint,uint>> &particle_pairs)
{
  // Make sure particles are of species A or B and organize them accordingly
  uint n_a(0), n_b(0);
  for (auto& p: particles) {
    uint s_i = p.first;
    uint p_i = p.second;
    if (s_i == species_a_i) {
      particles_a.push_back(p_i);
      n_a++;
    } else if (s_i == species_b_i) {
      particles_b.push_back(p_i);
      n_b++;
    }
  }
  if (n_a==0)
    if ((species_a_i==species_b_i) || (n_b==0))
      return;

  // Make std::vectors of other particles of species A
  std::vector<uint> other_particles_a;
  for (uint p_i=0; p_i<path.species_list[species_a_i]->n_part; ++p_i) {
    if (find(particles_a.begin(), particles_a.end(), p_i)==particles_a.end())
      other_particles_a.push_back(p_i);
  }
  // Make std::vectors of other particles of species B
  std::vector<uint> other_particles_b;
  for (uint p_i=0; p_i<path.species_list[species_b_i]->n_part; ++p_i) {
    if (find(particles_b.begin(), particles_b.end(), p_i)==particles_b.end())
      other_particles_b.push_back(p_i);
  }

  // Homologous
  if (species_a_i == species_b_i) {
    // Loop over A particles with other A particles
    for (auto& p: particles_a)
      for (auto& q: other_particles_a)
        particle_pairs.push_back(std::make_pair(p,q));
    // Loop over A particles with A particles
    for (uint p=0; p<n_a-1; ++p)
      for (uint q=p+1; q<n_a; ++q)
        particle_pairs.push_back(std::make_pair(particles_a[p],particles_a[q]));
  // Heterologous
  } else {
    // Loop over A particles with other B particles
    for (auto& p: particles_a)
      for (auto& q: other_particles_b)
        particle_pairs.push_back(std::make_pair(p,q));
    // Loop other A particles with B particles
    for (auto& p: other_particles_a)
      for (auto& q: particles_b)
        particle_pairs.push_back(std::make_pair(p,q));
    // Loop over A particles with B particles
    for (auto& p: particles_a)
      for (auto& q: particles_b)
        particle_pairs.push_back(std::make_pair(p,q));
  }

}

double PairAction::GetAction(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level)
{
  // Return zero if not relevant
  if (level > max_level || is_constant || species_a_i < 0 || species_b_i < 0)
    return 0.;

  // Generate particle pairs
  std::vector<uint> particles_a, particles_b;
  std::vector<std::pair<uint,uint>> particle_pairs;
  GenerateParticlePairs(particles, particles_a, particles_b, particle_pairs);
  if (particle_pairs.size() == 0)
    return 0.;

  // Sum up contributing terms
  uint skip = 1<<level;
  double tot = 0.;
  for (uint b_i=b0; b_i<b1; b_i+=skip) {
    uint b_j = b_i+skip;
    uint b_k = b_i-skip+path.n_bead;
    for (auto& particle_pair: particle_pairs) {
      double rMag, r_p_mag, r_r_p_mag;
      path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,particle_pair.first,particle_pair.second,rMag,r_p_mag,r_r_p_mag);
      tot += CalcU(rMag,r_p_mag,r_r_p_mag,level);
    }
  }

  // Add in long range part
  if (use_long_range) { // FIXME: currently this assumes level = 0
    if (path.species_list[species_a_i]->need_update_rho_k && path.GetMode()) {
      path.UpdateRhoKP(b0, b1, species_a_i, particles_a, level);
      path.species_list[species_a_i]->need_update_rho_k = false;
    }
    if (path.species_list[species_b_i]->need_update_rho_k && path.GetMode()) {
      path.UpdateRhoKP(b0, b1, species_b_i, particles_b, level);
      path.species_list[species_b_i]->need_update_rho_k = false;
    }
    tot += CalcULong(b0, b1, level);
  }

  return tot;
}

vec<double> PairAction::GetActionGradient(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level)
{
  // Return zero if not relevant
  vec<double> zero_vec;
  zero_vec.zeros(path.n_d);
  if (level > max_level || is_constant || species_a_i < 0 || species_b_i < 0)
    return zero_vec;

  // Generate particle pairs
  std::vector<uint> particles_a, particles_b;
  std::vector<std::pair<uint,uint>> particle_pairs;
  GenerateParticlePairs(particles, particles_a, particles_b, particle_pairs);
  if (particle_pairs.size() == 0)
    return zero_vec;

  // Sum up contributing terms
  uint skip = 1<<level;
  vec<double> tot(zero_vec);
  for (uint b_i=b0; b_i<b1; b_i+=skip) {
    uint b_j = b_i+skip;
    uint b_k = b_i-skip+path.n_bead;
    for (auto& particle_pair: particle_pairs) {
      tot += CalcGradientU(b_i,b_j,particle_pair.first,particle_pair.second,level);
      tot += CalcGradientU(b_i,b_k,particle_pair.first,particle_pair.second,level);
    }
  }

  // FIXME: Ignoring long range part for now
  //// Add in long range part
  //if (use_long_range) { // fixme: currently this assumes level = 0
  //  if (path.species_list[species_a_i]->need_update_rho_k && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, species_a_i, particles_a, level);
  //    path.species_list[species_a_i]->need_update_rho_k = false;
  //  }
  //  if (path.species_list[species_b_i]->need_update_rho_k && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, species_b_i, particles_b, level);
  //    path.species_list[species_b_i]->need_update_rho_k = false;
  //  }
  //  tot += CalcGradientULong(b0, b1, level);
  //}

  return tot;
}

double PairAction::GetActionLaplacian(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level)
{
  // Return zero if not relevant
  if (level > max_level || is_constant || species_a_i < 0 || species_b_i < 0)
    return 0.;

  // Generate particle pairs
  std::vector<uint> particles_a, particles_b;
  std::vector<std::pair<uint,uint>> particle_pairs;
  GenerateParticlePairs(particles, particles_a, particles_b, particle_pairs);
  if (particle_pairs.size() == 0)
    return 0.;

  // Sum up contributing terms
  uint skip = 1<<level;
  double tot = 0.;
  for (uint b_i=b0; b_i<b1; b_i+=skip) {
    uint b_j = b_i+skip;
    uint b_k = b_i-skip+path.n_bead;
    for (auto& particle_pair: particle_pairs) {
      tot += CalcLaplacianU(b_i,b_j,particle_pair.first,particle_pair.second,level)
          + CalcLaplacianU(b_i,b_k,particle_pair.first,particle_pair.second,level);
    }
  }

  // FIXME: Ignoring long range part for now
  //// Add in long range part
  //if (use_long_range) { // fixme: currently this assumes level = 0
  //  if (path.species_list[species_a_i]->need_update_rho_k && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, species_a_i, particles_a, level);
  //    path.species_list[species_a_i]->need_update_rho_k = false;
  //  }
  //  if (path.species_list[species_b_i]->need_update_rho_k && path.GetMode()) {
  //    path.UpdateRhoKP(b0, b1, species_b_i, particles_b, level);
  //    path.species_list[species_b_i]->need_update_rho_k = false;
  //  }
  //  tot += CalcLaplacianULong(b0, b1, level);
  //}

  return tot;
}

vec<double> PairAction::CalcGradientU(const uint b_i, const uint b_j, const uint p_i, const uint p_j, const uint level)
{
  // Store original position for particle i
  std::shared_ptr<Bead> b(path(species_a_i,p_i,b_i));
  vec<double> r0(b->r);

  // Numerical tolerance
  double eps = 1.e-7;

  // Calculate numerical gradient
  double rMag, r_p_mag, r_r_p_mag;
  vec<double> tot(path.n_d);
  for (uint d_i=0; d_i<path.n_d; ++d_i) {
    b->r(d_i) = r0(d_i) + eps;
    path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
    double f1 = CalcU(rMag,r_p_mag,r_r_p_mag,level);
    b->r(d_i) = r0(d_i) - eps;
    path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
    double f2 = CalcU(rMag,r_p_mag,r_r_p_mag,level);
    tot(d_i) = (f1-f2)/(2.*eps);
    b->r(d_i) = r0(d_i);
  }

  return tot;
}

double PairAction::CalcLaplacianU(const uint b_i, const uint b_j, const uint p_i, const uint p_j, const uint level)
{
  // Store original position for particle i
  std::shared_ptr<Bead> b(path(species_a_i,p_i,b_i));
  vec<double> r0(b->r);

  // Numerical tolerance
  double eps = 1.e-7;

  // Calculate numerical gradient
  double rMag, r_p_mag, r_r_p_mag;
  double tot = 0.;
  for (uint d_i=0; d_i<path.n_d; ++d_i) {
    path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
    double f0 = CalcU(rMag,r_p_mag,r_r_p_mag,level);
    b->r(d_i) = r0(d_i) + eps;
    path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
    double f1 = CalcU(rMag,r_p_mag,r_r_p_mag,level);
    b->r(d_i) = r0(d_i) - eps;
    path.DrDrpDrrp(b_i,b_j,species_a_i,species_b_i,p_i,p_j,rMag,r_p_mag,r_r_p_mag);
    double f2 = CalcU(rMag,r_p_mag,r_r_p_mag,level);
    tot += (f1+f2-2.*f0)/(eps*eps);
    b->r(d_i) = r0(d_i);
  }

  return tot;
}

double PairAction::ImportanceWeight()
{
  return is_importance_weight ? exp(DActionDBeta()/path.n_bead) : 1.;
}

void PairAction::Write() {}

void PairAction::Accept() // fixme: will accept even when unaffected
{
  if (use_long_range && species_a_i >= 0 && species_b_i >= 0) {
    path.species_list[species_a_i]->need_update_rho_k = true;
    path.species_list[species_b_i]->need_update_rho_k = true;
  }

}

void PairAction::Reject() // fixme: will accept even when unaffected, why is this true?
{
  if (use_long_range && species_a_i >= 0 && species_b_i >= 0) {
    path.species_list[species_a_i]->need_update_rho_k = true;
    path.species_list[species_b_i]->need_update_rho_k = true;
  }

}
