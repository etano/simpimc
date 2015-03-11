#include "pair_correlation_class.h"

void PairCorrelation::Init(Input &in)
{
  // Read in species info
  species_A = in.GetAttribute<std::string>("species_A");
  species_B = in.GetAttribute<std::string>("species_B");
  path.GetSpeciesInfo(species_A, species_A_i);
  path.GetSpeciesInfo(species_B, species_B_i);

  // Read in grid info
  double r_min = in.GetAttribute<double>("r_min",0.);
  double r_max = in.GetAttribute<double>("r_max",path.L/2.);
  uint n_r = in.GetAttribute<double>("n_r",1000);
  gr.x.CreateGrid(r_min,r_max,n_r);
  gr.y.zeros(n_r);

  // Compute rs
  vec<double> rs(gr.x.n_r-1);
  for (uint i=0; i<gr.x.n_r-1; i++) {
    double r1 = gr.x(i);
    double r2 = gr.x(i+1);
    if (path.n_d == 3)
      rs(i) = 0.75 * (r2*r2*r2*r2-r1*r1*r1*r1)/(r2*r2*r2-r1*r1*r1);
    else if (path.n_d == 2)
      rs(i) = (r2*r2*r2-r1*r1*r1)/(r2*r2-r1*r1); // fixme: Not sure if 2D and 1D are correct here
    else if (path.n_d == 1)
      rs(i) = 0.5*(r2-r1);
  }

  // Write things to file
  out.Write(prefix+"/species_A", species_A);
  out.Write(prefix+"/species_B", species_B);
  out.Write(prefix+"/r_min", r_min);
  out.Write(prefix+"/r_max", r_max);
  out.Write(prefix+"/n_r", n_r);
  out.Write(prefix+"/x", rs);

  Reset();
}

void PairCorrelation::Reset()
{
  n_measure = 0;
  gr.y.zeros();
}

void PairCorrelation::Accumulate()
{
  path.SetMode(1);
  // Homogeneous
  if (species_A_i == species_B_i) {
    for (uint b_i=0; b_i<path.n_bead; ++b_i) {
      for (uint p_i=0; p_i<path.species_list[species_A_i]->n_part-1; ++p_i) {
        for (uint jP=p_i+1; jP<path.species_list[species_A_i]->n_part; ++jP) {
          vec<double> dr(path.Dr(path(species_A_i,p_i,b_i),path(species_A_i,jP,b_i)));
          uint i = gr.x.ReverseMap(mag(dr));
          if (i < gr.x.n_r)
            gr.y(i) = gr.y(i) + 1.*path.sign*path.importance_weight;
        }
      }
    }
  // Homologous
  } else {
    for (uint b_i=0; b_i<path.n_bead; ++b_i) {
      for (uint p_i=0; p_i<path.species_list[species_A_i]->n_part; ++p_i) {
        for (uint jP=0; jP<path.species_list[species_B_i]->n_part; ++jP) {
          vec<double> dr(path.Dr(path(species_A_i,p_i,b_i),path(species_B_i,jP,b_i)));
          uint i = gr.x.ReverseMap(mag(dr));
          if (i < gr.x.n_r)
            gr.y(i) = gr.y(i) + 1.*path.sign*path.importance_weight;
        }
      }
    }
  }

  n_measure += 1;
}

void PairCorrelation::Write()
{
  if (n_measure > 0) {
    // Normalize histograms
    uint N_A = path.species_list[species_A_i]->n_part;
    uint N_B = path.species_list[species_B_i]->n_part;
    double vol = path.vol;
    double norm;
    if (species_A_i == species_B_i)
      norm = 0.5*n_measure*N_A*(N_A-1.)*path.n_bead/vol;
    else
      norm = n_measure*N_A*N_B*path.n_bead/vol;
    for (uint i=0; i<gr.x.n_r; i++) {
      double r1 = gr.x(i);
      double r2 = (i<(gr.x.n_r-1)) ? gr.x(i+1):(2.*gr.x(i)-gr.x(i-1));
      double r = 0.5*(r1+r2);
      double bin_vol;
      if (path.n_d == 3)
        bin_vol = 4.*M_PI/3. * (r2*r2*r2-r1*r1*r1);
      else if (path.n_d == 2)
        bin_vol = M_PI * (r2*r2-r1*r1);
      else if (path.n_d == 1)
        bin_vol = r2-r1;
      gr.y(i) = gr.y(i)/(bin_vol*norm);
      //gr.y(i) = gr.y(i)/(norm);
    }

    // Write to file
    if (first_time) {
      first_time = 0;
      out.CreateExtendableDataSet("/"+prefix, "y", gr.y);
    } else {
      out.AppendDataSet("/"+prefix, "y", gr.y);
    }

    Reset();
  }
}
