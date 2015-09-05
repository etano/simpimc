#ifndef SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_
#define SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_

#include "observable_class.h"

/// Record the pair correlation function between two species of particles
class PairCorrelation : public Observable
{
private:
  uint32_t species_a_i; ///< Index of first species
  uint32_t species_b_i; ///< Index of second species
  Histogram gr; ///< Histogram representing g(r)
  std::string species_a; ///< Name of first species
  std::string species_b; ///< Name of second species

  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);
    // Homogeneous
    int sign = path.CalcSign();
    if (species_a_i == species_b_i) {
      for (uint32_t b_i=0; b_i<path.n_bead; ++b_i) {
        for (uint32_t p_i=0; p_i<path.species_list[species_a_i]->n_part-1; ++p_i) {
          for (uint32_t p_j=p_i+1; p_j<path.species_list[species_a_i]->n_part; ++p_j) {
            vec<double> dr(path.Dr(path(species_a_i,p_i,b_i),path(species_a_i,p_j,b_i)));
            uint32_t i = gr.x.ReverseMap(mag(dr));
            if (i < gr.x.n_r)
              gr.y(i) = gr.y(i) + 1.*sign*path.importance_weight;
          }
        }
      }
    // Homologous
    } else {
      for (uint32_t b_i=0; b_i<path.n_bead; ++b_i) {
        for (uint32_t p_i=0; p_i<path.species_list[species_a_i]->n_part; ++p_i) {
          for (uint32_t p_j=0; p_j<path.species_list[species_b_i]->n_part; ++p_j) {
            vec<double> dr(path.Dr(path(species_a_i,p_i,b_i),path(species_b_i,p_j,b_i)));
            uint32_t i = gr.x.ReverseMap(mag(dr));
            if (i < gr.x.n_r)
              gr.y(i) = gr.y(i) + 1.*sign*path.importance_weight;
          }
        }
      }
    }

    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    n_measure = 0;
    gr.y.zeros();
  }

public:
  /// Constructor calls Init and sets output data_type
  PairCorrelation(Path &path, Input &in, IO &out)
    : Observable(path, in, out, "histogram")
  {
    // Read in species info
    species_a = in.GetAttribute<std::string>("species_a");
    species_b = in.GetAttribute<std::string>("species_b");
    path.GetSpeciesInfo(species_a, species_a_i);
    path.GetSpeciesInfo(species_b, species_b_i);

    // Read in grid info
    double r_min = in.GetAttribute<double>("r_min",0.);
    double r_max = in.GetAttribute<double>("r_max",path.L/2.);
    uint32_t n_r = in.GetAttribute<double>("n_r",1000);
    gr.x.CreateGrid(r_min,r_max,n_r);
    gr.y.zeros(n_r);

    // Compute rs
    vec<double> rs(gr.x.n_r-1);
    for (uint32_t i=0; i<gr.x.n_r-1; i++) {
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
    out.Write(prefix+"/species_a", species_a);
    out.Write(prefix+"/species_b", species_b);
    out.Write(prefix+"/r_min", r_min);
    out.Write(prefix+"/r_max", r_max);
    out.Write(prefix+"/n_r", n_r);
    out.Write(prefix+"/x", rs);

    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {
      // Normalize histograms
      uint32_t N_a = path.species_list[species_a_i]->n_part;
      uint32_t N_b = path.species_list[species_b_i]->n_part;
      double vol = path.vol;
      double norm;
      if (species_a_i == species_b_i)
        norm = 0.5*n_measure*N_a*(N_a-1.)*path.n_bead/vol;
      else
        norm = n_measure*N_a*N_b*path.n_bead/vol;
      for (uint32_t i=0; i<gr.x.n_r; i++) {
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

};

#endif // SIMPIMC_OBSERVABLES_PAIR_CORRELATION_CLASS_H_
