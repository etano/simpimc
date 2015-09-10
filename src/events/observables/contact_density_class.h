#ifndef SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
#define SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_

#include "observable_class.h"

namespace Contact_Density_Optimization_Functions {
    double f_simple(vec<double> ri, vec<double> RA){
        return 1;
    }
    
    vec<double> gradient_f_simple(vec<double> ri, vec<double> RA){
        return zeros<vec<double>>(path.GetND());
    }
    double laplace_f_simple(vec<double> ri, vec<double> RA){
        return 0;
    }
}

/// Measures the contact density between two species of particles. Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
class ContactDensity : public Observable
{
private:
  double total; ///< Running total
  uint32_t z_a; ///< Charge of ion-like particle
  uint32_t species_a_i; ///< Ion species index
  uint32_t species_b_i; ///< Other species index
  std::string species_a; ///< Name of ion species
  std::string species_b; ///< Name of other species
  std::vector<std::shared_ptr<Action>> action_list; ///< Vector of pointers to actions that involves species_a or species_b
  std::vector<std::shared_ptr<Action>> &full_action_list; ///< Vector of pointers to all actions
  double (*Function_f)(vec<double> ri, vec<double> RA); ///< Function pointer for the possible generalization
  vec<double> (*Function_gradient_f)(vec<double> ri, vec<double> RA)
  double (*Function_laplace_f)(vec<double> ri, vec<double> RA)
  
  vec<double> getRelevantNormalVector(vec<double> r1,vec<double> r2){
    vec<double> n(zeros<vec<double>>(path.GetND())); 
    double Compare_Measure=path.GetL()/2.0;
    for(int d=0;d<path.GetND();++d){
        if(r1[d]>0.9*path.GetL()&&r2[d]<0.1*path.GetL()) {//TODO check with etano if this is hacking or not
            n[d]=1;
        }
        else if(r1[d]<0.1*path.GetL()&&r2[d]>0.9*path.GetL()){
            n[d]=-1;
        }
    }
    return n;
  }


  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);

    // Form particle pairs
    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> particle_pairs;
    if (species_a_i == species_b_i) { // Homogeneous
      for (uint32_t p_i=0; p_i<path.species_list[species_a_i]->n_part-1; ++p_i) {
        for (uint32_t p_j=p_i+1; p_j<path.species_list[species_b_i]->n_part; ++p_j) {
          std::vector<std::pair<uint32_t,uint32_t>> particles;
          particles.push_back(std::make_pair(species_a_i,p_i));
          particles.push_back(std::make_pair(species_b_i,p_j));
          particle_pairs.push_back(particles);
        }
      }
    } else { // Homologous
      for (uint32_t p_i=0; p_i<path.species_list[species_a_i]->n_part; ++p_i) {
        for (uint32_t p_j=0; p_j<path.species_list[species_b_i]->n_part; ++p_j) {
          std::vector<std::pair<uint32_t,uint32_t>> particles;
          particles.push_back(std::make_pair(species_a_i,p_i));
          particles.push_back(std::make_pair(species_b_i,p_j));
          particle_pairs.push_back(particles);
        }
      }
    }

    // Add up contact probability
    // TODO: Currently only looking at origin
    double tot = 0.;
    size_t n_particle_pairs(particle_pairs.size());
    vec<bool> Checked(zeros<vec<bool>>(n_part_tot));//Make sure to only once measure the boundary terms
    #pragma omp parallel for collapse(2) reduction(+:tot)
    for (uint32_t pp_i=0; pp_i<n_particle_pairs; ++pp_i) {
      for (uint32_t b_i=0; b_i<path.GetNBead(); ++b_i) {

        // Set r's
        vec<double> RA = path(particle_pairs[pp_i][0].first,particle_pairs[pp_i][0].second,b_i)->r;
        vec<double> ri = path(particle_pairs[pp_i][1].first,particle_pairs[pp_i][1].second,b_i)->r;
        vec<double> ri_nextBead = path(particle_pairs[pp_i][1].first,particle_pairs[pp_i][1].second,b_i+1)->r;
        // Get differences
        vec<double> ri_RA(path.Dr(ri, RA));
        double mag_ri_RA = mag(ri_RA);
        double mag_Delta_ri = mag(path.Dr(ri,ri_nextBead));
        if(mag_ri_RA<10e-12)//possibly dividing by near zero, big numerical instabilities
            continue;
        // Compute functions
        double f= Function_f(ri, RA);
        vec<double> gradient_f=Function_gradient_f(ri, RA);
        double laplacian_f = Function_laplacian_f(ri, RA);
        //double f = 1; // TODO: Currently fixing f to 1
        //vec<double> gradient_f(zeros<vec<double>>(path.GetND()));
        //double laplacian_f = 0.;
        //double f = 1. + 2*z_a*(mag_ri_RA);
        //vec<double> gradient_f = 2*z_a*((ri_RA/mag_ri_RA));
        //double laplacian_f = 2*z_a*(path.GetND()-1)*((1./mag_ri_RA));

        // Sum over actions for ri
        std::vector<std::pair<uint32_t,uint32_t>> only_ri;
        only_ri.push_back(particle_pairs[pp_i][1]);
        vec<double> gradient_action(zeros<vec<double>>(path.GetND()));
        double laplacian_action = 0.;
        for (auto& action: action_list) {
          gradient_action += action->GetActionGradient(b_i,b_i+1,only_ri,0);
          laplacian_action += action->GetActionLaplacian(b_i,b_i+1,only_ri,0);
        }

        // Volume Term
        tot += (-1./mag_ri_RA*4.*M_PI)*(laplacian_f + f*(-laplacian_action + dot(gradient_action,gradient_action)) - 2.*dot(gradient_f,gradient_action));
        //Boundary Term
        int n_part_tot =path.species_list[species_a_i] + (species_a_i==species_b_i ? 0 : path.species_list[species_b_i]->n_part);
        if((!Checked[particle_pairs[pp_i][1].first])&&(mag_Delta_ri>0.8*path.GetL())) { //Boundary Event///TODO check because maybe I am summing here over all particles
            vec<double> NormalVector=getRelevantNormalVector(ri,ri_nextBead);
            for(int d=0;d<path.GetND();d++){//One has now to work with the picture of the particle in the other cell
                RA[d]+=NormalVector[d]*path.GetL();
            }
            ri_RA=path.Dr(ri,RA);
            mag_ri_RA=mag(ri_RA);
            if(mag_ri_RA<1e-5)//It acts in the 3 power in the following part, this can lead to numerical instabilities
                continue;
            vec<double> IntegrandVector=f*pow(mag_ri_RA,-3)*ri_RA+(f*gradient_action-gradient_f)/mag_ri_RA;//Compare calculation in "Calculation_Density_Estimator.pdf" Eq. (17)
            double VolumeFactor = path.GetVol()/path.GetSurface();//To correct the other measure
            tot+= VolumeFactor*dot(IntegrandVector,NormalVector);
        }
        Checked[particle_pairs[pp_i][1].first]=true;//Not to be checked anymore, because already calculated or not relevant
      }
    }
    total += tot;
    n_measure += 1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    total = 0;
    n_measure = 0;
  }

public:
  /// Constructor calls Init
  ContactDensity(Path &path, std::vector<std::shared_ptr<Action>>& t_action_list, Input &in, IO &out)
    : full_action_list(t_action_list), Observable(path, in, out, "scalar")
  {
    // Read in species info
    species_a = in.GetAttribute<std::string>("species_a");
    species_b = in.GetAttribute<std::string>("species_b");
    path.GetSpeciesInfo(species_a, species_a_i);
    path.GetSpeciesInfo(species_b, species_b_i);

    // Write things to file
    out.Write(prefix+"/species_a", species_a);
    out.Write(prefix+"/species_b", species_b);

    // Read in z_a
    z_a = in.GetAttribute<uint32_t>("z_a");

    // Generate action list
    std::vector<std::string> species_list;
    species_list.push_back(species_a);
    species_list.push_back(species_b);
    for (auto& action: full_action_list) {
      for (auto& sA: species_list) {
        if (std::find(action->species_list.begin(), action->species_list.end(), sA) != action->species_list.end()) {
          action_list.push_back(action);
          break;
        }
      }
    }
    //Set the improving function f
    //TODO find better f and also allow to choose differently
    Function_f=&Contact_Density_Optimization_Functions::f_simple; 
    Function_gradient_f=&Contact_Density_Optimization_Functions::gradient_f_simple;
    Function_laplace_f=&Contact_Density_Optimization_Functions::laplace_f_simple;

    Reset();
  }

  /// Write relevant information about an observable to the output
  virtual void Write()
  {
    if (n_measure > 0) {
      // Normalize
      uint32_t N_a = path.species_list[species_a_i]->n_part;
      uint32_t N_b = path.species_list[species_b_i]->n_part;
      double norm;
      if (species_a_i == species_b_i)
        norm = 0.5*n_measure*N_a*(N_a-1.)*path.GetNBead()/path.GetVol();
      else
        norm = n_measure*N_a*N_b*path.GetNBead()/path.GetVol();
      total /= norm;

      // Write to file
      if (first_time) {
        first_time = 0;
        out.CreateExtendableDataSet("/"+prefix, "x", total);
      } else {
        out.AppendDataSet("/"+prefix, "x", total);
      }

      Reset();
    }
  }
};

#endif // SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
