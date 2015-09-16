#ifndef SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
#define SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_

#include "observable_class.h"

namespace Contact_Density_Optimization_Functions {
    int ND=3;
    double f_simple(vec<double> ri, vec<double> RA){
        return 1;
    }
    
    vec<double> gradient_f_simple(vec<double> ri, vec<double> RA){
        return zeros<vec<double>>(ND);
    }
    double laplace_f_simple(vec<double> ri, vec<double> RA){
        return 0;
    }
}

/// Measures the contact density between two species of particles. Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
class ContactDensity : public Observable
{
private:
  Histogram gr; ///< Histogram representing g(r)
  uint32_t z_a; ///< Charge of ion-like particle
  double r_min; ///< minimal r in the histogram
  double r_max; ///< maximal r in the histogram
  int n_r; ///< number of bins in the histogram
  std::shared_ptr<Species> species_a; ///< ion species
  std::shared_ptr<Species> species_b; ///< other species 
  std::vector<std::shared_ptr<Action>> action_list; ///< Vector of pointers to actions that involves species_a or species_b
  std::vector<std::shared_ptr<Action>> &full_action_list; ///< Vector of pointers to all actions
  double (*Function_f)(vec<double> ri, vec<double> RA); ///< Function pointer for the possible generalization
  vec<double> (*Function_gradient_f)(vec<double> ri, vec<double> RA);
  double (*Function_laplace_f)(vec<double> ri, vec<double> RA);
  
  vec<double> getRelevantNormalVector(vec<double> r1,vec<double> r2){
    vec<double> n(zeros<vec<double>>(path.GetND())); 
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

  ///Get the relevant histogram running variable R
  inline double histR(double minimal_R, double maximal_R, int Number_R, int i){
    return minimal_R+(maximal_R-minimal_R)*(i*1.0/Number_R);
  }
  
  /// Accumulate the observable
  virtual void Accumulate()
  {
    path.SetMode(NEW_MODE);
    //TODO isn't the next part (forming the pairs) actually doable in the constructor, to speed up the calculation?
    // Form particle pairs
    std::vector<std::pair<uint32_t,uint32_t>> particle_pairs;
    if (species_a == species_b) { // Homogeneous
      for (uint32_t p_i=0; p_i<species_a->GetNPart()-1; ++p_i)
        for (uint32_t p_j=p_i+1; p_j<species_b->GetNPart(); ++p_j)
          particle_pairs.push_back(std::make_pair(p_i,p_j));
    } else { // Homologous
      for (uint32_t p_i=0; p_i<species_a->GetNPart(); ++p_i)
        for (uint32_t p_j=0; p_j<species_b->GetNPart(); ++p_j)
          particle_pairs.push_back(std::make_pair(p_i,p_j));
    }
    // Add up contact probability
    double * tot = new double[n_r]();//save here the values for the histogram
    assert(tot[0]==0); //Check if initialization works as assumed
    size_t n_particle_pairs(particle_pairs.size());
    #pragma omp parallel for 
    for (uint32_t pp_i=0; pp_i<n_particle_pairs; ++pp_i) {
      for (uint32_t b_i=0; b_i<path.GetNBead(); ++b_i) {
	    // Set r's
	    vec<double> RA = species_a->GetBead(particle_pairs[pp_i].first,b_i)->GetR();
	    vec<double> ri = species_b->GetBead(particle_pairs[pp_i].second,b_i)->GetR();
	    vec<double> ri_nextBead = species_b->GetBead(particle_pairs[pp_i].second,b_i+1)->GetR();
		//Histogram loop
        for (uint32_t i=0;i<n_r;++i){
            vec<double> Rhist(zeros<vec<double>>(path.GetND()));
            if(RA[0]<path.GetL()/2.0)//A bit of hacking, however this should be a small bit faster
                Rhist[0]=histR(r_min,r_max, n_r, i); //TODO check with etano if just measuring along x direction is fine or not
            else
                Rhist[0]=-histR(r_min,r_max, n_r, i);
            vec<double> R=Rhist+RA;
		    // Get differences
		    vec<double> ri_R(path.Dr(ri, R));
		    double mag_ri_R = mag(ri_R);
		    double mag_Delta_ri = mag(path.Dr(ri,ri_nextBead));
		    if(mag_ri_R<1e-12)//possibly dividing by near zero, big numerical instabilities //TODO check if this could lead to normalization problems, because for some R they are counted and for others not
		        continue;
		    // Compute functions
		    double f= Function_f(ri, R);
		    vec<double> gradient_f=Function_gradient_f(ri, R);
		    double laplacian_f = Function_laplace_f(ri, R);
		    //double f = 1; 
		    //vec<double> gradient_f(zeros<vec<double>>(path.GetND()));
		    //double laplacian_f = 0.;
		    //double f = 1. + 2*z_a*(mag_ri_RA);
		    //vec<double> gradient_f = 2*z_a*((ri_RA/mag_ri_RA));
		    //double laplacian_f = 2*z_a*(path.GetND()-1)*((1./mag_ri_RA));

		    // Sum over actions for ri
		    std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> only_ri{std::make_pair(species_a,particle_pairs[pp_i].second)};
		    vec<double> gradient_action(zeros<vec<double>>(path.GetND()));
		    double laplacian_action = 0.;
		    for (auto& action: action_list) {
		      gradient_action += action->GetActionGradient(b_i,b_i+1,only_ri,0);
		      laplacian_action += action->GetActionLaplacian(b_i,b_i+1,only_ri,0);
		    }

		    // Volume Term
            #pragma omp atomic
		    tot[i] += (-1./mag_ri_R*4.*M_PI)*(laplacian_f + f*(-laplacian_action + dot(gradient_action,gradient_action)) - 2.*dot(gradient_f,gradient_action));
		    //Boundary Term
		    if(mag_Delta_ri>0.8*path.GetL()) { //Boundary Event
		        vec<double> NormalVector=getRelevantNormalVector(ri,ri_nextBead);
		        for(int d=0;d<path.GetND();d++){//One has now to work with the picture of the particle in the other cell
		            R[d]+=NormalVector[d]*path.GetL();
		        }
		        ri_R=path.Dr(ri,R);
		        mag_ri_R=mag(ri_R);
		        if(mag_ri_R<1e-5)//It acts in the 3 power in the following part, this can lead to numerical instabilities
		            continue;
		        vec<double> IntegrandVector=f*pow(mag_ri_R,-3)*ri_R+(f*gradient_action-gradient_f)/mag_ri_R;//Compare calculation in "Calculation_Density_Estimator.pdf" Eq. (17)
		        double VolumeFactor = path.GetVol()/path.GetSurface();//To correct the other measurei
                #pragma omp atomic
		        tot[i]+= VolumeFactor*dot(IntegrandVector,NormalVector)/species_b->GetNPart();//if more then one ion is present, make sure to divide to normalize it correctly
		    }
		 }
      }
    }
    double cofactor = path.GetSign()*path.GetImportanceWeight();
    for (uint32_t i=0;i<n_r;++i){
        gr.y(i) += cofactor*tot[i];
    }
    delete[] tot;
    tot=0;
    n_measure+=1;
  }

  /// Reset the observable's counters
  virtual void Reset()
  {
    gr.y.zeros();
    n_measure = 0;
  }

public:
  /// Constructor calls Init
  ContactDensity(Path &path, std::vector<std::shared_ptr<Action>>& t_action_list, Input &in, IO &out)
    : full_action_list(t_action_list), Observable(path, in, out, "histogram")
  {
    // Read in species info
    std::string species_a_name = in.GetAttribute<std::string>("species_a");
    std::string species_b_name = in.GetAttribute<std::string>("species_b");
    species_a = path.GetSpecies(species_a_name);
    species_b = path.GetSpecies(species_b_name);

    // Write things to file
    out.Write(prefix+"/species_a", species_a);
    out.Write(prefix+"/species_b", species_b);

    // Read in grid info
    r_min = in.GetAttribute<double>("r_min",0.);
    r_max = in.GetAttribute<double>("r_max",path.GetL()/2.);
    n_r = in.GetAttribute<double>("n_r",1000);
    gr.x.CreateGrid(r_min,r_max,n_r);
    gr.y.zeros(n_r);

    //Write things to file
    out.Write(prefix+"/r_min", r_min);
    out.Write(prefix+"/r_max", r_max);
    out.Write(prefix+"/n_r", n_r);
    // Read in z_a
    z_a = in.GetAttribute<uint32_t>("z_a");

    // Generate action list
    for (auto& action: full_action_list)
        if ((std::find(action->species_list.begin(), action->species_list.end(), species_a)!=action->species_list.end()) or (std::find(action->species_list.begin(), action->species_list.end(), species_b)!=action->species_list.end()))
          action_list.push_back(action);

    //Set the improving function f
    //TODO find better f and also allow to choose differently
    Contact_Density_Optimization_Functions::ND=path.GetND();
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
      double norm;
      if (species_a == species_b)
        norm = 0.5*n_measure*species_a->GetNPart()*(species_a->GetNPart()-1.)*path.GetNBead()/path.GetVol();//TODO normalization issue mentioned above in Accumulate
      else
        norm = n_measure*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead()/path.GetVol();
      for (uint32_t i=0; i<gr.x.n_r; i++) {
        double r1 = gr.x(i);
        double r2 = (i<(gr.x.n_r-1)) ? gr.x(i+1):(2.*gr.x(i)-gr.x(i-1));
        double r = 0.5*(r1+r2);
        double bin_vol=-1;
        if (path.GetND() == 3)
          bin_vol = 4.*M_PI/3. * (r2*r2*r2-r1*r1*r1);
        else if (path.GetND() == 2)
          bin_vol = M_PI * (r2*r2-r1*r1);
        else if (path.GetND() == 1)
          bin_vol = r2-r1;
        assert(bin_vol!=-1);//Make sure one case is fulfilled
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

#endif // SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
