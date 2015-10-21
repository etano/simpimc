#ifndef SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
#define SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
#include "observable_class.h"
namespace Contact_Density_Optimization_Functions {
    extern int ND=3;
    extern int z_a=1;
    
    double f_simple(const double &mag_ri_RA){
        return 1;
    }
    
    vec<double> gradient_f_simple(const double &mag_ri_RA, const vec<double> &ri_RA){
        return zeros<vec<double>>(ND);
    }
    double laplace_f_simple(const double &mag_ri_RA){
        return 0;
    }

    double f_Assaraf(const double &mag_ri_RA){//TODO include the exponential factor for the long range part (compare to Assarafs paper eq 11)
        return 1. + 2*z_a*mag_ri_RA;
    }
    
    vec<double> gradient_f_Assaraf(const double &mag_ri_RA, const vec<double> &ri_RA){
        return 2*z_a*(ri_RA/mag_ri_RA);
    }
    double laplace_f_Assaraf(const double &mag_ri_RA){
        return 4*z_a/mag_ri_RA;
        return 0;
    }
    
    double f_Squared(const double &mag_ri_RA){
        return 1+2*z_a*mag_ri_RA*mag_ri_RA;
    }
    vec<double> gradient_f_Squared(const double &mag_ri_RA, const vec<double> &ri_RA){
        return 4*z_a*ri_RA;
    }
    double laplace_f_Squared(const double &mag_ri_RA){
        //return 2*z_a*(ND-1)/mag(ri-RA);
        return 3*4*z_a;
    }

}

/// Measures the contact density between two species of particles. Taken from Assaraf, Caffarel, and Scemma. Phys Rev E 75, 035701(R) (2007). http://journals.aps.org/pre/pdf/10.1103/PhysRevE.75.035701.
class ContactDensity : public Observable
{
private:
    Histogram gr_vol; ///< Histogram representing g(r) arising from volume terms
    Histogram gr_b; ///< Histogram representing g(r) arising from boundary terms
    vec<int> n_measure_vol; ///< How many times the volume term at the i'th position gets measured
    vec<int> n_measure_b; ///< How many times the boundary term at the i'th position gets measured
    uint32_t z_a; ///< Charge of ion-like particle
    double lambda_tau; ///< the typical length of a path between two beads
    std::shared_ptr<Species> species_a; ///< ion species
    std::shared_ptr<Species> species_b; ///< other species 
    std::vector<std::shared_ptr<Action>> action_list; ///< Vector of pointers to actions that involves species_a or species_b
    std::vector<std::shared_ptr<Action>> &full_action_list; ///< Vector of pointers to all actions
    std::vector<std::pair<uint32_t,uint32_t>> particle_pairs; ///< contains all the pairs of particles of species_a and species_b
    size_t n_particle_pairs;
    std::string Optimization_Strategy;
    double (*Function_f)(const double &mag_ri_RA); ///< Function pointer for the possible generalization
    vec<double> (*Function_gradient_f)(const double &mag_ri_RA, const vec<double> &ri_RA);
    double (*Function_laplace_f)(const double &mag_ri_RA);
  
    vec<double> getRelevantNormalVector(vec<double> r1,vec<double> r2){
        vec<double> n(zeros<vec<double>>(path.GetND())); 
        for(int d=0;d<path.GetND();++d){
            if(r1[d]>(path.GetL()-2*lambda_tau)&&r2[d]<2*lambda_tau) {
                  n[d]=1;
            }
        else if(r2[d]>(path.GetL()-2*lambda_tau)&&r1[d]<2*lambda_tau) {
                n[d]=-1;
            }
        }
        return n;
      }
  
    /// Accumulate the observable
    virtual void Accumulate()
    {
        path.SetMode(NEW_MODE);
        // Add up contact probability
        vec<double> tot_vol(zeros<vec<double>>(gr_vol.x.n_r));//Save the values for the volume terms in here 
        vec<double> tot_b(zeros<vec<double>>(gr_vol.x.n_r)); //Save the values for the boundary terms in here
        for (uint32_t pp_i=0; pp_i<n_particle_pairs; ++pp_i) {
            for (uint32_t b_i=0; b_i<path.GetNBead(); ++b_i) {
                // Set r's
                vec<double> RA = species_a->GetBead(particle_pairs[pp_i].first,b_i)->GetR();
                vec<double> ri = species_b->GetBead(particle_pairs[pp_i].second,b_i)->GetR();
                vec<double> ri_nextBead = species_b->GetBead(particle_pairs[pp_i].second,b_i)->GetNextBead(1)->GetR();
                // Sum over actions for ri
                std::vector<std::pair<std::shared_ptr<Species>,uint32_t>> only_ri{std::make_pair(species_b,particle_pairs[pp_i].second)};
                vec<double> gradient_action(zeros<vec<double>>(path.GetND()));
                double laplacian_action = 0.;
                for (auto& action: action_list) {
                    gradient_action += action->GetActionGradient(b_i,b_i+1,only_ri,0);
                    laplacian_action+= action->GetActionLaplacian(b_i,b_i+1,only_ri,0);
                }
                vec<double> Direction(path.GetND());
                Direction.randn();
                Direction=Direction/norm(Direction);
                //Histogram loop
#pragma omp parallel for
                for (uint32_t i=0;i<gr_vol.x.n_r;++i){
                    vec<double> Rhist=gr_vol.x.rs(i)*Direction;
                    vec<double> R=path.Dr(RA,Rhist);
                    // Get differences
                    vec<double> ri_R(ri-R);
                    double mag_ri_R = mag(ri_R);
                    double mag_Delta_ri = mag(ri - ri_nextBead);
                    if(mag_ri_R<1e-6){//possibly dividing by near zero, big numerical instabilities therefore skip 
                        continue;
                    }
                    // Compute functions
                    double f= Function_f(mag_ri_R);
                    vec<double> gradient_f=Function_gradient_f(mag_ri_R, ri_R);
                    double laplacian_f = Function_laplace_f(mag_ri_R);

                    // Volume Term
                    //double tmp1=(-1./(mag_ri_R*4.*M_PI)) * laplacian_f;
                    //double tmp2=(-1./(mag_ri_R*4.*M_PI)) * (-f*laplacian_action);
                    //double tmp3=(-1./(mag_ri_R*4.*M_PI)) * (f*dot(gradient_action,gradient_action));
                    //double tmp4=(-1./(mag_ri_R*4.*M_PI)) * (-2*dot(gradient_f,gradient_action));
                    #pragma omp atomic
                    tot_vol(i) +=(-1./(mag_ri_R*4.*M_PI))*(laplacian_f + f*(-laplacian_action + dot(gradient_action,gradient_action)) - 2.*dot(gradient_f,gradient_action));
                    n_measure_vol(i)++;
                    //if(i==0) std::cout << Optimization_Strategy<<"\t"<<tmp1<< "\t"<<tmp2<<"\t"<<tmp3<<"\t"<<tmp4<<"\t"<<tmp1+tmp2+tmp3+tmp4<<std::endl;
                    //Boundary Term
                    if((mag_Delta_ri>2*lambda_tau)&&path.GetPBC()) { //Boundary Event
                        vec<double> NormalVector=getRelevantNormalVector(ri,ri_nextBead);
                        R+=NormalVector*path.GetL();//One has now to work with the picture of the particle in the other cell
                        ri_R=ri-R;
                        mag_ri_R=mag(ri_R);
                        if(mag_ri_R<1e-5)//It acts in the 3 power in the following part, this can lead to numerical instabilities
                            continue;
                        vec<double> IntegrandVector=f*pow(mag_ri_R,-3)*ri_R+(f*gradient_action-gradient_f)/mag_ri_R;//Compare calculation in "Calculation_Density_Estimator.pdf" Eq. (17)
                        double VolumeFactor = path.GetVol()/path.GetSurface();//To correct the other measure
                        #pragma omp atomic
                        tot_b(i)+= VolumeFactor*dot(IntegrandVector,NormalVector)/species_b->GetNPart();//if more then one ion is present, make sure to divide to normalize it correctly
                        n_measure_b(i)++;
                    }
                }
            }
        }
        double cofactor = path.GetSign()*path.GetImportanceWeight();
        for (uint32_t i=0;i<gr_vol.x.n_r;++i){
            gr_vol.y(i) += cofactor*tot_vol[i];
            gr_b.y(i) += cofactor*tot_b[i];
        }
    }

    /// Reset the observable's counters
    virtual void Reset()
    {
        gr_vol.y.zeros();
        n_measure_vol = zeros<vec<int>>(gr_vol.x.n_r);
        gr_b.y.zeros();
        n_measure_b = zeros<vec<int>>(gr_vol.x.n_r);
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
        out.Write(prefix+"/species_a", species_a_name);
        out.Write(prefix+"/species_b", species_b_name);
        // Read in grid info
        double r_min = in.GetAttribute<double>("r_min",0.);
        double r_max = in.GetAttribute<double>("r_max",path.GetL()/2.);
        int n_r = in.GetAttribute<double>("n_r",1000);
        gr_vol.x.CreateGrid(r_min,r_max,n_r);
        gr_vol.y.zeros(n_r);
        gr_b.x.CreateGrid(r_min,r_max,n_r);
        gr_b.y.zeros(n_r);
        //Write things to file
        std::string data_type = "histogram";
        out.Write(prefix+"/r_min", gr_vol.x.r_min);
        out.Write(prefix+"/r_max", gr_vol.x.r_max);
        out.Write(prefix+"/n_r", gr_vol.x.n_r);
        out.Write(prefix+"/x", gr_vol.x.rs);
        out.CreateGroup(prefix+"volume");
        out.Write(prefix+"volume/r_min", gr_vol.x.r_min);
        out.Write(prefix+"volume/r_max", gr_vol.x.r_max);
        out.Write(prefix+"volume/n_r", gr_vol.x.n_r);
        out.Write(prefix+"volume/x", gr_vol.x.rs);
        out.Write(prefix+"volume/data_type",data_type);
        out.CreateGroup(prefix+"boundary");
        out.Write(prefix+"boundary/r_min", gr_vol.x.r_min);
        out.Write(prefix+"boundary/r_max", gr_vol.x.r_max);
        out.Write(prefix+"boundary/n_r", gr_vol.x.n_r);
        out.Write(prefix+"boundary/x", gr_vol.x.rs);
        out.Write(prefix+"boundary/data_type",data_type);
        // Read in z_a
        z_a = in.GetAttribute<uint32_t>("z_a");
        
        //Get lambda_tau
        lambda_tau=path.GetTau()*species_b->GetLambda(); 
        // Generate action list
        for (auto& action: full_action_list)
            if ((std::find(action->species_list.begin(), action->species_list.end(), species_a)!=action->species_list.end()) or (std::find(action->species_list.begin(), action->species_list.end(), species_b)!=action->species_list.end()))
              action_list.push_back(action);

        // Form particle pairs 
        if (species_a == species_b) { // Homogeneous
            for (uint32_t p_i=0; p_i<species_a->GetNPart()-1; ++p_i)
                for (uint32_t p_j=p_i+1; p_j<species_b->GetNPart(); ++p_j)
                    particle_pairs.push_back(std::make_pair(p_i,p_j));
        } 
        else { // Homologous
            for (uint32_t p_i=0; p_i<species_a->GetNPart(); ++p_i)
                for (uint32_t p_j=0; p_j<species_b->GetNPart(); ++p_j)
                    particle_pairs.push_back(std::make_pair(p_i,p_j));
        }
        n_particle_pairs = particle_pairs.size();
    
        //Choose the optimization stategy
        Optimization_Strategy = in.GetAttribute<std::string>("optimization_strategy");
        Contact_Density_Optimization_Functions::ND=path.GetND();
        Contact_Density_Optimization_Functions::z_a=z_a;
        if(Optimization_Strategy=="Simple"){
            Function_f=&Contact_Density_Optimization_Functions::f_simple; 
            Function_gradient_f=&Contact_Density_Optimization_Functions::gradient_f_simple;
            Function_laplace_f=&Contact_Density_Optimization_Functions::laplace_f_simple;
        }
        else if (Optimization_Strategy=="Assaraf"){
            Function_f=&Contact_Density_Optimization_Functions::f_Assaraf; 
            Function_gradient_f=&Contact_Density_Optimization_Functions::gradient_f_Assaraf;
            Function_laplace_f=&Contact_Density_Optimization_Functions::laplace_f_Assaraf;
        }
        else if (Optimization_Strategy=="Squared"){
            Function_f=&Contact_Density_Optimization_Functions::f_Squared; 
            Function_gradient_f=&Contact_Density_Optimization_Functions::gradient_f_Squared;
            Function_laplace_f=&Contact_Density_Optimization_Functions::laplace_f_Squared;
        }
        else {
            std::cout << "ERROR: unrecognized optimization strategy in contact density, please check again"<< std::endl;
            assert(false);
        }
        Reset();
    }

    /// Write relevant information about an observable to the output
    virtual void Write()
    {
        if (min(n_measure_vol+n_measure_b) > 0) {
            vec<double> tot(zeros<vec<double>>(gr_vol.x.n_r)); 
            // Normalize
            double norm_vol;
            double norm_b;
            //if (species_a == species_b)
            //  norm = 0.5*n_measure*species_a->GetNPart()*(species_a->GetNPart()-1.)*path.GetNBead()/path.GetVol();
            //else
            //  norm = n_measure*species_a->GetNPart()*species_b->GetNPart()*path.GetNBead()/path.GetVol();
            for (uint32_t i=0; i<gr_vol.x.n_r; i++) {
                //double r1 = gr_vol.x(i);
                //double r2 = (i<(gr_vol.x.n_r-1)) ? gr_vol.x(i+1):(2.*gr_vol.x(i)-gr_vol.x(i-1));
                //double r = 0.5*(r1+r2);
                //double bin_vol = r2-r1;
                //double bin_vol=-1;
                //if (path.GetND() == 3)
                //double bin_vol = 4.*M_PI/3. * (r2*r2*r2-r1*r1*r1);
                //else if (path.GetND() == 2)
                //    bin_vol = M_PI * (r2*r2-r1*r1);
                //else if (path.GetND() == 1)
                //    bin_vol = r2-r1;
                //assert(bin_vol!=-1);//Make sure one case is fulfilled
                norm_vol = n_measure_vol[i]/path.GetVol();
                //if(i==0) std::cout << Optimization_Strategy<<"\t"<<gr_vol.y(i)<< "\t" << bin_vol<< "\t" << norm_vol << "\t" << gr_vol.y(i)/(norm_vol) << std::endl;
                //gr_vol.y(i) = gr_vol.y(i)/(bin_vol*norm_vol);
                gr_vol.y(i) = gr_vol.y(i)/(norm_vol);
                if(path.GetPBC()){//If we do not have pbc then we just save 0, otherwise similar to volume term
                    norm_b = n_measure_b[i]/path.GetVol();
                    gr_b.y(i) = gr_b.y(i)/(norm_b);
                }
                tot(i)=gr_vol.y(i)+gr_b.y(i);
            }
            // Write to file
            if (first_time) {
                first_time = 0;
                out.CreateExtendableDataSet("/"+prefix,"y",tot);
                out.CreateExtendableDataSet("/"+prefix+"volume/", "y", gr_vol.y);
                out.CreateExtendableDataSet("/"+prefix+"boundary/", "y", gr_b.y);
            } else {
                out.AppendDataSet("/"+prefix,"y",tot);
                out.AppendDataSet("/"+prefix+"volume/", "y", gr_vol.y);
                out.AppendDataSet("/"+prefix+"boundary/", "y", gr_b.y);
              }

            Reset();
        }
    }
};

#endif // SIMPIMC_OBSERVABLES_CONTACT_DENSITY_CLASS_H_
