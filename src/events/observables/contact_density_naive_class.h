#ifndef SIMPIMC_OBSERVABLES_CONTACT_DENSITY_NAIVE_CLASS_H_
#define SIMPIMC_OBSERVABLES_CONTACT_DENSITY_NAIVE_CLASS_H_
#include "observable_class.h"

/// Record the pair correlation function between two species of particles
class ContactDensityNaive : public Observable {
   private:
    using particle_list_type = std::vector<std::pair<std::shared_ptr<Species>, uint32_t>>;
    double contact_density;                            ///< Contact density estimator
    std::shared_ptr<Species> species_a;                ///< First species
    std::shared_ptr<Species> species_b;                ///< Second species
    std::vector<particle_list_type> particle_pairs;    ///< Vector of particle pairs
    std::vector<std::shared_ptr<Action>> action_list;  ///< Vector of pointers to actions that involves species_a or species_b

    /// Accumulate the observable
    virtual void Accumulate() {
        double cofactor = path.GetSign() * path.GetImportanceWeight();
        for (uint32_t b_i = 0; b_i < path.GetNBead(); ++b_i) {
            for (const auto& p : particle_pairs) {
                // Move particle
                path.SetMode(NEW_MODE);
                vec<double> r_old = species_a->GetBead(p[0].second, b_i)->GetR();
                species_a->GetBead(p[0].second, b_i)->SetR(species_b->GetBead(p[1].second, b_i)->GetR());

                // Calculate action change
                double old_action = 0.;
                double new_action = 0.;
                for (auto& action : action_list) {
                    path.SetMode(OLD_MODE);
                    old_action += action->GetAction(b_i + path.GetNBead() - 1, b_i + path.GetNBead() + 1, p, 0);
                    path.SetMode(NEW_MODE);
                    new_action += action->GetAction(b_i + path.GetNBead() - 1, b_i + path.GetNBead() + 1, p, 0);
                }

                // Record exponential
                contact_density += cofactor * exp(-(new_action - old_action));

                // Move particle back
                path.SetMode(NEW_MODE);
                species_a->GetBead(p[0].second, b_i)->SetR(r_old);
            }
        }

        n_measure += 1;
    }

    /// Reset the observable's counters
    virtual void Reset() {
        n_measure = 0;
        contact_density = 0;
    }

   public:
    /// Constructor calls Init and sets output data_type
    ContactDensityNaive(Path& path, std::vector<std::shared_ptr<Action>>& full_action_list, Input& in, IO& out)
        : Observable(path, in, out, "scalar") {
        // Read in species info
        std::string species_a_name = in.GetAttribute<std::string>("species_a");
        std::string species_b_name = in.GetAttribute<std::string>("species_b");
        species_a = path.GetSpecies(species_a_name);
        species_b = path.GetSpecies(species_b_name);

        // Generate action list
        for (auto& action : full_action_list)
            if ((std::find(action->species_list.begin(), action->species_list.end(), species_a) != action->species_list.end()) or (std::find(action->species_list.begin(), action->species_list.end(), species_b) != action->species_list.end()))
                action_list.push_back(action);

        // Generate particle pairs
        if (species_a == species_b) {
            for (uint32_t p_i = 0; p_i < species_a->GetNPart() - 1; ++p_i) {
                for (uint32_t p_j = p_i + 1; p_j < species_b->GetNPart(); ++p_j) {
                    particle_list_type p;
                    p.push_back(std::make_pair(species_a, p_i));
                    p.push_back(std::make_pair(species_b, p_j));
                    particle_pairs.push_back(p);
                }
            }
        } else {
            for (uint32_t p_i = 0; p_i < species_a->GetNPart(); ++p_i) {
                for (uint32_t p_j = 0; p_j < species_b->GetNPart(); ++p_j) {
                    particle_list_type p;
                    p.push_back(std::make_pair(species_a, p_i));
                    p.push_back(std::make_pair(species_b, p_j));
                    particle_pairs.push_back(p);
                }
            }
        }

        // Write things to file
        out.Write(prefix + "/species_a", species_a_name);
        out.Write(prefix + "/species_b", species_b_name);

        Reset();
    }

    /// Write relevant information about an observable to the output
    virtual void Write() {
        if (n_measure > 0) {
            double norm = path.GetNBead() * n_measure * particle_pairs.size();
            contact_density /= norm;
            if (first_time) {
                first_time = 0;
                out.CreateExtendableDataSet(prefix, "x", contact_density);
            } else {
                out.AppendDataSet(prefix, "x", contact_density);
            }

            Reset();
        }
    }
};

#endif  // SIMPIMC_OBSERVABLES_CONTACT_DENSITY_NAIVE_CLASS_H_
