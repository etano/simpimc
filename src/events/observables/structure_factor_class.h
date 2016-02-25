#ifndef SIMPIMC_OBSERVABLES_STRUCTURE_FACTOR_CLASS_H_
#define SIMPIMC_OBSERVABLES_STRUCTURE_FACTOR_CLASS_H_

#include "observable_class.h"

/// Record the structure factor between two species of particles
class StructureFactor : public Observable {
   private:
    double k_cut;                        ///< Cutoff in k space
    std::shared_ptr<Species> species_a;  ///< First species
    std::shared_ptr<Species> species_b;  ///< Second species
    vec<double> sk;                      ///< Vector representing s(k)

    /// Accumulate the observable
    virtual void Accumulate() {
        path.SetMode(NEW_MODE);
        double cofactor = path.GetSign() * path.GetImportanceWeight();
        const auto &rho_k_a(species_a->GetRhoK());
        const auto &rho_k_b(species_b->GetRhoK());
        for (uint32_t k_i = 0; k_i < path.ks.mags.size(); k_i++) {
            if (path.ks.mags[k_i] < k_cut) {
                for (uint32_t b_i = 0; b_i < path.GetNBead(); ++b_i) {
                    sk(k_i) += cofactor * CMag2(rho_k_a(b_i)(k_i), rho_k_b(b_i)(k_i));
                }
            }
        }

        //if (species_b != species_a)
        //  sk = 2.*sk;

        n_measure += 1;
    }

    /// Reset the observable's counters
    virtual void Reset() {
        n_measure = 0;
        sk.zeros();
    }

   public:
    /// Constructor calls Init and sets output data_type
    StructureFactor(Path &path, Input &in, IO &out)
        : Observable(path, in, out, "histogram") {
        // Read in species info
        std::string species_a_name = in.GetAttribute<std::string>("species_a");
        std::string species_b_name = in.GetAttribute<std::string>("species_b");
        species_a = path.GetSpecies(species_a_name);
        species_b = path.GetSpecies(species_b_name);
        k_cut = in.GetAttribute<double>("k_cut", path.ks.cutoff);

        // Resize
        path.ks.Setup(k_cut);
        species_a->InitRhoK();
        species_b->InitRhoK();
        sk.zeros(path.ks.vecs.size());

        // Write things to file
        out.Write(prefix + "/species_a", species_a_name);
        out.Write(prefix + "/species_b", species_b_name);
        out.Write(prefix + "/k_cut", k_cut);
        out.CreateExtendableDataSet("/" + prefix + "/", "x", path.ks.mags[0]);
        for (uint32_t k_i = 1; k_i < path.ks.mags.size(); ++k_i)
            out.AppendDataSet("/" + prefix + "/", "x", path.ks.mags[k_i]);

        Reset();
    }

    /// Write relevant information about an observable to the output
    virtual void Write() {
        if (n_measure > 0) {
            // Normalize histograms
            double norm = n_measure * path.GetNBead() * species_a->GetNPart() * species_b->GetNPart();
            sk = sk / norm;

            // Write to file
            if (first_time) {
                first_time = 0;
                out.CreateExtendableDataSet("/" + prefix + "/", "y", sk);
            } else {
                out.AppendDataSet("/" + prefix + "/", "y", sk);
            }

            Reset();
        }
    }
};

#endif  // SIMPIMC_OBSERVABLES_STRUCTURE_FACTOR_CLASS_H_
