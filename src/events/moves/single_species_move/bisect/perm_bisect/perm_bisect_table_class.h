#ifndef SIMPIMC_MOVES_PERM_TABLE_BISECT_CLASS_H_
#define SIMPIMC_MOVES_PERM_TABLE_BISECT_CLASS_H_

#include "perm_bisect_class.h"

/// Permuting bisection move which samples a permutation table
class PermBisectTable : public PermBisect {
   private:
    /// Construct the permutation table
    double ConstructPermTable() {
        // Update t
        UpdatePermTable();

        // Run through permatation types
        double total_weight = 0.;
        for (uint32_t perm_index = 0; perm_index < all_cycles.size(); perm_index++) {
            Cycle &c = all_cycles(perm_index);
            c.weight = 1.;
            for (uint32_t p_i = 0; p_i < c.part.size(); p_i++)
                c.weight *= t(c.part(p_i), c.part(c.perm(p_i)));
            if (c.weight > epsilon) {
                total_weight += c.weight;
                c.contribution = total_weight;
                cycles.push_back(&c);
            }
        }
        return total_weight;
    }

    /// Update the permutation table
    void UpdatePermTable() {
        // Set initial and final beads
        field<std::shared_ptr<Bead>> b0(species->GetNPart()), b1(species->GetNPart());
        for (uint32_t p_i = 0; p_i < species->GetNPart(); p_i++) {
            b0(p_i) = species->GetBead(p_i, bead0);
            b1(p_i) = b0(p_i)->GetNextBead(n_bisect_beads);
        }

        // Construct t table
        for (uint32_t i = 0; i < species->GetNPart(); i++) {
            vec<double> dr_ii(path.Dr(b0(i), b1(i)));
            for (uint32_t j = 0; j < species->GetNPart(); j++) {
                vec<double> dr_ij(path.Dr(b0(i), b1(j)));
                double exponent = (-dot(dr_ij, dr_ij) + dot(dr_ii, dr_ii)) * i_4_lambda_tau_n_bisect_beads;
                if (exponent > log_epsilon)
                    t(i, j) = exp(exponent);
                else
                    t(i, j) = 0.;
            }
        }
    }

    /// Build all possible cycles
    void BuildCycles() {
        // Generate possible cycles of n_perm_part particles
        std::vector<uint32_t> tmp_poss_cycle;
        for (uint32_t i = 0; i < n_perm_part; ++i)
            tmp_poss_cycle.push_back(i);
        std::vector<std::vector<uint32_t>> tmp_poss_cycles;
        GenPerm(tmp_poss_cycles, tmp_poss_cycle);
        n_perm_type = tmp_poss_cycles.size();
        field<Cycle> possible_cycles;
        possible_cycles.set_size(n_perm_type);
        for (uint32_t i = 0; i < n_perm_type; ++i) {
            mat<uint32_t> p_i(zeros<mat<uint32_t>>(n_perm_part, n_perm_part));
            possible_cycles(i).perm.set_size(n_perm_part);
            for (uint32_t j = 0; j < n_perm_part; ++j) {
                possible_cycles(i).perm(j) = tmp_poss_cycles[i][j];
                for (uint32_t k = 0; k < n_perm_part; ++k)
                    if (tmp_poss_cycles[i][k] == j)
                        p_i(j, k) = 1;
            }
            vec<uint32_t> ic(n_perm_part);
            for (uint32_t j = 0; j < n_perm_part; ++j)
                ic(j) = tmp_poss_cycle[j];
            ic = p_i * ic;
            possible_cycles(i).i_perm = ic;
        }

        // Run through permatation types
        std::vector<uint32_t> tmp_cycle;
        for (uint32_t i = 0; i < species->GetNPart(); ++i)
            tmp_cycle.push_back(i);
        std::vector<std::vector<uint32_t>> tmp_cycles;
        GenCombPermK(tmp_cycles, tmp_cycle, n_perm_part, false, false);
        uint32_t n_cycle = tmp_cycles.size();
        all_cycles.set_size(n_perm_type * n_cycle);
        uint32_t perm_index = 0;
        for (uint32_t perm_type_i = 0; perm_type_i < n_perm_type; perm_type_i++) {
            for (uint32_t cycle_i = 0; cycle_i < n_cycle; ++cycle_i) {
                Cycle &c = all_cycles(perm_index);
                c.type = perm_type_i;
                c.index = perm_index;
                c.part.set_size(tmp_cycles[cycle_i].size());
                for (uint32_t i = 0; i < tmp_cycles[cycle_i].size(); ++i)
                    c.part(i) = tmp_cycles[cycle_i][i];
                c.perm = possible_cycles(perm_type_i).perm;
                c.i_perm = possible_cycles(perm_type_i).i_perm;
                perm_index += 1;
            }
        }
    }

    /// Select which cycle to perform a permutation with
    uint32_t SelectCycle(const double permTot) {
        double x = rng.UnifRand(0., permTot);
        uint32_t hi = cycles.size();
        uint32_t lo = 0;
        if (x < cycles[0]->contribution)
            return 0;
        while (hi - lo > 1) {
            uint32_t mid = (hi + lo) >> 1;
            if (x < cycles[mid]->contribution)
                hi = mid;
            else
                lo = mid;
        }
        return hi;
    }

    /// Attempt the move
    virtual bool Attempt() {
        bead0 = rng.UnifRand(species->GetNBead()) - 1;       // Pick first bead at random
        bead1 = bead0 + n_bisect_beads;                      // Set last bead in bisection
        bool roll_over = bead1 > (species->GetNBead() - 1);  // See if bisection overflows to next particle
        // Set up permutation
        cycles.clear();
        double perm_tot_0 = ConstructPermTable();  // Permutation weight table
        uint32_t cycle_index = SelectCycle(perm_tot_0);
        Cycle *c = cycles[cycle_index];  // TODO: Make a smart pointer

        // Set up pointers
        std::vector<std::pair<std::shared_ptr<Species>, uint32_t>> particles;
        n_perm_part = c->part.size();
        field<std::shared_ptr<Bead>> bead_i(n_perm_part), bead_fm1(n_perm_part), bead_f(n_perm_part);
        for (uint32_t i = 0; i < n_perm_part; i++) {
            bead_i(i) = species->GetBead(c->part(i), bead0);
            bead_fm1(i) = bead_i(i)->GetNextBead(n_bisect_beads - 1);
            bead_f(i) = bead_fm1(i)->GetNextBead(1);
            particles.push_back(std::make_pair(species, c->part(i)));
        }

        // Permute particles
        PermuteBeads(bead_fm1, bead_f, *c);

        // Note affected beads
        field<std::shared_ptr<Bead>> bead_a(n_perm_part);
        affected_beads.clear();
        for (uint32_t i = 0; i < n_perm_part; i++) {
            for (bead_a(i) = bead_i(i); bead_a(i) != bead_f(i); bead_a(i) = bead_a(i)->GetNextBead(1))
                affected_beads.push_back(bead_a(i));
        }

        // Perform the bisection (move exactly through kinetic action)
        field<std::shared_ptr<Bead>> bead_b(n_perm_part), bead_c(n_perm_part);
        double prev_action_change = -log(c->weight);
        for (int level_i = n_level - 1; level_i >= 0; level_i -= 1) {
            // Level specific quantities
            uint32_t skip = 1 << level_i;
            double levelTau = path.GetTau() * skip;
            double sigma2 = species->GetLambda() * levelTau;
            double sigma = sqrt(sigma2);

            // Calculate sampling probability
            double old_log_sample_prob = 0.;
            double new_log_sample_prob = 0.;
            for (uint32_t i = 0; i < n_perm_part; i++) {
                bead_a(i) = bead_i(i);
                while (bead_a(i) != bead_f(i)) {
                    // Old sampling
                    path.SetMode(OLD_MODE);
                    bead_b(i) = bead_a(i)->GetNextBead(skip);
                    bead_c(i) = bead_b(i)->GetNextBead(skip);
                    vec<double> delta_old(path.Dr(bead_b(i), path.RBar(bead_c(i), bead_a(i))));
                    old_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_old);

                    // New sampling
                    path.SetMode(NEW_MODE);
                    bead_b(i) = bead_a(i)->GetNextBead(skip);
                    bead_c(i) = bead_b(i)->GetNextBead(skip);
                    vec<double> delta_new(path.GetND());
                    rng.NormRand(delta_new, 0., sigma);
                    path.PutInBox(delta_new);
                    bead_b(i)->SetR(path.RBar(bead_c(i), bead_a(i)) + delta_new);
                    new_log_sample_prob += rho_free_splines[level_i].GetLogRhoFree(delta_new);

                    bead_a(i) = bead_c(i);
                }
            }

            // Calculate action change
            double old_action = 0.;
            double new_action = 0.;
            for (auto &action : action_list) {
                // Old action
                path.SetMode(OLD_MODE);
                old_action += action->GetAction(bead0, bead1, particles, level_i);

                // New action
                path.SetMode(NEW_MODE);
                new_action += action->GetAction(bead0, bead1, particles, level_i);
            }

            // Calculate acceptance ratio
            double log_sample_ratio = -new_log_sample_prob + old_log_sample_prob;
            double current_action_change = new_action - old_action;
            double log_accept_probablity = log_sample_ratio - current_action_change + prev_action_change;

            // Metropolis step
            if (log_accept_probablity < log(rng.UnifRand()))
                return 0;

            prev_action_change = current_action_change;
        }

        // Have to check this if neighborhood has changed
        if (n_perm_part < species->GetNPart()) {
            // Construct Permutation Table
            path.SetMode(NEW_MODE);
            double perm_tot_1 = ConstructPermTable();

            // Decide whether or not to accept whole bisection
            if ((perm_tot_0 / perm_tot_1) < rng.UnifRand())
                return 0;
        }

        return 1;
    }

   public:
    /// Constructor instantiates parent class, sets the number of particles involved in a permutation, and builds all possible cycles
    PermBisectTable(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
        : PermBisect(path, rng, action_list, in, out) {
        // Read in things
        n_perm_part = in.GetAttribute<uint32_t>("n_perm_part");

        // Generate all cycles
        BuildCycles();
    }
};

#endif  // SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
