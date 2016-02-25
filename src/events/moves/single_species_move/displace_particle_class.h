#ifndef SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
#define SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_

#include "single_species_move_class.h"

/// This moves displaces whole particles by a defined step size
class DisplaceParticle : public SingleSpeciesMove {
   private:
    double step_size;                                   ///< Distanced displaced during move
    std::vector<std::shared_ptr<Bead>> affected_beads;  ///< Beads affected by the move

    /// Accept the move
    virtual void Accept() {
        // Move Accepted, so copy new coordinates
        for (auto& b : affected_beads) {
            b->StoreR();
            b->StoreRhoK();
        }
        for (uint32_t b_i = 0; b_i < species->GetNBead(); ++b_i)
            species->StoreRhoK(b_i);

        // Call accept for each action
        for (auto& action : action_list)
            action->Accept();
    }

    /// Attempt the move
    virtual bool Attempt() {
        // Set which particles are affected by the move
        uint32_t p_i = rng.UnifRand(species->GetNPart()) - 1;  // Pick particle at random
        std::vector<std::pair<std::shared_ptr<Species>, uint32_t>> particles;
        particles.push_back(std::make_pair(species, p_i));

        // New sampling
        path.SetMode(NEW_MODE);
        vec<double> dr(path.GetND());
        rng.UnifRand(dr, step_size);

        // Set which beads are affected by the move
        // and move them
        affected_beads.clear();
        std::shared_ptr<Bead> beadA(species->GetBead(p_i, 0));
        std::shared_ptr<Bead> beadF(beadA->GetNextBead(species->GetNBead() - 1));
        while (beadA != beadF) {
            affected_beads.push_back(beadA);
            beadA->SetR(beadA->GetR() + dr);
            beadA = beadA->GetNextBead(1);
        }
        affected_beads.push_back(beadF);
        beadF->SetR(beadF->GetR() + dr);

        // Calculate action change
        double old_action = 0.;
        double new_action = 0.;
        for (auto& action : action_list) {
            // Old action
            path.SetMode(OLD_MODE);
            old_action += action->GetAction(0, species->GetNBead(), particles, 0);

            // New action
            path.SetMode(NEW_MODE);
            new_action += action->GetAction(0, species->GetNBead(), particles, 0);
        }

        double log_accept_probablity = old_action - new_action;

        // Metropolis reject step
        if (log_accept_probablity < log(rng.UnifRand()))
            return 0;
        else
            return 1;
    }

    /// Rejects the move
    virtual void Reject() {
        // Move rejected, so return old coordinates
        for (auto& b : affected_beads) {
            b->RestoreR();
            b->RestoreRhoK();
        }
        for (uint32_t b_i = 0; b_i < species->GetNBead(); ++b_i)
            species->RestoreRhoK(b_i);

        // Call reject for each action
        for (auto& action : action_list)
            action->Reject();
    }

   public:
    /// Constructor instantiates parent and sets the step size
    DisplaceParticle(Path& path, RNG& rng, std::vector<std::shared_ptr<Action>>& action_list, Input& in, IO& out)
        : SingleSpeciesMove(path, rng, action_list, in, out) {
        step_size = in.GetAttribute<double>("step_size", path.GetL() / 10.);
    }
};

#endif  // SIMPIMC_MOVES_DISPLACE_PARTICLE_CLASS_H_
