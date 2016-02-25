#ifndef SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_
#define SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_

#include "nodal_class.h"
#include "../../free_spline_class.h"

/// Bare free particle nodal action class
class FreeNodal : public Nodal {
   private:
    std::vector<FreeSpline> rho_free_splines;  ///< Holds the splined action for every time slice

    /// Creates splined action for all time slices
    virtual void SetupSpline() {
        // Create splines
        uint32_t n_spline = species->GetNBead() / 2 + (species->GetNBead() % 2) + 1;
        rho_free_splines.resize(n_spline);
#pragma omp parallel for
        for (uint32_t spline_i = 0; spline_i < n_spline; ++spline_i)
            rho_free_splines[spline_i] = FreeSpline(path.GetL(), n_images, species->GetLambda(), path.GetTau() * (spline_i + 1), false);
    }

    /// Returns the value of g_ij
    virtual double GetGij(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff) {
        return rho_free_splines[slice_diff - 1].GetRhoFree(path.Dr(b_i, b_j));
    }

    /// Returns the spatial derivative of g_ij
    virtual double GetGijDGijDr(const std::shared_ptr<Bead> &b_i, const std::shared_ptr<Bead> &b_j, const uint32_t slice_diff, vec<double> &dgij_dr) {
        return rho_free_splines[slice_diff - 1].GetGradRhoFree(path.Dr(b_i, b_j), dgij_dr);
    }

   public:
    // Constructor calls Init
    FreeNodal(Path &path, Input &in, IO &out)
        : Nodal(path, in, out) {
        // Setup splines
        SetupSpline();

        // Test
        bool init_good = TestNodes();
        out.Write("Actions/" + name + "/init_good", init_good);
    }
};

#endif  // SIMPIMC_ACTIONS_FREE_NODAL_CLASS_H_
