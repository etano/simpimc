#ifndef SIMPIMC_K_SPACE_CLASS_H_
#define SIMPIMC_K_SPACE_CLASS_H_

/// Container for all k space data
struct KSpace {
    double cutoff;                         ///< Largest cutoff used for k vectors
    double L;                              ///< Box side length
    uint32_t n_d;                          ///< Number of physical dimensions
    std::vector<double> mags;              ///< Vector holding magnitudes of k vectors
    std::vector<vec<double>> vecs;         ///< Vector holding k vectors
    std::vector<vec<int>> indices;         ///< Vector holding k index vectors
    vec<double> box;                       ///< Box defined by Brillouin zone
    vec<int> max_index;                    ///< Maximum k index used for each physical dimension
    field<vec<std::complex<double>>> c_k;  ///< Constant defined for each charge density

    /// Whether or not to include a k vector
    bool Include(const vec<double>& k, const double k_cut) {
        double k2 = dot(k, k);
        if (k2 < k_cut * k_cut && k2 != 0.) {
            if (k(0) > 0.)
                return true;
            else if (n_d >= 2 && (k(0) == 0. && k(1) > 0.))
                return true;
            else if (n_d >= 3 && (k(0) == 0. && k(1) == 0. && k(2) > 0.))
                return true;
            else
                return false;
        } else
            return false;
    }

    /// Setup k vectors
    void Setup(const double k_cut) {
        if (k_cut <= cutoff)
            return;
        else {
            vecs.clear();
            indices.clear();
            mags.clear();
            cutoff = k_cut;
        }

        // Calculate k box
        box.set_size(n_d);
        for (uint32_t d_i = 0; d_i < n_d; d_i++)
            box(d_i) = 2. * M_PI / L;

        // Calculate max k index based on k cutoff
        max_index.set_size(n_d);
        for (uint32_t d_i = 0; d_i < n_d; d_i++)
            max_index(d_i) = (uint32_t)ceil(1.1 * k_cut / box(d_i));

        // Set up C vector
        c_k.set_size(n_d);
        for (uint32_t d_i = 0; d_i < n_d; d_i++)
            c_k(d_i).set_size(2 * max_index(d_i) + 1);

        // Generate all possible combinations and permutations of k indices
        std::vector<int> is;
        for (uint32_t d_i = 0; d_i < n_d; d_i++)
            for (int i = -max_index(d_i); i <= max_index(d_i); i++)
                is.push_back(i);
        std::vector<std::vector<int>> t_indices;
        GenCombPermK(t_indices, is, n_d, false, true);

        // Iterate through indices, form k vectors, and include only those that should be included
        for (auto& ik : t_indices) {
            vec<double> k(n_d);
            vec<int> ki(n_d);
            for (uint32_t d_i = 0; d_i < n_d; d_i++) {
                k(d_i) = ik[d_i] * box(d_i);
                ki(d_i) = max_index(d_i) + ik[d_i];
            }
            if (Include(k, k_cut)) {
                vecs.push_back(k);
                indices.push_back(ki);
                mags.push_back(sqrt(dot(k, k)));
            }
        }
    }

    /// Calculate constant for each charge density
    void CalcC(const vec<double>& r) {
        for (uint32_t d_i = 0; d_i < box.size(); d_i++) {
            std::complex<double> tmp_c_k;
            double phi = r(d_i) * box(d_i);
            tmp_c_k = std::complex<double>(cos(phi), sin(phi));
            c_k(d_i)(max_index(d_i)) = 1.;
            for (uint32_t k_i = 1; k_i <= max_index(d_i); k_i++) {
                c_k(d_i)(max_index(d_i) + k_i) = tmp_c_k * c_k(d_i)(max_index(d_i) + k_i - 1);
                c_k(d_i)(max_index(d_i) - k_i) = conj(c_k(d_i)(max_index(d_i) + k_i));
            }
        }
    }
};

#endif  // SIMPIMC_K_SPACE_CLASS_H_
