#pragma once

#include <vector>
#include <functional>
#include <complex>
#include <fftw3.h>
#include <omp.h>

#include "h_nessi/function.hpp"
#include "h_nessi/herm_matrix_hodlr.hpp"

using namespace h_nessi;

/**
 * @class Hubb_diss
 * @brief Implements the Second Born (2B) self-energy for the Hubbard model
 * using optimized FFTW executions and thread-local buffering.
 */
class Hubb_diss_all {
public:
    /**
     * @brief Constructor
     * @param U Interaction strength
     * @param L Linear dimension of the square lattice
     * @param Nk Number of irreducible k-points
     * @param nthreads Number of OpenMP threads for local execution
     */
    Hubb_diss_all(function &U);

    /**
     * @brief High-level dispatcher for self-energy calculation.
     * Manages MPI communication and routes to specific time-domain branches.
     */
    void Sigma_2b(int tstp, 
                     herm_matrix_hodlr &G, 
                     herm_matrix_hodlr &Sigma);


    void Sigma_hf(int tstp, 
                  herm_matrix_hodlr &G,
                  function &ellG);

private:
    // Constants
    function &U_; 

};
