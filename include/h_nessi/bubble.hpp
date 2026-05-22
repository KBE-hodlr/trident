
/**
 * @file bubble.hpp
 * @brief Declares functions for taking products of Keldysh components
 *
 */
#ifndef BUBBLE_DECL
#define BUBBLE_DECL

#include "utils.hpp"
#include "herm_matrix_hodlr.hpp"

/**
 * @namespace h_nessi
 * @brief Namespace for hierarchical matrix algorithms and data structures.
 */
namespace h_nessi {

  /**
   * @brief C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b2,b1}(t2,t1)
   */
  void Bubble1(int tstp, herm_matrix_hodlr &C, int c1, int c2, herm_matrix_hodlr &A, int a1, int a2, herm_matrix_hodlr &B, int b1, int b2, dlr_info &dlr); 

  /**
   * @brief C_{c1,c2}(t1,t2) = ii * A_{a1,a2}(t1,t2) * B_{b1,b2}(t1,t2)
   */
  void Bubble2(int tstp, herm_matrix_hodlr &C, int c1, int c2, herm_matrix_hodlr &A, int a1, int a2, herm_matrix_hodlr &B, int b1, int b2); 

}
#endif
