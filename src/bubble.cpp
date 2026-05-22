#ifndef BUBBLE_IMPL
#define BUBBLE_IMPL

#include "h_nessi/bubble.hpp"

namespace h_nessi {

  void Bubble1(int tstp, herm_matrix_hodlr &C, int c1, int c2, herm_matrix_hodlr &A, int a1, int a2, herm_matrix_hodlr &B, int b1, int b2, dlr_info &dlr) {
    if(tstp == -1) {
      // C_{c1,c2}(tau) = -    A_{a1,a2}(tau) * B_{b2,b1}(-tau)
      //                = -\xi A_{a1,a2}(tau) * B_{b2,b1}(beta-tau)

      int ntau = B.ntau();
      int sig = -1 * B.sig();
      DMatrix BM_reversed(ntau, 1);
      B.get_mat_reversed(dlr, BM_reversed, b2, b1);

      for(int tau = 0; tau < ntau; tau++) {
        C.map_mat(tau)(c1,c2) = sig * A.map_mat(tau)(a1,a2) * BM_reversed(tau,0);
      }
    }
    else {
      cplx cplxi(0.,1.);

      // C^\rceil_{c1,c2}(t,tau) = iA^\rceil_{a1,a2}(t,tau) * B^\lceil_{b2,b1}(tau,t)
      int ntau = C.ntau();
      ZMatrix Bvt(ntau, 1);
      B.get_vt(tstp, dlr, Bvt, b2, b1);
    
      for(int tau = 0; tau < ntau; tau++) {
        C.map_tv(tstp,tau)(c1,c2) = cplxi * A.map_tv(tstp,tau)(a1,a2) * Bvt(tau,0);
      }

      for(int t = 0; t <= tstp; t++) {
        cplx ALtT12 = A.map_les_curr(t,tstp)(a1,a2);
        cplx BGTt21 = B.map_ret_curr(tstp,t)(b2,b1) - std::conj(B.map_les_curr(t,tstp)(b1,b2));
        
        cplx ALtT21 = A.map_les_curr(t,tstp)(a2,a1);
        cplx BLtT21 = B.map_les_curr(t,tstp)(b2,b1);
        cplx ARTt12 = A.map_ret_curr(tstp,t)(a1,a2);
        cplx BRTt12 = B.map_ret_curr(tstp,t)(b1,b2);

        C.map_les_curr(t,tstp)(c1,c2) = cplxi * ALtT12 * BGTt21;
        C.map_ret_curr(tstp,t)(c1,c2) = cplxi * ARTt12 * BLtT21 - cplxi * std::conj(ALtT21) * std::conj(BRTt12);
      }
    }
  }

  void Bubble2(int tstp, herm_matrix_hodlr &C, int c1, int c2, herm_matrix_hodlr &A, int a1, int a2, herm_matrix_hodlr &B, int b1, int b2) {
    if(tstp == -1) {
      // C_{c1,c2}(tau) = - A_{a1,a2}(tau) * B_{b1,b2}(tau)

      int ntau = B.ntau();
      for(int tau = 0; tau < ntau; tau++) {
        C.map_mat(tau)(c1,c2) = - A.map_mat(tau)(a1,a2) * B.map_mat(tau)(b1,b2);
      }
    }
    else {
      cplx cplxi(0.,1.);

      // C^\rceil_{c1,c2}(t,tau) = i A^\rceil_{a1,a2}(t,tau) * B^\rceil_{b1,b2}(t,tau)
      int ntau = C.ntau();
      for(int tau = 0; tau < ntau; tau++) {
        C.map_tv(tstp,tau)(c1,c2) = cplxi * A.map_tv(tstp,tau)(a1,a2) * B.map_tv(tstp,tau)(b1,b2);
      }

      for(int t = 0; t <= tstp; t++) {
        cplx ALtT12 = A.map_les_curr(t,tstp)(a1,a2);
        cplx BLtT12 = B.map_les_curr(t,tstp)(b1,b2);
        
        cplx ARTt12 = A.map_ret_curr(tstp,t)(a1,a2);
        cplx BRTt12 = B.map_ret_curr(tstp,t)(b1,b2);
        cplx ALtT21 = A.map_les_curr(t,tstp)(a2,a1);
        cplx BLtT21 = B.map_les_curr(t,tstp)(b2,b1);

        C.map_les_curr(t,tstp)(c1,c2) = cplxi * ALtT12 * BLtT12;
        C.map_ret_curr(tstp,t)(c1,c2) = cplxi * ARTt12 * BRTt12 - cplxi * std::conj(ALtT21) * BRTt12 - cplxi * ARTt12 * std::conj(BLtT21);
      }
    }
  }

}//namespace
#endif
