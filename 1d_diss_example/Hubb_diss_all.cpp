#include "Hubb_diss_all.hpp"
#include <iostream>

// Constructor
Hubb_diss_all::Hubb_diss_all(function &U) 
    : U_(U) 
{}

void Hubb_diss_all::Sigma_2b(int tstp, 
                          herm_matrix_hodlr &G, 
                          herm_matrix_hodlr &Sigma) 
{
  int nao = G.size1();
  ZMatrix GL_tT(nao, nao);
  ZMatrix GL_Tt(nao, nao);
  ZMatrix GG_tT(nao, nao);
  ZMatrix GG_Tt(nao, nao);

  ZMatrix SG_tT(nao, nao);

  for(int t = 0; t <= tstp; t++) {
    SG_tT.setZero();
    cplx U2 = U_(tstp,0,0) * U_(t,0,0);
    GL_tT = G.map_les_curr(t,tstp);
    GL_Tt = -GL_tT.adjoint();
    GG_Tt = G.map_ret_curr(tstp,t) + GL_Tt;
    GG_tT = -GG_Tt.adjoint();

    Sigma.map_les_curr(t,tstp) =     U2 * GL_tT.cwiseProduct(GL_tT.cwiseProduct(GG_Tt.transpose()));

    SG_tT                      =    -U2 * GG_tT.cwiseProduct(GG_tT.cwiseProduct(GL_Tt.transpose()));
    SG_tT                     +=  2.*U2 * GL_tT.cwiseProduct(GL_tT.cwiseProduct(GL_Tt.transpose()));

    Sigma.map_ret_curr(tstp,t) = -SG_tT.adjoint() + Sigma.map_les_curr(t,tstp).adjoint();
  }

}

void Hubb_diss_all::Sigma_hf(int tstp, 
                         herm_matrix_hodlr &G,
                         function &ellG) 
{
  int nao = G.size1();
  ZMatrix rho(nao, nao);

  G.get_les_curr(tstp, tstp, rho);
  rho *= cplx(0.,-1.);

  ellG.get_map(tstp) = U_(tstp,0,0) * rho.diagonal().asDiagonal();
}
