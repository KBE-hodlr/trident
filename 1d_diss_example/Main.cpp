#include <filesystem>

#include "h_nessi/read_inputfile.hpp"
#include "h_nessi/dyson.hpp"
#include "Hubb_diss_all.hpp"

bool isPowerOfTwo(int n) {
    // A power of 2 must be strictly positive, and
    // the bitwise AND of n and (n - 1) must equal 0.
    return (n > 0) && ((n & (n - 1)) == 0);
}

using namespace h_nessi;

int main(int argc, char *argv[]) {

  int SolverOrder = 5;

  int Nt;
  int nlvl, r;
  double svdtol, h;

  int BootMaxIter, StepMaxIter;
  double BootMaxErr, StepMaxErr;
  int rho_version = 0;

  double mu = 0;
  double J = 1;
  double U_max;
  int size = 4;
  int xi = -1;

  char* flin;
  flin=argv[1];

  find_param(flin, "__nt=",     Nt);
  find_param(flin, "__nlevel=", nlvl);
  find_param(flin, "__svdtol=", svdtol);
  find_param(flin, "__h=",      h);

  find_param(flin, "__BootMaxIter=", BootMaxIter);
  find_param(flin, "__BootMaxErr=",  BootMaxErr);
  find_param(flin, "__StepMaxIter=", StepMaxIter);
  find_param(flin, "__StepMaxErr=",  StepMaxErr);

  find_param(flin, "__Umax=", U_max);

  //============================================================================
  //                           Initialize Objects
  //============================================================================

  function U(Nt, 1, 1);
  for(int t = 0; t < Nt; t++) U(t,0,0) = U_max * (std::erf(4*(t*h-2))+1)/2;

  dlr_info dlr(r, 0.9, 0.9, 1, size, xi);
  
  std::cout << "r=" << r << std::endl;

  Integration::Integrator I(SolverOrder);

  dyson dyson_sol(Nt, size, SolverOrder, dlr, rho_version, false);

  herm_matrix_hodlr G(Nt, r, nlvl, svdtol, size, size, xi, SolverOrder);
  herm_matrix_hodlr Sigma(Nt, r, nlvl, svdtol, size, size, xi, SolverOrder);

  int nblock = std::pow(2, nlvl)-1;
  for(int i = 0; i < nblock; i++) {
    std::cout << G.blkr2(i)-G.blkr1(i)+1 << " " << G.blkc2(i)-G.blkc1(i)+1 << std::endl;
    if(!isPowerOfTwo(G.blkr2(i)-G.blkr1(i)+1) or !isPowerOfTwo(G.blkc2(i)-G.blkc1(i)+1)) std::cout << "FALSE" << std::endl;
  }

  function H(Nt,size,size);
  function ellG(Nt,size,size);
  function ellL(Nt,size,size);

  for(int t = -1; t < Nt; t++) {
    for(int i = 0; i < size; i++) {
      H(t,i,(i+1)%size) = -J;
      H(t,i,(i-1+size)%size) = -J;
    }
  }

  // Coefficients for the |phi> orbital localized on the even sub-lattice (sites 0 and 2)
  // This explicitly assumes a lattice of at least N=4.
  int N = size;
  std::vector<double> a(N, 0.0);
  if (N >= 4) {
      a[0] = 1.0 / std::sqrt(2.0);
      a[2] = 1.0 / std::sqrt(2.0);
  }
  ZMatrix rho0(size,size);
  for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
          // Since a[i] are purely real, a_i * a_j^* is just a[i] * a[j]
          double val = a[i] * a[j];
          rho0(i,j) = cplx(val, 0.0);
      }
  }

  // G^< = -i\xi\rho
  G.map_les_curr(0,0) = cplx(0.,1.) * rho0;

  Hubb_diss_all se_eval(U);

  //============================================================================
  //                           Bootstrap
  //============================================================================

  for(int iter = 0; iter < BootMaxIter; iter++) {
    for(int tstp = 0; tstp <= SolverOrder; tstp++) {
      se_eval.Sigma_hf(tstp, G, ellG);
      se_eval.Sigma_2b(tstp, G, Sigma);
    }

    double err = dyson_sol.dyson_start_2leg_diss(G, mu, H, ellL, ellG, Sigma, I, h);

    std::cout  << "boot " << iter << " " << err << std::endl;

    if(err < BootMaxErr) break;
  }

  //============================================================================
  //                           Timestep
  //============================================================================

  for(int tstp = SolverOrder+1; tstp < Nt; tstp++) {
    G.update_blocks(I);
    Sigma.update_blocks(I);

    dyson_sol.extrapolate(G,I);

    for(int iter = 0; iter < StepMaxIter; iter++) {
      se_eval.Sigma_hf(tstp, G, ellG);
      se_eval.Sigma_2b(tstp, G, Sigma);

      double err = dyson_sol.dyson_timestep_2leg_diss(tstp, G, mu, H, ellL, ellG, Sigma, I, h);

      std::cout << tstp << " " << iter << " " << err << std::endl;

      if(err < StepMaxErr) break;
    }
  }

  h5e::File out_file(argv[2], h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
  G.write_to_hdf5(out_file, "G/");
  G.write_rho_to_hdf5(out_file, "/rho/", dlr);


  return 0;
}
