#include <iostream>
#include <random>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::identity_matrix<double> ubimatrix;

class monte_carlo {

  public:

  monte_carlo() {
    traj_length = 100;
    n_accept = 0;
    n_reject = 0;
  }

  virtual ~monte_carlo() {}

  // Numver of iterations
  size_t n_iters;

  // Number of parameters
  size_t n_params;

	// Positions and momenta for current and next points
  ubvector pos_curr, pos_next;
  ubvector mom_curr, mom_next;

  // Weights of current and next points
  double wgt_curr, wgt_next;
  
  // Potential and kinetic energies
  double pot_curr, pot_next;
  double kin_curr, kin_next;
  
  // Gradient of potential energy 
  ubvector grad_pot;
  
  // Step sizes for each parameter
  ubvector stepsize;

  // Trajectory length
  size_t traj_length;

  // Number of acceptance
  int n_accept;

  // Number of rejection
  int n_reject;

  // HMC base function
  int hmc(std::function<double(ubvector)>, ubvector);

};


class prob_density {
  
  public:
  
  prob_density() {
    mean_x = 0.0;
    mean_y = 0.0;
    width_x = 1.0;
    width_y = 1.0;
    rho = 0.0;
  }
  
  virtual ~prob_density() {}

  // Center and width of the distribution
  double mean_x, mean_y; 
  double width_x, width_y;
  
  // Correlation coefficient
  double rho;
  
  // Bivariate normal distribution function
  double binorm_pdf(ubvector);

};


class gradient {

  public:

  gradient() {
    h_rel = 1.0e-6;
    h_min = 1.0e-15;
  }

  virtual ~gradient() {}

  // Relative and minimum step sizes
  double h_rel, h_min;

  void grad_potential(std::function<double(ubvector)>, ubvector, ubvector);

};
