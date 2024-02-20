#include <fstream>
#include <cmath>
#include <chrono>
#include <functional>
#include "hmc.h"

using namespace std;

double prob_density::binorm_pdf(ubvector u) {
  double shift_x = u(0)-mean_x;
  double shift_y = u(1)-mean_y;
  double x = pow(shift_x/width_x, 2.0);
  double y = pow(shift_y/width_y, 2.0);
  double z = 2.0*rho*shift_x*shift_y/(width_x*width_y);
  double r = 1.0-rho*rho;
  double c = 2.0*M_PI*width_x*width_y*sqrt(r);
  
  return 1.0/c*exp(-0.5*(x+y-z)/r);
}


void gradient::grad_potential(std::function<double(ubvector)> f, ubvector x, ubvector g) {
  double fv1, fv2, h;
  
  fv1=f(x);
  
  for(size_t i=0; i<x.size(); i++) {
	  
    h=h_rel*fabs(x(i));
	  if (fabs(h)<=h_min) h=h_rel;
	  
    x(i)+=h;
	  fv2=f(x);
	  x(i)-=h;
	  g(i)=(fv2-fv1)/h;
  }
}


int monte_carlo::hmc(std::function<double(ubvector)> f, ubvector x) {

  std::random_device seed;
  std::mt19937 gen(seed());
  std::uniform_real_distribution<> unif(0, 1);
  std::normal_distribution<> norm(0, 1);

  pos_curr.resize(n_params);
  mom_curr.resize(n_params);
  pos_next.resize(n_params);
  mom_next.resize(n_params);
  grad_pot.resize(n_params);

  pos_curr=x;

  ubmatrix mass_inv=ubimatrix(n_params);
  mass_inv*=0.1;

  bool done=false;
  size_t i_iters=0;

  ofstream file;
  file.open("binormal.dat");

  // Get current time at start
  auto start_time = chrono::high_resolution_clock::now();

  while (!done) {

    // Print status report
    cout << "Iteration " << i_iters << " of " << n_iters << ": "
    << "n_accept=" << n_accept << ", n_reject=" << n_reject << endl;

    // Get current weight
    wgt_curr=f(pos_curr);

    // Write current positions and weight to file
    file << std::scientific;
    for (size_t i=0; i<n_params; i++) {
      file << pos_curr(i) << "\t ";
    }
    file << wgt_curr << endl;

    // Set current momentum
    for (size_t i=0; i<n_params; i++) {
      mom_curr(i) = norm(gen);
    }

    // Rescale momentum by stepsize
    mom_curr=element_prod(mom_curr, stepsize);

    // Set next position and momentum
    pos_next=pos_curr;
    mom_next=mom_curr;

    // Compute gradient of potential energy
    gradient g;
    g.grad_potential(f, pos_next, grad_pot);

    // Make a half step for momentum at the beginning
    mom_next-=0.5*element_prod(stepsize, grad_pot);

    // Leapfrog updates: Full steps for position and momentum
    for (size_t i=1; i<=traj_length; i++) {

      // Make a full step for the position
      pos_next+=element_prod(stepsize, mom_next);

      // Compute gradient of potential energy
      g.grad_potential(f, pos_next, grad_pot);

      // Make a full step for the momentum except at the end
      if (i!=traj_length) {
        mom_next-=element_prod(stepsize, grad_pot);
      }
    }

    // Make a half step for momentum at the end
    mom_next-=0.5*element_prod(stepsize, grad_pot);

    // Negate momentum to make the proposal symmetric
    mom_next=-mom_next;

    // Evaluate potential energies for current and next points
    pot_curr=-log(0.5*f(pos_curr));
    pot_next=-log(0.5*f(pos_next));

    // Evaluate kinetic energies for current and next points
    ubvector mom_curr_t=trans(mom_curr);
    ubvector mom_next_t=trans(mom_next);
    kin_curr=0.5*inner_prod(mom_curr, prod(mass_inv, mom_curr_t));
    kin_next=0.5*inner_prod(mom_next, prod(mass_inv, mom_next_t));

    // Accept or reject the proposal
    double r=unif(gen);
    
    if (r<exp(pot_curr-pot_next+kin_curr-kin_next)) {

      // Accept: Update current position
      n_accept++;
      pos_curr=pos_next;
      
    } else {
      
      // Reject: Keep current position
      n_reject++;
    }

    // Check whether to terminate, else increment iteration counter
    if (i_iters==n_iters) {
      done=true;
      cout << "HMC terminated: Maximum iterations reached." << endl;
    } else {
      i_iters++;
    }
  }

  file.close();

  // Get current time at end
  auto end_time = chrono::high_resolution_clock::now();
  
  // Compute elapsed time and average time per iteration
  auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);
  auto avg_time = duration.count()/n_iters;

  // Print final status report
  cout << "n_accept=" << n_accept << ", n_reject=" << n_reject
       << ", accept_rate=" << (double)n_accept/(n_accept+n_reject) << endl;

  cout << "Time elapsed: " << duration.count() << " us" << endl;
  cout << "Time per iteration: " << avg_time << " us" << endl;

  return 0;
}


int main() {
  std::random_device seed;
  std::mt19937 gen(seed());
  std::uniform_real_distribution<> unif(-1, 1);

  monte_carlo mc;
  mc.n_params=2;
  mc.traj_length=20;
  mc.n_iters=10000;
  mc.stepsize.resize(mc.n_params);

  ubvector x(mc.n_params);
  for (size_t i=0; i<mc.n_params; i++) {
    mc.stepsize(i) = 0.18;
    x(i) = unif(gen);
  }

  prob_density pd;
  std::function<double(ubvector)> f = std::bind(&prob_density::binorm_pdf, 
                                      &pd, std::placeholders::_1);
  return mc.hmc(f, x);
}