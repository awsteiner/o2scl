/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2026, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_INTERPM_LIBTORCH_H
#define O2SCL_INTERPM_LIBTORCH_H

#if defined(O2SCL_SET_LIBTORCH) || defined(DOXYGEN)

#include <torch/torch.h>
#include <o2scl/interpm_base.h>
#include <o2scl/tensor.h>
#include <o2scl/set_libtorch.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/lib_settings.h>
#include <o2scl/set_openmp.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif

/** \file interpm_libtorch.h
 */

namespace o2scl {
  
  /** \brief A simple MLP network with hidden layers
   */
  struct torch_mlp : torch::nn::Module {

    /// Hidden layers
    torch::nn::ModuleList hidden_layers{nullptr};
  
    /// Output layer
    torch::nn::Linear output_layer{nullptr};

    /// Number of input dimensions (set in constructor)
    size_t input_dim;
  
    /// Number of output dimensions (set in constructor)
    size_t output_dim;

    /// Activation function (default "silu")
    std::string act_func;
    
    /** \brief Create an empty network
     */
    torch_mlp() {
      input_dim=0;
      output_dim=0;
      act_func="silu";
    }
    
    /** \brief Create an MLP network
     */
    torch_mlp(size_t in_dim, const std::vector<size_t> &hidden_dims,
              size_t out_dim) {
      act_func="silu";
      init(in_dim,hidden_dims,out_dim);
    }
    
    /** \brief Create an MLP network
     */
    void init(size_t in_dim, const std::vector<size_t> &hidden_dims,
              size_t out_dim) {

      input_dim=in_dim;
      output_dim=out_dim;
    
      hidden_layers=register_module("hidden_layers",
                                    torch::nn::ModuleList());
      
      int prev_dim=(int)input_dim;

      // Create hidden layers
      for (size_t i=0; i<hidden_dims.size(); ++i) {
        torch::nn::Linear layer=torch::nn::Linear(prev_dim,
                                                  (int)hidden_dims[i]);
        hidden_layers->push_back(layer);
        prev_dim=(int)(hidden_dims[i]);
      }

      // Output layer
      output_layer=register_module("output_layer",
                                   torch::nn::Linear(prev_dim,
                                                     (int)output_dim));
      
      return;
    }

    /** \brief Forward propagation
     */
    torch::Tensor forward(torch::Tensor x) {
    
      for (std::shared_ptr<torch::nn::Module> &layer_ptr : *hidden_layers) {
        torch::nn::Linear layer=
	  std::dynamic_pointer_cast<torch::nn::LinearImpl>(layer_ptr);
        if (act_func=="tanh") {
          x=torch::tanh(layer->forward(x));
        } else if (act_func=="relu") {
          x=torch::nn::functional::relu(layer->forward(x));
        } else if (act_func=="silu") {
          x=torch::nn::functional::silu(layer->forward(x));
        } else if (act_func=="gelu") {
          x=torch::nn::functional::gelu(layer->forward(x));
        } else if (act_func=="softplus") {
          x=torch::nn::functional::softplus(layer->forward(x));
        } else if (act_func=="leaky_relu") {
          x=torch::nn::functional::leaky_relu(layer->forward(x));
        } else {
          O2SCL_ERR("Invalid activation.",o2scl::exc_einval);
        }
      }
    
      return output_layer->forward(x);
    }
  
  };

  /** \brief Multidimensional interpolation with libtorch (experimental)

      \verbatim embed:rst
      See also the :ref:`Higher-dimensional Interpolation` 
      section of the User's guide. 
      \endverbatim
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_x_t=o2scl::const_matrix_view_table<>,
           class mat_y_t=o2scl::matrix_view_table<>>
  class interpm_libtorch :
    public interpm_base<vec_t,mat_x_t,mat_y_t> {

  protected:
    
    /// MLP network
    std::shared_ptr<torch_mlp> model;

    /// Input dimensions
    size_t in_dim;

    /// Output dimensions
    size_t out_dim;
    
  public:

    interpm_libtorch() : device(torch::kCPU) {
      epochs=1000;
      patience=500;
      hidden_size={64,64,64};
      test_size=0.0;
      verbose=1;
      adam_lr=1.0e-3;
      seed=0;
      epoch_step=50;
    }
    
    /// Number of epochs (default 1000)
    int epochs;

    /// Patience (default 500)
    int patience;

    /// The most recent training loss
    double train_loss;

    /// The most recent validation loss
    double val_loss;

    /// Fraction of data to use for validation (default 0.0)
    double test_size;
    
    /// Size of the hidden layers (default {64,64,64})
    std::vector<size_t> hidden_size;

    /// Verbosity (default 1)
    int verbose;

    /// Adam learning rate (default 1.0e-3)
    double adam_lr;

    /// Device
    torch::Device device;

    /// The last device used
    std::string dev_str;

    /** \brief The random seed used to initialize the model's
        weights (and, if \ref test_size is nonzero, the
        train/validation split), or 0 for an automatic,
        thread-unique seed (default 0)

        If nonzero, \ref set_data_tensor() calls
        <tt>torch::manual_seed(seed)</tt> before constructing the
        model, so every call with the same nonzero \ref seed value
        starts from the same initial weights. If left at the
        default of 0, each call instead gets its own automatically
        generated seed via \ref o2scl::rng_set_seed(), which is
        guaranteed unique across concurrent OpenMP threads (unlike
        LibTorch's own unseeded ambient state, which is simply
        whatever the single global generator's state happens to be
        at that moment) -- see the comments in \ref
        set_data_tensor() for why this matters when this class is
        used from multiple threads at once, e.g. via \ref
        o2scl::nucmass_fit_iso::fit_interp_pass() with its
        n_threads greater than 1.
    */
    int seed;

    /// Number of steps between sparse verbose output (default 50)
    int epoch_step;
    
    /** \brief Set the data to be interpolated
     */
    virtual int set_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &user_x, mat_y_t &user_y) {
      
      o2scl::tensor<> tin, tout;
      std::vector<size_t> in_size={n_pts,n_in}, out_size={n_pts,n_out};
      tin.resize(2,in_size);
      tout.resize(2,out_size);
      for(size_t j=0;j<n_pts;j++) {
        std::vector<size_t> ix;
        for(size_t k=0;k<n_in;k++) {
          ix={j,k};
          tin.get(ix)=user_x(j,k);
        }
        for(size_t k=0;k<n_out;k++) {
          ix={j,k};
          tout.get(ix)=user_y(j,k);
        }
      }

      return set_data_tensor(n_in,n_out,n_pts,tin,tout);
    }

    /** \brief Set the data to be interpolated (tensor form)
     */
    int set_data_tensor(size_t n_in, size_t n_out, size_t n_pts,
                        const o2scl::tensor<> &params,
                        const o2scl::tensor<> &outputs) {

      device=(torch::cuda::is_available() ? torch::kCUDA : torch::kCPU);
              
      in_dim=n_in;
      out_dim=n_out;

      std::ostringstream oss;
      oss << device;
      dev_str=oss.str();
      if (verbose>0) {
        std::cout << "interpm_libtorch::set_data_tensor(): "
                  << "Device: " << dev_str << std::endl;
      }
      
      // ────────────────────────────────────────────────────────────
      // Convert input o2scl tensor objects to torch tensor objects
      // (shouldn't require a copy except for the move to the GPU)
      
      void *vpp=(void *)(params.get_data().data());
      torch::Tensor tp=torch::from_blob(vpp,{((long int)n_pts),
					     ((long int)n_in)},
	torch::kFloat64).to(device);
      void *vpo=(void *)(outputs.get_data().data());
      torch::Tensor to=torch::from_blob(vpo,{((long int)n_pts),
					     ((long int)n_out)},
	torch::kFloat64).to(device);

      // ────────────────────────────────────────────────────────────
      // Seed LibTorch's global random number generator, then
      // construct (and thus randomly initialize) the model and, if
      // used, draw the train/validation split permutation. These are
      // the only steps in this function that consume LibTorch's
      // global (not thread-local) default generator -- weight
      // initialization happens inside the \ref torch_mlp constructor
      // via that ambient generator, and torch::randperm() below
      // draws from the same one. If this is called concurrently from
      // multiple OpenMP threads (e.g. via \ref
      // o2scl::nucmass_fit_iso::fit_interp_pass() with its
      // n_threads greater than 1), reseeding and then drawing from
      // that shared generator must happen atomically per call, or
      // one thread's reseed could land in the middle of another
      // thread's draws and silently corrupt both. The training loop
      // below draws no further randomness, so only this setup
      // section needs to be serialized -- the (expensive) per-epoch
      // work stays fully parallel. This mirrors the
      // OpenMP-critical-block pattern used by \ref
      // o2scl::rng_set_seed(), which is also used here (when \ref
      // seed is left at its default of 0) to give each call its own
      // automatic seed unique across threads, rather than every
      // thread drawing from the same unseeded ambient state.

      int total_samples=n_pts;
      int train_size=static_cast<int>(test_size*total_samples);
      torch::Tensor perm;

#ifdef O2SCL_SET_OPENMP
#pragma omp critical (o2scl_interpm_libtorch_seed)
#endif
      {
        if (seed!=0) {
          torch::manual_seed(seed);
        } else {
          o2scl::rng<> r;
          o2scl::rng_set_seed(r);
          torch::manual_seed(r.get_seed());
        }

        model=std::make_shared<torch_mlp>(n_in,hidden_size,n_out);

        if (test_size>0.0) {
          perm=torch::randperm(total_samples,torch::kLong).to(device);
        }
      }

      model->to(device,torch::kFloat64);

      torch::optim::Adam optimizer(model->parameters(),
                                   torch::optim::AdamOptions(adam_lr));

      //const bool enforce_exact=false;
      //const double lambda_exact=1000.0;

      // ────────────────────────────────────────────────────────────
      // Train/validation split (index selection only -- no further
      // random draws, so this doesn't need the critical block above)

      torch::Tensor train_inputs, train_targets, val_inputs, val_targets;

      if (test_size>0.0) {

        torch::Tensor train_idx=perm.index
          ({torch::indexing::Slice(0,train_size)}).to(device);
        torch::Tensor val_idx=perm.index({torch::indexing::Slice
            (train_size,torch::indexing::None)}).to(device);

        train_inputs=tp.index_select(0,train_idx).to(device);
        train_targets=to.index_select(0,train_idx).to(device);

        val_inputs=tp.index_select(0,val_idx).to(device);
        val_targets=to.index_select(0,val_idx).to(device);

      } else {

        train_inputs=tp.to(device);
        train_targets=to.to(device);

      }

      // ────────────────────────────────────────────────────────────
      // Early stopping parameters

      int patience_counter=0;

      double best_loss=std::numeric_limits<double>::infinity();

      // Best model parameters
      std::vector<torch::Tensor> best_params;

      // ────────────────────────────────────────────────────────────
      // Training

      for (int epoch=0;epoch<epochs;epoch++) {

        model->train();
        optimizer.zero_grad();

        torch::Tensor pred=model->forward(train_inputs);
        torch::Tensor loss=torch::mse_loss(pred,train_targets);

        //if (enforce_exact) {
	//torch::Tensor pred_interp=model->forward(interp_nodes);
	//torch::Tensor loss_interp=torch::mse_loss(pred_interp,
        //interp_targets);
	//loss=loss+lambda_exact*loss_interp;
	//}

        loss.backward();
        optimizer.step();

        // ────────────────────────────────────────────────────────────
        // Validation

        model->eval();

        if (test_size==0.0) {
          
          train_loss=loss.item<double>();
          
        } else {
          
          torch::NoGradGuard no_grad;
          
          torch::Tensor val_pred=model->forward(val_inputs);
          torch::Tensor val_loss_tensor=torch::mse_loss
            (val_pred,val_targets);
          
          val_loss=val_loss_tensor.item<double>();
          train_loss=loss.item<double>();

        }

        // ────────────────────────────────────────────────────────────
        // Early stopping logic

        if ((test_size>0.0 && val_loss<best_loss) ||
            (test_size==0.0 && train_loss<best_loss)) {

          if (test_size>0.0) {
            best_loss=val_loss;
          } else {
            best_loss=train_loss;
          }

          patience_counter=0;

          // Save best parameters
          best_params.clear();
      
          std::vector<torch::Tensor> params=model->parameters();
          for (size_t i = 0; i < params.size(); ++i) {
            best_params.push_back(params[i].detach().clone());
          }
      
        } else {
	  
          patience_counter++;
	  
        }

        if (verbose>0) {
          if (epoch%epoch_step==0 || verbose>1) {
            std::cout << "interpm_libtorch::set_data_tensor(): ";
            std::cout << "Epoch ";
            std::cout.width(4);
            std::cout << epoch
                      << " Train loss: " << train_loss;
            if (test_size>0.0) {
              std::cout << " Val loss: "   << val_loss;
            }
            std::cout << std::endl;
          }
        }
	
        if (patience_counter >= patience) {
          break;
        }
      }
      
      // ────────────────────────────────────────────────────────────
      // Restore best model
      
      for (size_t i=0; i<best_params.size(); ++i) {
        model->parameters()[i].data().copy_(best_params[i]);
      }

      if (verbose>0) {
        if (test_size>0.0) {
          std::cout << "Best validation loss: "
                    << best_loss << std::endl;
        } else {
          std::cout << "Best loss: "
                    << best_loss << std::endl;
        }
      }
      
      return 0;
    }

    /** \brief Evaluate the interpolation at point \c x,
        returning \c y
    */
    virtual int eval(const vec_t &x, vec_t &y) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                  "interpm_libtorch::eval().",o2scl::exc_einval);
      }

      model->eval();
      torch::NoGradGuard no_grad;
 
      torch::TensorOptions opts=torch::TensorOptions().dtype
        (torch::kFloat64).device(torch::kCPU);
      torch::Tensor xt=torch::empty({(long int)1,(long int)in_dim},opts);
          
      double *xt_ptr=xt.data_ptr<double>();
      for(size_t i=0;i<x.size();i++) {
        xt_ptr[i]=x[i];
      }

      // Move to the GPU (if necessary), evaluate the model,
      // then move back
      torch::Tensor xt_dev=xt.to(device);
      torch::Tensor pred=model->forward(xt_dev);
      torch::Tensor pred_dev=pred.to(torch::kCPU);
      double *pred_ptr=pred_dev.data_ptr<double>();

      for(size_t j=0;j<out_dim;j++) {
        y[j]=pred_ptr[j];
      }

      return 0;
    }

    /** \brief Evaluate the interpolation at a list of points in \c x,
        returning a list of results in \c y (tensor form)
    */
    virtual int eval_list_tensor(const o2scl::tensor<> &x,
                                 o2scl::tensor<> &y) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                  "interpm_libtorch::eval_list_tensor().",
                   o2scl::exc_einval);
      }

      model->eval();
      torch::NoGradGuard no_grad;

      size_t n_pts=x.get_size(0);
      
      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt=torch::from_blob(vpp,{((long int)n_pts),
					     ((long int)in_dim)},
	torch::kFloat64).to(device);

      // Evaluate the model, then move back to the CPU

      torch::Tensor pred=model->forward(xt);
      torch::Tensor pred_dev=pred.to(torch::kCPU);
      torch::TensorAccessor<double,2> pred_acc=
        pred_dev.accessor<double,2>();
      
      std::vector<size_t> ix(2);
      for(size_t j=0;j<out_dim;j++) {
        ix[1]=j;
        for(size_t k=0;k<n_pts;k++) {
          ix[0]=k;
          y.get(ix)=pred_acc[k][j];
        }
      }

      return 0;
    }

    /** \brief Set the data to be interpolated
     */
    virtual int eval_list(size_t n_pts,
                          mat_x_t &user_x, mat_y_t &user_y) {
      
      o2scl::tensor<> tin, tout;
      
      std::vector<size_t> in_size={n_pts,in_dim};
      std::vector<size_t> out_size={n_pts,out_dim};
      
      tin.resize(2,in_size);
      tout.resize(2,out_size);
      
      for(size_t j=0;j<n_pts;j++) {
        std::vector<size_t> ix;
        for(size_t k=0;k<in_dim;k++) {
          ix={j,k};
          tin.get(ix)=user_x(j,k);
        }
      }

      int ret=eval_list_tensor(tin,tout);
      
      for(size_t j=0;j<n_pts;j++) {
        std::vector<size_t> ix;
        for(size_t k=0;k<out_dim;k++) {
          ix={j,k};
          tout.get(ix)=user_y(j,k);
        }
      }

      return ret;
    }

    /** \brief Evaluate the derivative with respect to variable \c ix
        at point \c x, returning \c y
    */
    virtual int deriv(const vec_t &x, vec_t &y, size_t ix) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_libtorch::deriv().",o2scl::exc_einval);
      }

      model->eval();

      torch::TensorOptions opts=torch::TensorOptions().dtype
        (torch::kFloat64).device(torch::kCPU);
      
      torch::Tensor xt=torch::empty({(long int)1,(long int)in_dim},opts);

      double *xt_ptr=xt.data_ptr<double>();
      for(size_t i=0;i<x.size();i++) {
        xt_ptr[i]=x[i];
      }

      // Move to the GPU (if necessary), evaluate the model,
      // then move back
      torch::Tensor xt_dev=xt.clone().to(device);
      xt_dev.requires_grad_(true);
      torch::Tensor pred=model->forward(xt_dev);

      // We have to compute the derivative separately, one call to
      // autograd::grad() for each output dimension
      for(size_t j=0;j<out_dim;j++) {
        
        // Fourth argument is "retain_graph", and fifth is
        // "create_graph". 
        torch::Tensor grad_out=torch::zeros({((long int)1),
            ((long int)out_dim)},torch::kFloat64).to(device);
        grad_out[0][j]=1.0;
        
        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt_dev},{grad_out},true,false)[0].to(torch::kCPU);
        double *grads_ptr=grads.data_ptr<double>();
        
        y[j]=grads_ptr[ix];
      }

      return 0;
    }
    
    /** \brief Desc
     */
    virtual int deriv_list_tensor
    (const o2scl::tensor<> &x, o2scl::tensor<> &y, size_t ix) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_libtorch::deriv_list_tensor().",
                   o2scl::exc_einval);
      }

      model->eval();

      size_t n_pts=x.get_size(0);
      
      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt_cpu=torch::from_blob
        (vpp,{((long int)n_pts),((long int)in_dim)},
         torch::kFloat64);
      
      torch::Tensor xt=xt_cpu.clone().to(device).requires_grad_(true);
      torch::Tensor pred=model->forward(xt);

      // We have to compute the derivative separately, one call to
      // autograd::grad() for each output dimension
      for(size_t j=0;j<out_dim;j++) {

        // Fourth argument is "retain_graph", and fifth is
        // "create_graph". 
        torch::Tensor grad_out=torch::zeros({((long int)n_pts),
            ((long int)out_dim)},torch::kFloat64).to(device);
        for(size_t k=0;k<n_pts;k++) {
          grad_out[k][j]=1.0;
        }
        
        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt},{grad_out},true,false)[0].to(torch::kCPU);
        
        torch::TensorAccessor<double,2> grads_acc=
          grads.accessor<double,2>();
        
        std::vector<size_t> index(2);
        index[1]=j;
        for(size_t k=0;k<n_pts;k++) {
          index[0]=k;
          y.get(index)=grads_acc[k][ix];
        }
        
      }

      return 0;
    }

    /** \brief Compute the second derivative of output j with respect
        to inputs ix and ix2 at point x, storing results in y
    */
    virtual int deriv2(const vec_t &x, vec_t &y,
                       size_t ix, size_t ix2) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_libtorch::deriv2().",o2scl::exc_einval);
      }
      
      model->eval();
      
      torch::TensorOptions opts=torch::TensorOptions().dtype
        (torch::kFloat64).device(torch::kCPU);
      
      torch::Tensor xt=torch::empty({(long int)1,(long int)in_dim},opts);
      
      double *xt_ptr=xt.data_ptr<double>();
      for(size_t i=0;i<x.size();i++) {
        xt_ptr[i]=x[i];
      }
      
      // Move to the GPU (if necessary)
      torch::Tensor xt_dev=xt.clone().to(device);
      xt_dev.requires_grad_(true);
      torch::Tensor pred=model->forward(xt_dev);
      
      // One call to autograd::grad() for each output dimension
      for(size_t j=0;j<out_dim;j++) {
        
        // First derivative, use create_graph=true so we can
        // differentiate again
        torch::Tensor grad_out=torch::zeros({((long int)1),
            ((long int)out_dim)},torch::kFloat64).to(device);
        grad_out[0][j]=1.0;
        
        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt_dev},{grad_out},true,true )[0];
        
        // Derivative with respect to ix has shape [1,in_dim], 
        // and grads_ix has shape [1,1].
        torch::Tensor grads_ix=grads.index
          ({torch::indexing::Slice(),((long int)ix)}).unsqueeze(1);
        
        // Second derivative
        torch::Tensor grad_out2=torch::ones_like(grads_ix);
        torch::Tensor grads2=torch::autograd::grad
          ({grads_ix},{xt_dev},{grad_out2},true,false
           )[0].to(torch::kCPU);
        
        double *grads2_ptr=grads2.data_ptr<double>();
        
        // Variable grads2 has shape [1,in_dim], extract component ix2
        y[j]=grads2_ptr[ix2];
      }
      
      return 0;
    }    
    
    /** \brief Compute the second derivative of outputs with respect to
        inputs ix and ix2 at a list of points x, storing results in y
    */
    virtual int deriv2_list_tensor
    (const o2scl::tensor<> &x, o2scl::tensor<> &y,
     size_t ix, size_t ix2) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_libtorch::deriv2_list_tensor().",
                   o2scl::exc_einval);
      }
      
      model->eval();
      
      size_t n_pts=x.get_size(0);
      
      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt_cpu=torch::from_blob
        (vpp,{((long int)n_pts),((long int)in_dim)},
         torch::kFloat64);
      
      torch::Tensor xt=xt_cpu.clone().to(device).requires_grad_(true);
      torch::Tensor pred=model->forward(xt);
      
      // One call to autograd::grad() for each output dimension
      for(size_t j=0;j<out_dim;j++) {
        
        // First derivative, use create_graph=true so we can
        // differentiate again
        torch::Tensor grad_out=torch::zeros({((long int)n_pts),
            ((long int)out_dim)},torch::kFloat64).to(device);
        for(size_t k=0;k<n_pts;k++) {
          grad_out[k][j]=1.0;
        }
        
        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt},{grad_out},true,true)[0];
        
        // Derivative with respect to ix has shape [1,in_dim], 
        // and grads_ix has shape [1,1].
        torch::Tensor grads_ix=grads.index
          ({torch::indexing::Slice(),((long int)ix)}).unsqueeze(1);
        
        // Second derivative
        torch::Tensor grad_out2=torch::ones_like(grads_ix);
        
        torch::Tensor grads2=torch::autograd::grad
          ({grads_ix},{xt},{grad_out2},true,false)[0].to(torch::kCPU);
        
        // Variable grads2 has shape [1,in_dim], extract component ix2
        torch::TensorAccessor<double,2> grads2_acc=
          grads2.accessor<double,2>();
        
        std::vector<size_t> index(2);
        index[1]=j;
        for(size_t k=0;k<n_pts;k++) {
          index[0]=k;
          y.get(index)=grads2_acc[k][ix2];
        }
      }
      
      return 0;
    }
    
    /** \brief Save the interpolator to an HDF5 file
     */
    int hdf_output(o2scl_hdf::hdf_file &hf, std::string name) const {
      
      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_libtorch::save().",
                   o2scl::exc_einval);
      }

      // Start group
      hid_t top=hf.get_current_id();
      hid_t group=hf.open_group(name);
      hf.set_current_id(group);

      // Add typename
      hf.sets_fixed("o2scl_type","interpm_libtorch");
      
      // Save scalar metadata
      hf.set_szt("in_dim",in_dim);
      hf.set_szt("out_dim",out_dim);
      hf.sets("act_func",model->act_func);
      
      // Save hidden layer sizes
      std::vector<size_t> hs(hidden_size.begin(),hidden_size.end());
      hf.set_szt_vec("hidden_sizes",hs);
      
      // Save each named parameter as a flat double vector
      // plus its shape as an int vector
      for (const auto &entry : model->named_parameters()) {

        std::string pname=entry.key();
        torch::Tensor t=entry.value().detach()
          .to(torch::kCPU).to(torch::kFloat64).contiguous();

        std::vector<int> shape;
        for (int i=0;i<t.dim();i++) {
          shape.push_back((int)t.size(i));
        }
        hf.seti_vec(pname+"_shape",shape);

        std::vector<double> data(t.data_ptr<double>(),
                                 t.data_ptr<double>()+t.numel());
        hf.setd_vec(pname+"_data",data);
      }

      // Save buffers
      for (const auto &entry : model->named_buffers()) {

        std::string pname=entry.key();
        torch::Tensor t=entry.value().detach()
          .to(torch::kCPU).to(torch::kFloat64).contiguous();
        
        std::vector<int> shape;
        for (int i=0;i<t.dim();i++) {
          shape.push_back((int)t.size(i));
        }
        hf.seti_vec(pname+"_buf_shape",shape);
        
        std::vector<double> data(t.data_ptr<double>(),
                                 t.data_ptr<double>()+t.numel());
        hf.setd_vec(pname+"_buf_data",data);
      }

      // Close group
      hf.close_group(group);

      // Return location to previous value
      hf.set_current_id(top);
      
      return 0;
    }

    /** \brief Load the interpolator with a given name from an HDF5 file
     */
    void hdf_input(o2scl_hdf::hdf_file &hf, std::string name) {
      hdf_input_n(hf,name);
      return;
    }
    
    /** \brief Load the first interpolator from an HDF5 file
        and return the name
     */
    int hdf_input_n(o2scl_hdf::hdf_file &hf, std::string &name) {

      // If no name specified, find name of first group of specified type
      if (name.length()==0) {
        hf.find_object_by_type("interpm_libtorch",name);
        if (name.length()==0) {
          O2SCL_ERR3("No object of type interpm_libtorch found in ",
                     "o2scl_hdf::hdf_input_n(hdf_file &,",
                     "interpm_libtorch &,string &).",exc_efailed);
        }
      }

      // Open main group
      hid_t top=hf.get_current_id();
      hid_t group=hf.open_group(name);
      hf.set_current_id(group);

      // Check typename
      std::string type2;
      hf.gets_fixed("o2scl_type",type2);
      if (type2!="interpm_libtorch") {
        O2SCL_ERR2("Typename in HDF group does not match ",
                   "class in table::hdf_input().",exc_einval);
      }
      
      // Load scalar metadata
      std::string act_func_str;
      hf.get_szt("in_dim",in_dim);
      hf.get_szt("out_dim",out_dim);
      hf.gets("act_func",act_func_str);
      
      // Load hidden layer sizes (key name must match hdf_output()'s
      // "hidden_sizes", not the unrelated per-instance "hidden_size"
      // member name)
      std::vector<size_t> hs;
      hf.get_szt_vec("hidden_sizes",hs);
      hidden_size=hs;
      
      // Reconstruct model architecture before loading weights
      model=std::make_shared<torch_mlp>(in_dim,hidden_size,out_dim);
      model->act_func=act_func_str;
      model->to(device,torch::kFloat64);
      
      // Load each named parameter
      for (const auto &entry : model->named_parameters()) {

        std::string pname=entry.key();

        std::vector<int> shape;
        hf.geti_vec(pname+"_shape",shape);

        std::vector<double> data;
        hf.getd_vec(pname+"_data",data);

        std::vector<long int> tshape(shape.begin(),shape.end());
        torch::Tensor t=torch::from_blob
          (data.data(),
           torch::IntArrayRef(tshape.data(),(size_t)tshape.size()),
           torch::kFloat64).clone().to(device);

        torch::NoGradGuard no_grad;
        entry.value().copy_(t);
      }

      // Load buffers
      for (const auto &entry : model->named_buffers()) {

        std::string pname=entry.key();

        std::vector<int> shape;
        hf.geti_vec(pname+"_buf_shape",shape);

        std::vector<double> data;
        hf.getd_vec(pname+"_buf_data",data);

        std::vector<long int> tshape(shape.begin(),shape.end());
        torch::Tensor t=torch::from_blob
          (data.data(),
           torch::IntArrayRef(tshape.data(),(size_t)tshape.size()),
           torch::kFloat64).clone().to(device);

        torch::NoGradGuard no_grad;
        entry.value().copy_(t);
      }
      
      hf.close_group(group);
      
      // Return location to previous value
      hf.set_current_id(top);
      
      return 0;
    }
    
  };
  
}

#endif

#endif



