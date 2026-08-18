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

#ifdef O2SCL_LIBTORCH

#include <torch/torch.h>
#include <o2scl/interpm_base.h>
#include <o2scl/tensor.h>

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
    int input_dim;
  
    /// Number of output dimensions (set in constructor)
    int output_dim;

    /** \brief Create an empty network
     */
    torch_mlp() {
      input_dim=0;
      output_dim=0;
    }
    
    /** \brief Create an MLP network
     */
    torch_mlp(int in_dim, const std::vector<int> &hidden_dims,
              int out_dim) {
      init(in_dim,hidden_dims,out_dim);
    }
    
    /** \brief Create an MLP network
     */
    void init(int in_dim, const std::vector<int> &hidden_dims,
              int out_dim) {

      input_dim=in_dim;
      output_dim=out_dim;
    
      hidden_layers=register_module("hidden_layers",
                                    torch::nn::ModuleList());
      
      int prev_dim=input_dim;

      // Create hidden layers
      for (size_t i=0; i<hidden_dims.size(); ++i) {
        torch::nn::Linear layer=torch::nn::Linear(prev_dim,hidden_dims[i]);
        hidden_layers->push_back(layer);
        prev_dim=hidden_dims[i];
      }

      // Output layer
      output_layer=register_module("output_layer",
                                   torch::nn::Linear(prev_dim,output_dim));
      
      return;
    }

    /** \brief Forward propagation
     */
    torch::Tensor forward(torch::Tensor x) {
    
      for (std::shared_ptr<torch::nn::Module> &layer_ptr : *hidden_layers) {
        torch::nn::Linear layer=
	  std::dynamic_pointer_cast<torch::nn::LinearImpl>(layer_ptr);
        x=torch::nn::functional::silu(layer->forward(x));
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
    std::vector<int> hidden_size;

    /// Verbosity (default 1)
    int verbose;

    /// Adam learning rate (default 1.0e-3)
    double adam_lr;

    /// Device
    torch::Device device;

    /// The last device used
    std::string dev_str;

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

      device=(torch::cuda::is_available()
              ? torch::kCUDA : torch::kCPU);
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
      //torch::manual_seed(0);

      model=std::make_shared<torch_mlp>(n_in,hidden_size,n_out);
      
      model->to(device,torch::kFloat64);

      torch::optim::Adam optimizer(model->parameters(),
                                   torch::optim::AdamOptions(adam_lr));
  
      //const bool enforce_exact=false;
      //const double lambda_exact=1000.0;

      // ────────────────────────────────────────────────────────────
      // Train/validation split
      
      int total_samples=n_pts;

      torch::Tensor train_inputs, train_targets, val_inputs, val_targets;
      
      if (test_size>0.0) {
        
        int train_size=static_cast<int>(test_size*total_samples);
        
        torch::Tensor perm=torch::randperm
          (total_samples,torch::kLong).to(device);
        
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
          if (epoch%50==0 || verbose>1) {
            std::cout << "interpm_libtorch::set_data_tensor(): ";
            std::cout << "Epoch " << epoch
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

    /** \brief Desc
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
      torch::Tensor xt_dev=xt.to(device);
      xt_dev=xt_dev.detach().requires_grad_(true);
      torch::Tensor pred=model->forward(xt_dev);

      // We have to compute the derivative separately, one call to
      // autograd::grad() for each output dimension
      for(size_t j=0;j<out_dim;j++) {
        
        // Fourth argument is "retain_graph", and fifth is
        // "create_graph". 
        //torch::Tensor grad_out=torch::zeros({((long int)1),
        //((long int)out_dim)},torch::kFloat64).to(device);
        //grad_out[0][j]=1.0;
        
        torch::Tensor grads=torch::autograd::grad
          ({pred[0][j]},{xt_dev},{},true,true)[0].to(torch::kCPU);
        torch::TensorAccessor<double,2> acc=
          grads.accessor<double,2>();
        
        y[j]=acc[0][ix];
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
      std::cout << "n_pts: " << n_pts << std::endl;
      
      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt_cpu=torch::from_blob
        (vpp,{((long int)n_pts),((long int)in_dim)},
         torch::kFloat64).clone();
      
      torch::Tensor xt=xt_cpu.to(device);
      xt=xt.detach().requires_grad_(true);
      torch::Tensor pred=model->forward(xt);

      for(size_t j=0;j<out_dim;j++) {
      // We have to compute the derivative separately, one call to
      // autograd::grad() for each output dimension

        // Fourth argument is "retain_graph", and fifth is
        // "create_graph". 
        //torch::Tensor grad_out=torch::zeros({((long int)n_pts),
        //((long int)out_dim)},torch::kFloat64).to(device);
        
        //torch::Tensor xt_k=xt[k].unsqueeze(0).clone().detach().requires_grad_(true);
        //torch::Tensor pred_k=model->forward(xt_k);
        
        //grad_out[k][j]=1.0;

        //std::cout << "Here1." << std::endl;
        //std::cout << "Here2." << std::endl;
        
        //torch::Tensor grads=torch::autograd::grad
        //({pred_k[0][j]},{xt_k},{},true,false)[0];
          
        for(size_t k=0;k<n_pts;k++) {
          
        torch::Tensor grads=torch::autograd::grad
          ({pred[k][j]},
           {xt},{},true,true)[0];
        
        torch::Tensor grads_cpu=grads.to(torch::kCPU).contiguous();
        torch::TensorAccessor<double,2> grads_acc=
          grads_cpu.accessor<double,2>();
          
          std::vector<size_t> index={k,j};
            //index[0]=k;
            //index[1]=j;
          y.get(index)=grads_acc[0][ix];
        }
        
      }

      return 0;
    }
    
    /*
    virtual int deriv2(const vec_t &x, vec_t &y, size_t i, size_t j) const {
      
      model->eval();
      
      void *vpp=(void *)(x.data());
      torch::Tensor xt=torch::from_blob(vpp,{((long int)1),
					     ((long int)n_in)},
	torch::kFloat64).to(device);
      xt.set_requires_grad(true);
      
      torch::Tensor pred=model->forward(test_pts).squeeze();

      torch::Tensor grads=torch::autograd::grad
        ({out},{test_pts},{torch::ones_like(out)},true,true)[0];

      Tensor d1_pred=grads.index({torch::indexing::Slice(),i});

      torch::Tensor gradx2=torch::autograd::grad
        ({d1_pred},{test_pts},{torch::ones_like(out)},true,true)[0];

      Tensor d2_pred=gradx2.index({torch::indexing::Slice(),j});
      
      for(size_t j=0;j<n_out;j++) {
        y[j]=dx_pred(0,j);
      }

      return 0;
    }
    */
    
  };
  
}

#endif

#endif



