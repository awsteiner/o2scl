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
#ifndef O2SCL_INTERPM_RBF_LIBTORCH_H
#define O2SCL_INTERPM_RBF_LIBTORCH_H

#if defined(O2SCL_SET_LIBTORCH) || defined(DOXYGEN)

#include <torch/torch.h>
#include <o2scl/interpm_base.h>
#include <o2scl/tensor.h>
#include <o2scl/set_libtorch.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

/** \file interpm_rbf_libtorch.h
 */
namespace o2scl {

  /** \brief An RBF network with learned centers and widths
   */
  struct torch_rbf : torch::nn::Module {

    /// RBF centers, shape [n_centers, input_dim]
    torch::Tensor centers;

    /// RBF log-widths (log for positivity), shape [n_centers]
    torch::Tensor log_widths;

    /// Output linear layer
    torch::nn::Linear output_layer{nullptr};

    /// Number of input dimensions
    size_t input_dim;

    /// Number of output dimensions
    size_t output_dim;

    /// Number of RBF centers
    size_t n_centers;

    /// RBF kernel type (default "gaussian")
    std::string rbf_type;

    /** \brief Create an empty network
     */
    torch_rbf() {
      input_dim=0;
      output_dim=0;
      n_centers=0;
      rbf_type="gaussian";
    }

    /** \brief Create an RBF network
     */
    torch_rbf(size_t in_dim, size_t n_cen, size_t out_dim) {
      rbf_type="gaussian";
      init(in_dim,n_cen,out_dim);
    }

    /** \brief Initialize the RBF network
     */
    void init(size_t in_dim, size_t n_cen, size_t out_dim) {

      input_dim=in_dim;
      output_dim=out_dim;
      n_centers=n_cen;

      // Initialize centers randomly in [-1,1]
      centers=register_parameter
        ("centers",torch::randn({(int)n_centers,(int)input_dim},
                                torch::kFloat64));
      
      // Initialize log-widths to zero (width=1)
      log_widths=register_parameter
        ("log_widths",torch::zeros({(int)n_centers},
                                   torch::kFloat64));
      
      // Output layer maps from n_centers to output_dim
      output_layer=register_module
        ("output_layer",
         torch::nn::Linear(n_centers,(int)output_dim));

      return;
    }

    /** \brief Compute RBF activations
     */
    torch::Tensor rbf_activations(torch::Tensor x) {

      // x shape: [batch, input_dim]
      // centers shape: [n_centers, input_dim]
      // Expand for broadcasting:
      // x_exp shape: [batch, 1, input_dim]
      // c_exp shape: [1, n_centers, input_dim]
      torch::Tensor x_exp=x.unsqueeze(1);
      torch::Tensor c_exp=centers.unsqueeze(0);

      // Squared distances: [batch, n_centers]
      torch::Tensor diff=x_exp-c_exp;
      torch::Tensor dist2=(diff*diff).sum(2);

      // Widths: [1, n_centers]
      torch::Tensor widths=torch::exp(log_widths).unsqueeze(0);
      torch::Tensor widths2=widths*widths;

      if (rbf_type=="gaussian") {
        return torch::exp(-dist2/widths2);
      } else if (rbf_type=="multiquadric") {
        return torch::sqrt(dist2/widths2+1.0);
      } else if (rbf_type=="inverse_multiquadric") {
        return 1.0/torch::sqrt(dist2/widths2+1.0);
      } else if (rbf_type=="linear") {
        return torch::sqrt(dist2/widths2+1.0e-12);
      } else if (rbf_type=="thin_plate") {
        torch::Tensor r=torch::sqrt(dist2/widths2+1.0e-12);
        return r*r*torch::log(r+1.0e-12);
      } else {
        O2SCL_ERR("Invalid RBF type.",o2scl::exc_einval);
      }
      // Unreachable, but silences compiler warning
      return torch::zeros({1});
    }

    /** \brief Forward propagation
     */
    torch::Tensor forward(torch::Tensor x) {
      torch::Tensor phi=rbf_activations(x);
      return output_layer->forward(phi);
    }

  };

  /** \brief Multidimensional interpolation with an RBF network
      and libtorch (experimental)

      \verbatim embed:rst
      See also the :ref:`Higher-dimensional Interpolation`
      section of the User's guide.
      \endverbatim
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_x_t=o2scl::const_matrix_view_table<>,
           class mat_y_t=o2scl::matrix_view_table<>>
  class interpm_rbf_libtorch :
    public interpm_base<vec_t,mat_x_t,mat_y_t> {

  protected:

    /// RBF network
    std::shared_ptr<torch_rbf> model;

    /// Input dimensions
    size_t in_dim;

    /// Output dimensions
    size_t out_dim;

  public:

    interpm_rbf_libtorch() : device(torch::kCPU) {
      epochs=1000;
      patience=500;
      n_centers=64;
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

    /// Number of RBF centers (default 64)
    size_t n_centers;

    /// Verbosity (default 1)
    int verbose;

    /// Adam learning rate (default 1.0e-3)
    double adam_lr;

    /// Device
    torch::Device device;

    /// The last device used
    std::string dev_str;

    /// The random seed, or 0 for automatic (default 0)
    int seed;

    /// Number of steps between sparse verbose output (default 50)
    int epoch_step;

    /** \brief Set the data to be interpolated
     */
    virtual int set_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &user_x, mat_y_t &user_y) {

      o2scl::tensor<> tin, tout;
      std::vector<size_t> in_size={n_pts,n_in};
      std::vector<size_t> out_size={n_pts,n_out};
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
        std::cout << "interpm_rbf_libtorch::set_data_tensor(): "
                  << "Device: " << dev_str << std::endl;
      }

      // ──────────────────────────────────────────────────────────
      // Convert input o2scl tensor objects to torch tensor objects

      void *vpp=(void *)(params.get_data().data());
      torch::Tensor tp=torch::from_blob
        (vpp,{((long int)n_pts),((long int)n_in)},
         torch::kFloat64).to(device);
      void *vpo=(void *)(outputs.get_data().data());
      torch::Tensor to=torch::from_blob
        (vpo,{((long int)n_pts),((long int)n_out)},
         torch::kFloat64).to(device);

      // ──────────────────────────────────────────────────────────

      if (seed!=0) {
        torch::manual_seed(seed);
      }

      model=std::make_shared<torch_rbf>(n_in,n_centers,n_out);

      // Initialize centers from a random subset of training data
      // so they start in a sensible region of input space
      {
        torch::Tensor perm=torch::randperm
          ((long int)n_pts,torch::kLong);
        torch::Tensor idx=perm.index
          ({torch::indexing::Slice(0,(long int)n_centers)});
        torch::Tensor init_centers=tp.index_select(0,idx.to(device))
          .detach().clone().to(torch::kCPU);
        torch::NoGradGuard no_grad;
        model->centers.copy_(init_centers);
      }

      model->to(device,torch::kFloat64);

      torch::optim::Adam optimizer(model->parameters(),
                                   torch::optim::AdamOptions(adam_lr));

      // ──────────────────────────────────────────────────────────
      // Train/validation split

      int total_samples=n_pts;

      torch::Tensor train_inputs, train_targets;
      torch::Tensor val_inputs, val_targets;

      if (test_size>0.0) {

        int train_size=static_cast<int>
          ((1.0-test_size)*total_samples);

        torch::Tensor perm=torch::randperm
          (total_samples,torch::kLong).to(device);

        torch::Tensor train_idx=perm.index
          ({torch::indexing::Slice(0,train_size)}).to(device);
        torch::Tensor val_idx=perm.index
          ({torch::indexing::Slice
              (train_size,torch::indexing::None)}).to(device);

        train_inputs=tp.index_select(0,train_idx).to(device);
        train_targets=to.index_select(0,train_idx).to(device);

        val_inputs=tp.index_select(0,val_idx).to(device);
        val_targets=to.index_select(0,val_idx).to(device);

      } else {

        train_inputs=tp.to(device);
        train_targets=to.to(device);

      }

      // ──────────────────────────────────────────────────────────
      // Early stopping parameters

      int patience_counter=0;
      double best_loss=std::numeric_limits<double>::infinity();
      std::vector<torch::Tensor> best_params;

      // ──────────────────────────────────────────────────────────
      // Training

      for (int epoch=0;epoch<epochs;epoch++) {

        model->train();
        optimizer.zero_grad();

        torch::Tensor pred=model->forward(train_inputs);
        torch::Tensor loss=torch::mse_loss(pred,train_targets);

        loss.backward();
        optimizer.step();

        // ────────────────────────────────────────────────────────
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

        // ────────────────────────────────────────────────────────
        // Early stopping logic

        if ((test_size>0.0 && val_loss<best_loss) ||
            (test_size==0.0 && train_loss<best_loss)) {

          if (test_size>0.0) {
            best_loss=val_loss;
          } else {
            best_loss=train_loss;
          }

          patience_counter=0;

          best_params.clear();
          std::vector<torch::Tensor> params=model->parameters();
          for (size_t i=0;i<params.size();++i) {
            best_params.push_back(params[i].detach().clone());
          }

        } else {

          patience_counter++;

        }

        if (verbose>0) {
          if (epoch%epoch_step==0 || verbose>1) {
            std::cout
              << "interpm_rbf_libtorch::set_data_tensor(): ";
            std::cout << "Epoch ";
            std::cout.width(4);
            std::cout << epoch
                      << " Train loss: " << train_loss;
            if (test_size>0.0) {
              std::cout << " Val loss: " << val_loss;
            }
            std::cout << std::endl;
          }
        }

        if (patience_counter>=patience) {
          break;
        }
      }

      // ──────────────────────────────────────────────────────────
      // Restore best model

      for (size_t i=0;i<best_params.size();++i) {
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
                   "interpm_rbf_libtorch::eval().",
                   o2scl::exc_einval);
      }

      model->eval();
      torch::NoGradGuard no_grad;

      torch::TensorOptions opts=torch::TensorOptions().dtype
        (torch::kFloat64).device(torch::kCPU);
      torch::Tensor xt=torch::empty
        ({(long int)1,(long int)in_dim},opts);

      double *xt_ptr=xt.data_ptr<double>();
      for(size_t i=0;i<x.size();i++) {
        xt_ptr[i]=x[i];
      }

      torch::Tensor xt_dev=xt.to(device);
      torch::Tensor pred=model->forward(xt_dev);
      torch::Tensor pred_cpu=pred.to(torch::kCPU);
      double *pred_ptr=pred_cpu.data_ptr<double>();

      for(size_t j=0;j<out_dim;j++) {
        y[j]=pred_ptr[j];
      }

      return 0;
    }

    /** \brief Evaluate the interpolation at a list of points,
        returning results in \c y (tensor form)
    */
    virtual int eval_list_tensor(const o2scl::tensor<> &x,
                                 o2scl::tensor<> &y) const {

      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_rbf_libtorch::eval_list_tensor().",
                   o2scl::exc_einval);
      }

      model->eval();
      torch::NoGradGuard no_grad;

      size_t n_pts=x.get_size(0);

      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt=torch::from_blob
        (vpp,{((long int)n_pts),((long int)in_dim)},
         torch::kFloat64).to(device);

      torch::Tensor pred=model->forward(xt);
      torch::Tensor pred_cpu=pred.to(torch::kCPU);
      torch::TensorAccessor<double,2> pred_acc=
        pred_cpu.accessor<double,2>();

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

    /** \brief Evaluate the derivative with respect to variable
        \c ix at point \c x, returning \c y
    */
    virtual int deriv(const vec_t &x, vec_t &y, size_t ix) const {

      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_rbf_libtorch::deriv().",
                   o2scl::exc_einval);
      }

      model->eval();

      torch::TensorOptions opts=torch::TensorOptions().dtype
        (torch::kFloat64).device(torch::kCPU);
      torch::Tensor xt=torch::empty
        ({(long int)1,(long int)in_dim},opts);

      double *xt_ptr=xt.data_ptr<double>();
      for(size_t i=0;i<x.size();i++) {
        xt_ptr[i]=x[i];
      }

      torch::Tensor xt_dev=xt.clone().to(device);
      xt_dev.requires_grad_(true);
      torch::Tensor pred=model->forward(xt_dev);

      for(size_t j=0;j<out_dim;j++) {

        torch::Tensor grad_out=torch::zeros({((long int)1),
            ((long int)out_dim)},torch::kFloat64).to(device);
        grad_out[0][j]=1.0;

        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt_dev},{grad_out},
           true,false)[0].to(torch::kCPU);
        double *grads_ptr=grads.data_ptr<double>();

        y[j]=grads_ptr[ix];
      }

      return 0;
    }

    /** \brief Evaluate the derivative with respect to variable
        \c ix at a list of points, storing results in \c y
    */
    virtual int deriv_list_tensor
    (const o2scl::tensor<> &x, o2scl::tensor<> &y,
     size_t ix) const {

      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_rbf_libtorch::deriv_list_tensor().",
                   o2scl::exc_einval);
      }

      model->eval();

      size_t n_pts=x.get_size(0);

      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt_cpu=torch::from_blob
        (vpp,{((long int)n_pts),((long int)in_dim)},
         torch::kFloat64);

      torch::Tensor xt=xt_cpu.clone().to(device);
      xt.requires_grad_(true);
      torch::Tensor pred=model->forward(xt);

      for(size_t j=0;j<out_dim;j++) {

        torch::Tensor grad_out=torch::zeros({((long int)n_pts),
            ((long int)out_dim)},torch::kFloat64).to(device);
        for(size_t k=0;k<n_pts;k++) {
          grad_out[k][j]=1.0;
        }

        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt},{grad_out},
           true,false)[0].to(torch::kCPU);

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

    /** \brief Compute the second derivative of outputs with respect
        to inputs ix and ix2 at point x, storing results in y
    */
    virtual int deriv2(const vec_t &x, vec_t &y,
                       size_t ix, size_t ix2) const {

      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_rbf_libtorch::deriv2().",
                   o2scl::exc_einval);
      }

      model->eval();

      torch::TensorOptions opts=torch::TensorOptions().dtype
        (torch::kFloat64).device(torch::kCPU);
      torch::Tensor xt=torch::empty
        ({(long int)1,(long int)in_dim},opts);

      double *xt_ptr=xt.data_ptr<double>();
      for(size_t i=0;i<x.size();i++) {
        xt_ptr[i]=x[i];
      }

      torch::Tensor xt_dev=xt.clone().to(device);
      xt_dev.requires_grad_(true);
      torch::Tensor pred=model->forward(xt_dev);

      for(size_t j=0;j<out_dim;j++) {

        torch::Tensor grad_out=torch::zeros({((long int)1),
            ((long int)out_dim)},torch::kFloat64).to(device);
        grad_out[0][j]=1.0;

        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt_dev},{grad_out},true,true)[0];

        torch::Tensor grads_ix=grads.index
          ({torch::indexing::Slice(),
            ((long int)ix)}).unsqueeze(1);

        torch::Tensor grad_out2=torch::ones_like(grads_ix);
        torch::Tensor grads2=torch::autograd::grad
          ({grads_ix},{xt_dev},{grad_out2},
           true,false)[0].to(torch::kCPU);

        double *grads2_ptr=grads2.data_ptr<double>();
        y[j]=grads2_ptr[ix2];
      }

      return 0;
    }

    /** \brief Compute the second derivative of outputs with respect
        to inputs ix and ix2 at a list of points, storing results
        in y
    */
    virtual int deriv2_list_tensor
    (const o2scl::tensor<> &x, o2scl::tensor<> &y,
     size_t ix, size_t ix2) const {

      if (!model) {
        O2SCL_ERR2("Model not specified in ",
                   "interpm_rbf_libtorch::deriv2_list_tensor().",
                   o2scl::exc_einval);
      }

      model->eval();

      size_t n_pts=x.get_size(0);

      void *vpp=(void *)(x.get_data().data());
      torch::Tensor xt_cpu=torch::from_blob
        (vpp,{((long int)n_pts),((long int)in_dim)},
         torch::kFloat64);

      torch::Tensor xt=xt_cpu.clone().to(device);
      xt.requires_grad_(true);
      torch::Tensor pred=model->forward(xt);

      for(size_t j=0;j<out_dim;j++) {

        torch::Tensor grad_out=torch::zeros({((long int)n_pts),
            ((long int)out_dim)},torch::kFloat64).to(device);
        for(size_t k=0;k<n_pts;k++) {
          grad_out[k][j]=1.0;
        }

        torch::Tensor grads=torch::autograd::grad
          ({pred},{xt},{grad_out},true,true)[0];

        torch::Tensor grads_ix=grads.index
          ({torch::indexing::Slice(),
            ((long int)ix)}).unsqueeze(1);

        torch::Tensor grad_out2=torch::ones_like(grads_ix);
        torch::Tensor grads2=torch::autograd::grad
          ({grads_ix},{xt},{grad_out2},
           true,false)[0].to(torch::kCPU);

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
                   "interpm_libtorch_rbf::save().",
                   o2scl::exc_einval);
      }
      
      // Start group
      hid_t top=hf.get_current_id();
      hid_t group=hf.open_group(name);
      hf.set_current_id(group);

      // Add typename
      hf.sets_fixed("o2scl_type","interpm_rbf_libtorch");
      
      // Save scalar metadata
      hf.set_szt("in_dim",in_dim);
      hf.set_szt("out_dim",out_dim);
      hf.set_szt("n_centers",n_centers);
      hf.sets("rbf_type",model->rbf_type);
      
      // Save each named parameter as a flat double vector
      // plus its shape as an int vector
      for (std::pair<std::string,torch::Tensor> entry :
             model->named_parameters()) {
        
        std::string pname=entry.first;
        torch::Tensor t=entry.second.detach()
          .to(torch::kCPU).to(torch::kFloat64).contiguous();
        
        // Save the shape
        std::vector<int> shape;
        for (int i=0;i<t.dim();i++) {
          shape.push_back((int)t.size(i));
        }
        hf.seti_vec(pname+"_shape",shape);
        
        // Save the flat data
        std::vector<double> data(t.data_ptr<double>(),
                                 t.data_ptr<double>()+t.numel());
        hf.setd_vec(pname+"_data",data);
      }
      
      // Also save buffers (e.g. batch norm running stats if present)
      for (std::pair<std::string,torch::Tensor> entry :
             model->named_buffers()) {
        
        std::string pname=entry.first;
        torch::Tensor t=entry.second.detach()
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
        hf.find_object_by_type("interpm_rbf_libtorch",name);
        if (name.length()==0) {
          O2SCL_ERR3("No object of type interpm_rbf_libtorch found in ",
                     "o2scl_hdf::hdf_input(hdf_file &,",
                     "interpm_rbf_libtorch &,string &).",exc_efailed);
        }
      }

      // Open main group
      hid_t top=hf.get_current_id();
      hid_t group=hf.open_group(name);
      hf.set_current_id(group);

      // Check typename
      std::string type2;
      hf.gets_fixed("o2scl_type",type2);
      if (type2!="interpm_rbf_libtorch") {
        O2SCL_ERR2("Typename in HDF group does not match ",
                   "class in table::hdf_input().",exc_einval);
      }
      
      // Load scalar metadata and reconstruct model
      std::string rbf_type_str;
      hf.get_szt("in_dim",in_dim);
      hf.get_szt("out_dim",out_dim);
      hf.get_szt("n_centers",n_centers);
      hf.gets("rbf_type",rbf_type_str);
      
      // Reconstruct model architecture before loading weights
      model=std::make_shared<torch_rbf>(in_dim,n_centers,out_dim);
      model->rbf_type=rbf_type_str;
      model->to(device,torch::kFloat64);

      // Load each named parameter
      for (std::pair<std::string,torch::Tensor> entry :
             model->named_parameters()) {
        
        std::string pname=entry.first;
        
        std::vector<int> shape;
        hf.geti_vec(pname+"_shape",shape);
        
        std::vector<double> data;
        hf.getd_vec(pname+"_data",data);
        
        // Reconstruct the shape for torch
        std::vector<long int> tshape(shape.begin(),shape.end());
        torch::Tensor t=torch::from_blob
          (data.data(),
           torch::IntArrayRef(tshape.data(),(size_t)tshape.size()),
           torch::kFloat64).clone().to(device);
        
        torch::NoGradGuard no_grad;
        entry.second.copy_(t);
      }
      
      // Load buffers
      for (std::pair<std::string,torch::Tensor> entry :
             model->named_buffers()) {
        
        std::string pname=entry.first;
        
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
        entry.second.copy_(t);
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
