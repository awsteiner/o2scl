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
/** \file transform.h
    \brief File for transformations
*/
#ifndef O2SCL_TRANSFORM_H
#define O2SCL_TRANSFORM_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <gsl/gsl_cdf.h>

#include <o2scl/tensor.h>
#include <o2scl/rng.h>

namespace o2scl {

  /** \brief Split a matrix-like data set into separate matrices for
      testing and training
   */
  template<class mat_t=o2scl::tensor2<>>
  int train_test_split(const mat_t &data, double frac_test,
                       mat_t &train, mat_t &test) {
    o2scl::rng<> r;
    r.clock_seed();
    return train_test_split_rng(data,frac_test,r,train,test);
  }

  /** \brief Split a matrix-like data set into separate matrices for
      testing and training (with RNG)
   */
  template<class mat_t=o2scl::tensor2<>,
           class mat_col_t=o2scl::const_matrix_column_gen<o2scl::tensor2<> > >
  int train_test_split_rng(const mat_t &data, double frac_test,
                           o2scl::rng<> &r, mat_t &train, mat_t &test) {
                           
    
    if (frac_test<0.0 || frac_test>1.0 || !std::isfinite(frac_test)) {
      O2SCL_ERR2("Invalid fraction in ",
                 "train_test_split_rng().",o2scl::exc_einval);
    }
    
    size_t n_pts=data.size1();
    size_t n_out=data.size2();
    size_t n_test=n_pts*frac_test;
    size_t n_train=n_pts-n_test;
    if (n_test==0 || n_train==0) {
      O2SCL_ERR2("Fraction is too large or small to result in ",
                "non-zero testing set.",o2scl::exc_einval);
    }

    train.resize(n_pts-n_test,n_out);
    test.resize(n_test,n_out);
    
    std::vector<size_t> ix(n_pts);
    for(size_t i=0;i<n_pts;i++) {
      ix[i]=i;
    }
    o2scl::vector_shuffle<std::vector<size_t>,size_t>(r,n_pts,ix);
    
    for(size_t i=0;i<n_out;i++) {
      mat_col_t col(data,i);
      for(size_t j=0;j<n_test;j++) {
        test(j,i)=data(ix[j],i);
      }
      size_t j2=0;
      for(size_t j=n_test;j<n_pts;j++,j2++) {
        train(j2,i)=data(ix[j],i);
      }
    }
    
    return 0;
  }
  
  /** \brief Base class for data transformations
   */
  template<class mat_t=o2scl::tensor2<>>
  class transform_base {
    
  public:

    /** \brief Compute the transformation from the matrix \c t
     */
    virtual void fit(const mat_t &t)=0;

    /** \brief Perform an already computed transformation to the
        matrix \c t
    */
    virtual void transform(mat_t &t) const=0;

    /** \brief Compute the transformation from the matrix \c t
        and then perform the transformation
    */
    virtual mat_t fit_transform(const mat_t &t) {
      fit(t);
      mat_t out=t;
      transform(out);
      return out;
    }
      
    /** \brief Perform the inverse of the transformation on the
        matrix \c t
    */
    virtual void inverse_transform(mat_t &t) const=0;
    
  };

  /** \brief Quantile transformation
      
      Supports two output distributions:
      - "uniform" : maps each feature to [0, 1]
      - "normal" : maps each feature to N(0,1) via the probit function
 
      The transformer operates on the *last* axis of the tensor (i.e.
      axis rank-1), treating every earlier axis as the sample/batch
      dimensions. This mirrors sklearn's convention where axis-0 is
      samples and axis-1 is features for a 2-D array.
  */
  template<class mat_t=o2scl::tensor2<>>
  class transform_quantile : public transform_base<mat_t> {

  public:

    /** \brief Create a quantile transformer

        \param output_distribution "uniform" or "normal"
        \param n_quantiles  Number of quantile reference points to store
        per feature (sklearn default: 1000).
        \param subsample  Max samples used to estimate quantiles
        (0=use all).
        \param clip   If true, clamp output to [0,1] (uniform) or
        the probit of [eps, 1-eps] (normal) to avoid
        infinities at the boundaries.
    */
    explicit transform_quantile(std::string output_distribution="normal",
                                size_t n_quantiles=1000,
                                size_t subsample=0,
                                bool clip=true) :
      output_dist(std::move(output_distribution)),
      n_quantiles_(n_quantiles),
      subsample_(subsample),
      clip_(clip),
      fitted_(false) {

      if (output_dist !="uniform" && output_dist !="normal") {
        O2SCL_ERR(
                  "output_distribution must be 'uniform' or 'normal'",
                  o2scl::exc_einval);
      }
      if (n_quantiles_<2) {
        O2SCL_ERR("n_quantiles must be >=2",o2scl::exc_einval);
      }
    }

    /** \brief Compute the transformation from the matrix \c t
 
        The last axis of \p t is the feature axis; all other axes are treated
        as batch/sample dimensions and are flattened together.
        
        \param t Input tensor (not modified).
    */
    void fit(const mat_t &t) {

      const size_t rank=2;
      const size_t n_features=t.get_size(rank-1);
      const size_t n_samples=t.total_size()/n_features;

      // Build reference quantile positions [0/(q-1) ... 1/(q-1) ... 1]
      size_t q=n_quantiles_;
      refs_.resize(q);
      for (size_t i=0;i<q;++i) {
        refs_[i]=
          static_cast<double>(i)/static_cast<double>(q-1);
      }

      quantiles_.assign(n_features,std::vector<double>());

      for (size_t f=0;f<n_features;f++) {

        std::vector<double> col=
          extract_feature(t,f,n_samples,n_features);

        if (subsample_ > 0 && col.size() > subsample_) {
          for (size_t i=0;i<subsample_;i++) {
            size_t j=i+static_cast<size_t>(std::rand())
              %(col.size()-i);
            std::swap(col[i],col[j]);
          }
          col.resize(subsample_);
        }

        std::sort(col.begin(),col.end());
        quantiles_[f]=build_quantile_table(col);
      }

      fitted_=true;
      n_features_=n_features;

      return;
    }

    /** \brief Perform an already computed transformation to the
        matrix \c t
        
        Each element v of feature f is replaced by:
        -(uniform) the interpolated CDF value in [0, 1]
        -(normal) probit(CDF value) ~ N(0,1)
    */
    virtual void transform(mat_t &t) const {
      apply_transform(t,false);
    }

    /** \brief Perform the inverse of the transformation on the
        matrix \c t
        
        For "normal" output: first apply Phi (normal CDF) to get [0,1], then
        invert the quantile mapping.
    */
    virtual void inverse_transform(mat_t &t) const {
      apply_transform(t,true);
    }

    size_t n_quantiles() const { return n_quantiles_; }

    /// Return the quantile table for feature \p f (length==n_quantiles_).
    const std::vector<double> &quantile_table(size_t f) const {
      if (f >=quantiles_.size()) {
        O2SCL_ERR("feature index out of range",o2scl::exc_einval);
      }
      return quantiles_[f];
    }

  protected:

    /// Desc
    std::string output_dist;
    
    /// Number of quantiles requested
    size_t n_quantiles_;
    
    /// Desc
    size_t subsample_;
    
    /// Desc
    bool clip_;
    
    /// If true, a transformation has been computed (default false)
    bool fitted_;
    
    /// Desc
    size_t n_features_{0};

    /// quantiles_[f][q]=value at quantile refs_[q] for feature f
    std::vector<std::vector<double>> quantiles_;
    /// Uniform reference positions in [0,1] with length n_quantiles_
    std::vector<double> refs_;

    /** \brief Desc
        
        Extract all values of the last-axis feature index \p f into a flat
        vector. 
    */
    std::vector<double> extract_feature(const mat_t &t,
                                        size_t f,
                                        size_t n_samples,
                                        size_t n_features) const {
      std::vector<double> col;
      col.reserve(n_samples);
      for (size_t s=0;s<n_samples;++s) {
        col.push_back(static_cast<double>(t(s,f)));
      }
      // suppress unused parameter warning
      (void)n_features;
      return col;
    }

    /** \brief Desc
     */
    std::vector<double> build_quantile_table
    (const std::vector<double> &sorted_col) const {

      const size_t n=sorted_col.size();
      std::vector<double> qtab(n_quantiles_);
      for (size_t qi=0;qi<n_quantiles_;++qi) {
        // refs_ is valid here because fit() no longer clears it
        double p=refs_[qi];
        double pos=p*static_cast<double>(n-1);
        size_t lo=static_cast<size_t>(std::floor(pos));
        size_t hi=static_cast<size_t>(std::ceil(pos));
        if (hi >=n) hi=n-1;
        double frac=pos-static_cast<double>(lo);
        qtab[qi]=
          sorted_col[lo]*(1.0-frac)+sorted_col[hi]*frac;
      }
      return qtab;
    }

    /** \brief Desc Forward quantile lookup (value → CDF)
        Given a value \p v for feature \p f, return the interpolated uniform
        CDF in [0, 1] using the stored quantile table.
    */
    double value_to_uniform(double v, size_t f) const {
      const std::vector<double> &qtab=quantiles_[f];
      const size_t q=qtab.size();

      std::vector<double>::const_iterator it=
        std::lower_bound(qtab.begin(),qtab.end(),v);

      if (it==qtab.begin()) return refs_[0];
      if (it==qtab.end()) return refs_[q-1];

      size_t hi=static_cast<size_t>(it-qtab.begin());
      size_t lo=hi-1;

      double v_lo=qtab[lo], v_hi=qtab[hi];
      double dv=v_hi-v_lo;
      double frac=(dv==0.0) ? 0.5 : (v-v_lo)/dv;
      frac=std::clamp(frac, 0.0, 1.0);

      return refs_[lo]+frac*(refs_[hi]-refs_[lo]);
    }

    /** \brief Inverse quantile lookup (CDF → value)
        Given a uniform quantile \p u in [0,1] for feature \p f, return the
        interpolated original-scale value.
    */
    double uniform_to_value(double u, size_t f) const {
      u=std::clamp(u, 0.0, 1.0);
      const std::vector<double> &qtab=quantiles_[f];
      const size_t q=qtab.size();

      std::vector<double>::const_iterator it=
        std::lower_bound(refs_.begin(),refs_.end(),u);

      if (it==refs_.begin()) return qtab[0];
      if (it==refs_.end()) return qtab[q-1];

      size_t hi=static_cast<size_t>(it-refs_.begin());
      size_t lo=hi-1;

      double r_lo=refs_[lo], r_hi=refs_[hi];
      double dr=r_hi-r_lo;
      double frac=(dr==0.0) ? 0.5 : (u-r_lo)/dr;
      frac=std::clamp(frac, 0.0, 1.0);

      return qtab[lo]+frac*(qtab[hi]-qtab[lo]);
    }

    /** \brief Core transform loop
     */
    void apply_transform(mat_t &t, bool inverse) const {

      const size_t rank=2;
      const size_t n_features=t.get_size(rank-1);

      if (n_features !=n_features_) {
        O2SCL_ERR(
                  (((std::string)"transform_quantile: tensor's ")+
                   "last dimension ("+
                   std::to_string(n_features)+
                   ") does not match fitted n_features ("+
                   std::to_string(n_features_)+")").c_str(),
                  o2scl::exc_einval);
      }

      const size_t total=t.total_size();
      const size_t n_samp=total/n_features;
      const double eps=1e-7;

      for (size_t s=0;s<n_samp;++s) {
        for (size_t f=0;f<n_features;++f) {

          double v=static_cast<double>(t(s,f));

          if (!inverse) {
            double u=value_to_uniform(v,f);

            if (output_dist=="uniform") {
              if (clip_) u=std::clamp(u, 0.0, 1.0);
              t.get(s,f)=u;
            } else {
              if (clip_) u=std::clamp(u, eps, 1.0-eps);
              t.get(s,f)=gsl_cdf_ugaussian_Pinv(u);
            }
          } else {
            double u;
            if (output_dist=="uniform") {
              u=std::clamp(v, 0.0, 1.0);
            } else {
              u=gsl_cdf_ugaussian_P(v);
              u=std::clamp(u, 0.0, 1.0);
            }
            t.get(s,f)=uniform_to_value(u,f);
          }
        }
      }
    }

  };  

  /** /brief Rescale data to a new minimum and maximum value
      
      Each output quantity (each column) is independently rescaled.
      The forward transform is
      \f[
      x^{\prime}=(x-x_{\mathrm{min}ad_rmf.h
      eos_h})/(x_{mathrm{max}}-x_{\mathrm{min}})
      (r_{\mathrm{max}}-r_{\mathrm{min}})+r_{\mathrm{min}} \, .
      \f]
      This transformation is invertible.

      Constant quantities with \f$x_{\mathrm{max}}=
      x_{\mathrm{min}}\f$ are mapped to \f$ r_{\mathrm{min}} \f$,
      matching <tt>sklearn</tt>'s behaviour.
  */
  template<class mat_t=o2scl::tensor2<>>
  class transform_minmax : public transform_base<mat_t> {
    
  public:
    
    /** \brief Create a transformation object from
        \c range_min to \c range_max
    */
    explicit transform_minmax(double range_min=0.0, double range_max=1.0) {
      if (range_min_>=range_max_) {
        O2SCL_ERR2("transform_minmax: range_min must ",
                   "be strictly less than range_max",o2scl::exc_einval);
      }
                  
      range_min_=range_min;
      range_max_=range_max;
      fitted_=false;
    }

    /** \brief Compute the transformation from the matrix \c t
     */
    void fit(const mat_t &t) {
      
      const size_t n_features=t.size2();
      const size_t n_samples=t.size1();

      data_min_.assign(n_features,std::numeric_limits<double>::infinity());
      data_max_.assign(n_features,-std::numeric_limits<double>::infinity());

      for (size_t s=0;s<n_samples;++s) {
        for (size_t f=0;f<n_features;++f) {
          double v=t(s,f);
          if (!std::isfinite(v)) {
            O2SCL_ERR2("Non-finite value in ",
                       "transform_minmax::fit().",o2scl::exc_einval);
          }
          if (v<data_min_[f]) data_min_[f]=v;
          if (v > data_max_[f]) data_max_[f]=v;
        }
      }

      // Pre-compute per-feature scale and shift so transform() is
      // O(1) per element. x_scaled=x*scale_[f]+shift_[f]
      scale_.resize(n_features);
      shift_.resize(n_features);
      
      const double range_span=range_max_-range_min_;
      for (size_t f=0;f<n_features;++f) {
        double span=data_max_[f]-data_min_[f];
        if (span==0.0) {
          // Constant feature: map to range_min (mirrors sklearn)
          scale_[f]=0.0;
          shift_[f]=range_min_;
        } else {
          scale_[f]=range_span/span;
          shift_[f]=range_min_-data_min_[f]*scale_[f];
        }
      }

      fitted_=true;
      n_features_=n_features;
    }

    /** \brief Perform an already computed transformation to the
        matrix \c t
    */
    virtual void transform(mat_t &t) const {

      if (!fitted_) {
        O2SCL_ERR2("Called transform() before fit() in ",
                   "transform_minmax::transform()",
                   o2scl::exc_einval);
      }

      const size_t n_features=t.size2();
      const size_t n_samples=t.size1();
      
      for (size_t s=0;s<n_samples;++s) {
        for (size_t f=0;f<n_features;++f) {
          t.get(s,f)=t.get(s,f)*scale_[f]+shift_[f];
        }
      }
      return;
    }

    /** \brief Perform the inverse of the transformation on the
        matrix \c t
        
        Constant features are restored
        to x_{\mathrm{min}}.
    */
    virtual void inverse_transform(mat_t &t) const {
      
      if (!fitted_) {
        O2SCL_ERR2("Called inverse_transform() before fit() in ",
                   "transform_minmax::inverse_transform()",
                   o2scl::exc_einval);
      }

      const double range_span=range_max_-range_min_;
      const size_t n_features=t.size2();
      const size_t n_samples=t.size1();
      
      for (size_t s=0;s<n_samples;++s) {
        for (size_t f=0;f<n_features;++f) {
          
          double v=t(s,f);
          double x;
          if (scale_[f]==0.0) {
            x=data_min_[f]; // constant feature
          } else {
            x=(v-range_min_)/range_span*
              (data_max_[f]-data_min_[f])+data_min_[f];
          }
          t.get(s,f)=x;
        }
      }

      return;
    }

    double feature_min(size_t f) const { return data_min_.at(f); }
    double feature_max(size_t f) const { return data_max_.at(f); }
    double scale(size_t f) const { return scale_.at(f); }

  protected:
 
    double range_min_, range_max_;
    
    /// If true, a transformation has been computed (default false)
    bool fitted_;
    size_t n_features_{0};

    std::vector<double> data_min_, data_max_;
    std::vector<double> scale_, shift_; 

    void check_features(mat_t &t) const {
      const size_t rank=2;
      if (t.get_size(rank-1) !=n_features_)
        O2SCL_ERR((((std::string)"transform_minmax: tensor's ")+
                   "last dimension ("+
                   std::to_string(t.get_size(rank-1))+
                   ") does not match fitted n_features ("+
                   std::to_string(n_features_)+")").c_str(),
                  o2scl::exc_einval);
    }

  };
 
  /** \brief Transform to data to zero mean and unit variance
 
      The forward transform is:
 
      x_scaled=(x-mean)/std  (with_std=true, default)
      x_scaled=(x-mean)  (with_std=false)
      x_scaled=x/std   (with_mean=false, with_std=true)
 
      Features with zero variance are left as is. (std=1
      internally).
 
      The transformer operates on the *last* axis of the tensor,
      treating every earlier axis as the sample/batch dimensions
      (row-major flat storage).
  */
  template<class mat_t=o2scl::tensor2<>>
  class transform_standard : public transform_base<mat_t> {
 
  protected:

    /// If true, include a shift of the mean
    bool with_mean_;

    /// If true, Include a rescaling to unit standard deviation
    bool with_std_;

    /// If true, a transformation has been computed (default false)
    bool fitted_;

    /// The number of data sets (the number of matrix columns)
    size_t n_features_;

    /// The vector of means
    std::vector<double> mean_;

    /// The vector of standard deviations
    std::vector<double> std_;

    /** \brief Compute the transformation
     */
    void apply(mat_t &t, bool inverse) const {

      const size_t n_features=t.size2();
      const size_t n_samples=t.size1();

      for (size_t s=0;s<n_samples;s++) {
        for (size_t f=0;f<n_features;f++) {
          double v=t(s,f);

          if (!inverse) {
            if (with_mean_) v-=mean_[f];
            if (with_std_) v/=std_[f];
          } else {
            if (with_std_) v*=std_[f];
            if (with_mean_) v+=mean_[f];
          }

          t(s,f)=v;
        }
      }
      return;
    }

  public:

    /** \brief Create a standard transformation
        
        \param with_mean If true (default), subtract the per-feature
        mean.
        \param with_std If true (default), divide by the
        per-feature population standard deviation. Zero-variance
        features keep std=1.
    */
    explicit transform_standard(bool with_mean=true,
                                bool with_std=true) {
      with_mean_=with_mean;
      with_std_=with_std;
      fitted_=false;
      n_features_=0;
    }
    
    /** \brief Compute the transformation from the matrix \c t
        
        Compute means and standard deviations from \p t (not
        modified). 
    */
    void fit(const mat_t &t) {

      const size_t n_features=t.size2();
      const size_t n_samples=t.size1();

      mean_.assign(n_features,0.0);
      // default=1 (no scaling for zero-var)
      std_.assign(n_features,1.0); 

      if (with_mean_ || with_std_) {
 
        // Welford's single-pass algorithm for mean and M2 (sum of
        // squared diffs)
        std::vector<double> M2(n_features, 0.0);
        for (size_t s=0;s<n_samples;++s) {
          for (size_t f=0;f<n_features;++f) {
            double v=t(s,f);
            double delta=v-mean_[f];
            mean_[f]+=delta/(s+1);
            double delta2=v-mean_[f];
            M2[f]+=delta*delta2;
          }
        }

        if (with_std_) {
          for (size_t f=0;f<n_features;++f) {
            // population
            double var=M2[f]/(n_samples); 
            std_[f]=(var==0.0) ? 1.0 : std::sqrt(var);
          }
        }
      }

      fitted_=true;
      n_features_=n_features;
    }

    /** \brief Perform an already computed transformation to the
        matrix \c t
    */
    virtual void transform(mat_t &t) const {
      if (!fitted_) {
        O2SCL_ERR2("Called fit() before transform() in ",
                   "transform_standard::transform()",o2scl::exc_einval);
      }
      apply(t,false);
    }

    /** \brief Perform the inverse of the transformation on the
        matrix \c t
    */
    virtual void inverse_transform(mat_t &t) const {
      if (!fitted_) {
        O2SCL_ERR2("Called fit() before transform() in ",
                   "transform_standard::inverse_transform()",
                   o2scl::exc_einval);
      }
      apply(t,true);
    }

    /** \brief Return the vector of means
     */
    std::vector<double> get_means() {
      return mean_;
    }
    
    /** \brief Return the vector of standard deviations
     */
    std::vector<double> get_stds() {
      return std_;
    }
    
  };
 
}

#endif
