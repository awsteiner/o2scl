/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

  -------------------------------------------------------------------
*/
/* monte/vegas.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_MCARLO_VEGAS_H
#define O2SCL_MCARLO_VEGAS_H

/** \file mcarlo_vegas.h
    \brief File defining \ref o2scl::mcarlo_vegas
*/

#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_monte_vegas.h>

#include <o2scl/mcarlo.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Multidimensional integration using Vegas Monte Carlo (GSL)

      The output options are a little different than the original GSL
      routine. The default setting of \ref mcarlo::verbose is 0,
      which turns off all output. A verbose value of 1 prints summary
      information about the weighted average and final result, while a
      value of 2 also displays the grid coordinates. A value of 3
      prints information from the rebinning procedure for each
      iteration.

      Some original documentation from GSL:

      \verbatim 
      The input coordinates are x[j], with upper and lower limits
      xu[j] and xl[j]. The integration length in the j-th direction is
      delx[j]. Each coordinate x[j] is rescaled to a variable y[j] in
      the range 0 to 1. The range is divided into bins with boundaries
      xi[i][j], where i=0 corresponds to y=0 and i=bins to y=1. The
      grid is refined (ie, bins are adjusted) using d[i][j] which is
      some variation on the squared sum. A third parameter used in
      defining the real coordinate using random numbers is called z.
      It ranges from 0 to bins. Its integer part gives the lower index
      of the bin into which a call is to be placed, and the remainder
      gives the location inside the bin.

      When stratified sampling is used the bins are grouped into
      boxes, and the algorithm allocates an equal number of function
      calls to each box.

      The variable alpha controls how "stiff" the rebinning algorithm
      is. alpha = 0 means never change the grid. Alpha is typically
      set between 1 and 2. 
      \endverbatim

      \verbatim embed:rst
      .. todo:: 
      
         CLass mcarlo_vegas: Mode = importance only doesn't give the
         same answer as GSL yet.

      \endverbatim

      \future Prettify the verbose output

      \future Allow the user to get information about the how
      the sampling was done, possibly by converting the bins
      and boxes into a structure or class.

      \future Allow the user to change the maximum number of bins.

      \verbatim embed:rst
      Based on [Lepage78]_
      \endverbatim

      The current version of the algorithm
      was described in the Cornell preprint CLNS-80/447 of March,
      1980. The GSL code follows most closely the C version by D. R.
      Yennie, coded in 1984.
  */
  template<class func_t=multi_funct, 
    class vec_t=boost::numeric::ublas::vector<double>,
    class rng_t=rng>
    class mcarlo_vegas : public mcarlo<func_t,vec_t,rng_t> {
    
    public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    typedef boost::numeric::ublas::vector<int> ubvector_int;

    /// \name Integration mode (default is mode_importance)
    //@{
    int mode;
    static const int mode_importance=1;
    static const int mode_importance_only=0;
    static const int mode_stratified=-1;
    //@}

    /// Result from last iteration
    double result;

    /// Uncertainty from last iteration
    double sigma;

    /** \brief The stiffness of the rebinning algorithm (default 1.5)

        This usual range is between 1 and 2.
    */
    double alpha;

    /// Set the number of iterations (default 5)
    unsigned int iterations;

    /** \brief The chi-squared per degree of freedom for the weighted
        estimate of the integral
        
        After an integration, this should be close to 1. If it is not,
        then this indicates that the values of the integral from
        different iterations are inconsistent, and the error may be
        underestimated. Further iterations of the algorithm may enable
        one to obtain more reliable results.
    */
    double chisq;
    
    /// The output stream to send output information (default \c std::cout)
    std::ostream *outs;

#ifndef DOXYGEN_INTERNAL

    protected:

    /// Maximum number of bins
    static const size_t bins_max=50;
    
    /// Number of dimensions
    size_t dim;
    /// Number of bins
    unsigned int bins;
    /// Number of boxes
    unsigned int boxes; 
    /// Boundaries for each bin
    ubvector xi;
    /// Storage for grid refinement
    ubvector xin;
    /// The iteration length in each direction
    ubvector delx;
    /// The weight for each bin
    ubvector weight;
    /// The volume of the current bin
    double vol;

    /// The bins for each direction
    ubvector_int bin;
    /// The boxes for each direction
    ubvector_int box;

    /// Distribution 
    ubvector d;

    /** \name Scratch variables preserved between calls to 
        vegas_minteg_err()
    */
    //@{
    double jac;
    double wtd_int_sum;
    double sum_wgts;
    double chi_sum;
    //@}
    
    /// The starting iteration number
    unsigned int it_start;
    /// The total number of iterations
    unsigned int it_num;
    /// Number of samples for computing chi squared
    unsigned int samples;
    /// Number of function calls per box
    unsigned int calls_per_box;
    
    /// Initialize box coordinates
    virtual void init_box_coord(ubvector_int &boxt) {
      size_t i;
      for (i=0;i<dim;i++) {
        boxt[i]=0;
      }
      return;
    }
  
    /** \brief Change box coordinates
        
        Steps through the box coordinates, e.g.
        \verbatim
        {0,0}, {0,1}, {0,2}, {0,3}, {1,0}, {1,1}, {1,2}, ...
        \endverbatim
        
        This is among the member functions that is not virtual
        because it is part of the innermost loop.
    */
    int change_box_coord(ubvector_int &boxt) {
      int j=dim-1;
      int ng=boxes;
      
      while (j >= 0) {

        boxt[j]=(boxt[j]+1) % ng;
        if (boxt[j] != 0) {
          return 1;
        }
        j--;

      }
      
      return 0;
    }
 
    /// Initialize grid
    virtual void init_grid(const vec_t &xl, const vec_t &xu, size_t ldim) {
      size_t j;
      vol=1.0;
      bins=1;
      
      for (j=0;j<ldim;j++) {
        double dx=xu[j]-xl[j];
        delx[j]=dx;
        vol *= dx;
        
        xi[j]=0.0;
        xi[dim+j]=1.0;
      }
      
      return;
    }

    /// Reset grid values
    virtual void reset_grid_values() {
      size_t i, j;
      
      for (i=0;i<bins;i++) {
        for (j=0;j<dim;j++) {
          d[i*dim+j]=0.0;
        }
      }
      return;
    }

    /** \brief Add the most recently generated result to the distribution

        This is among the member functions that is not virtual
        because it is part of the innermost loop.
    */
    void accumulate_distribution(ubvector_int &lbin, double y) {
      size_t j;
      
      for (j=0;j<dim;j++) {
        int i=lbin[j];
        d[i*dim+j]+=y;
      }
      return;
    }
    
    /** \brief Generate a random position in a given box
        
        Use the random number generator to return a random position x
        in a given box. The value of bin gives the bin location of the
        random position (there may be several bins within a given box)
        
        This is among the member functions that is not virtual
        because it is part of the innermost loop.
    */
    void random_point(vec_t &lx, ubvector_int &lbin, double &bin_vol,
                      const ubvector_int &lbox, const vec_t &xl, 
                      const vec_t &xu) {

      double lvol=1.0;

      size_t j;

      for (j=0;j<dim;++j) {

        // The equivalent of gsl_rng_uniform_pos()
        double rdn;
        do { 
          rdn=this->rng_dist(this->rng);
        } while (rdn==0);
        
        /* lbox[j] + ran gives the position in the box units, while z
           is the position in bin units.  */
        double z=((lbox[j]+rdn)/boxes)*bins;
        
        int k=(int)z;
        
        double y, bin_width;
        
        lbin[j]=k;

        if (k == 0) {
          bin_width=xi[dim+j];
          y=z*bin_width;
        } else {
          bin_width=xi[(k+1)*dim+j]-xi[k*dim+j];
          y=xi[k*dim+j]+(z-k)*bin_width;
        }
        
        lx[j]=xl[j]+y*delx[j];
        
        lvol *= bin_width;
      }
      
      bin_vol=lvol;

      return;
    }

    /// Resize the grid
    virtual void resize_grid(unsigned int lbins) {
      size_t j, k;

      /* weight is ratio of bin sizes */

      double pts_per_bin=(double) bins/(double) lbins;

      for (j=0;j<dim;j++) {
        double xold;
        double xnew=0;
        double dw=0;
        int i=1;

        for (k=1;k <= bins;k++) {
          dw+=1.0;
          xold=xnew;
          xnew=xi[k*dim+j];

          for (;dw > pts_per_bin;i++) {
            dw -= pts_per_bin;
            (xin[i])=xnew-(xnew-xold)*dw;
          }
        }

        for (k=1 ;k<lbins;k++) {
          xi[k*dim+j]=xin[k];
        }

        xi[lbins*dim+j]=1;
      }
      
      bins=lbins;

      return;
    }

    /// Refine the grid
    virtual void refine_grid() {

      size_t i, j, k;

      for (j=0;j<dim;j++) {
        double grid_tot_j, tot_weight;

        double oldg=d[j];
        double newg=d[dim+j];

        d[j]=(oldg+newg)/2;
        grid_tot_j=d[j];

        /* This implements gs[i][j]=(gs[i-1][j]+gs[i][j]+gs[i+1][j])/3 */

        for (i=1;i<bins-1;i++) {
          double rc=oldg+newg;
          oldg=newg;
          newg=d[(i+1)*dim+j];
          d[i*dim+j]=(rc+newg)/3;
          grid_tot_j+=d[i*dim+j];
        }
        d[(bins-1)*dim+j]=(newg+oldg)/2;

        grid_tot_j+=d[(bins-1)*dim+j];

        tot_weight=0;

        for (i=0;i<bins;i++) {
          weight[i]=0;

          if (d[i*dim+j] > 0) {
            oldg=grid_tot_j/d[i*dim+j];
            /* damped change */
            weight[i]=pow (((oldg-1)/oldg/log (oldg)), alpha);
          }
              
          tot_weight+=weight[i];
        }

        {
          double pts_per_bin=tot_weight/bins;

          double xold;
          double xnew=0;
          double dw=0;
          i=1;

          for (k=0;k<bins;k++) {
            dw+=weight[k];
            xold=xnew;
            xnew=xi[(k+1)*dim+j];

            for (;dw > pts_per_bin;i++) {
              dw -= pts_per_bin;
              (xin[i])=xnew-(xnew-xold)*dw/weight[k];
            }
          }

          for (k=1 ;k<bins ;k++) {
            xi[k*dim+j]=xin[k];
          }
            
          xi[bins*dim+j]=1;
        }
      }
      return;
    }

    /// Print limits of integration
    virtual void print_lim(const vec_t &xl, const vec_t &xu, 
                           unsigned long ldim) {
                   
      unsigned long j;
      
      (*outs) << "The limits of integration are:\n" << std::endl;
      for (j=0;j<ldim;++j) {
        (*outs) << "xl[" << j << "]=" << xl[j] << "    xu[" << j 
                <<  "]=" << xu[j] << std::endl;
      }
      (*outs) << std::endl;

      return;
    }

    /// Print header
    virtual void print_head(unsigned long num_dim, unsigned long calls,
                            unsigned int lit_num, unsigned int lbins, 
                            unsigned int lboxes)  {

      (*outs) << "num_dim=" << num_dim << ", calls=" << calls 
      << ", it_num=" << lit_num << ", max_it_num=" 
      << iterations << ", verb=" << this->verbose << ", alph=" << alpha 
      << ",\n mode=" << mode << ", boxes=" << lboxes
      << "\n\n    single.......iteration               "
      << "accumulated......results\n"
      << "iter   integral     sigma            integral   "
      << "  sigma        chi-sq/it" << std::endl;
      
      return;
    }
    
    /// Print results
    virtual void print_res(unsigned int itr, double res, 
                           double err, double cum_res, double cum_err,
                           double chi_sq) {
      (*outs).precision(5);
      (*outs) << itr << "     ";
      outs->setf(std::ios::showpos);
      (*outs) << res << "  ";
      outs->unsetf(std::ios::showpos);
      (*outs) << err << "     ";
      outs->setf(std::ios::showpos);
      (*outs) << cum_res << "  ";
      outs->unsetf(std::ios::showpos);
      (*outs) << cum_err << "  " << chi_sq << std::endl;
      (*outs).precision(6);
      return;
    }
 
    /// Print distribution
    virtual void print_dist(unsigned long ldim) {
      unsigned long i, j;

      if (this->verbose<2) {
        return;
      }
      
      for (j=0;j<ldim;++j) {
        (*outs) << "\n Axis " << j << std::endl;
        (*outs) << "       x     g" << std::endl;
        outs->setf(std::ios::showpos);
        for (i=0;i<bins;i++) {
          (*outs) << "weight [ " << (xi[(i)*dim+(j)]) << " , " 
                  << xi[(i+1)*dim+j] << " ] = ";
          (*outs) << " " << (d[(i)*dim+(j)]) << std::endl;
        }
        outs->unsetf(std::ios::showpos);
        (*outs) << std::endl;
      }
      (*outs) << std::endl;
      return;
    }
 
    /// Print grid
    virtual void print_grid(unsigned long ldim) {

      if (this->verbose<2) {
        return;
      }
      
      unsigned long i, j;
      for (j=0;j<ldim;++j) {
        (*outs) << "\n Axis " << j << std::endl;
        (*outs) << "      x     " << std::endl;
        outs->setf(std::ios::showpos);
        for (i=0;i <= bins;i++)  {
          (*outs) << (xi[(i)*dim+(j)]) << " ";
          if (i % 5 == 4) {
            (*outs) << std::endl;
          }
        }
        (*outs) << std::endl;
        outs->unsetf(std::ios::showpos);
      }
      (*outs) << std::endl;
      return;
    }

    /// Point for function evaluation
    vec_t x;

#endif

    public:

    mcarlo_vegas() {
      this->verbose=0;
      outs=&std::cout;
      alpha=1.5;
      iterations=5;
      mode=mode_importance;
      chisq=0;
      bins=bins_max;
      dim=0;
    }
    
    /// Allocate memory
    virtual int allocate(size_t ldim) {

      delx.resize(ldim);
      d.resize(bins_max*ldim);
      xi.resize((bins_max+1)*ldim);
      xin.resize(bins_max+1);
      weight.resize(bins_max);
      box.resize(ldim);
      bin.resize(ldim);
      x.resize(ldim);

      dim=ldim;

      return 0;
    }

    /** \brief Integrate function \c func from x=a to x=b.

        Original documentation from GSL:
        
        Normally, <tt>stage = 0</tt> which begins with a new uniform
        grid and empty weighted average. Calling vegas with <tt>stage
        = 1</tt> retains the grid from the previous run but discards
        the weighted average, so that one can "tune" the grid using a
        relatively small number of points and then do a large run with
        <tt>stage = 1</tt> on the optimized grid. Setting <tt>stage =
        2</tt> keeps the grid and the weighted average from the
        previous run, but may increase (or decrease) the number of
        histogram bins in the grid depending on the number of calls
        available. Choosing <tt>stage = 3</tt> enters at the main
        loop, so that nothing is changed, and is equivalent to
        performing additional iterations in a previous call.

        \verbatim embed:rst
        .. todo:: 

           Function mcarlo_vegas::vegas_minteg_err(): 

           - Should stage be passed by reference?
           - There was an update between gsl-1.12 and 1.15 which
             has not been implemented here yet.

        \endverbatim
    */
    virtual int vegas_minteg_err(int stage, func_t &func, size_t ndim, 
                                 const vec_t &xl, const vec_t &xu, 
                                 double &res, double &err) {

      size_t calls=this->n_points;

      double cum_int, cum_sig;
      size_t i, k, it;
        
      for (i=0;i<dim;i++) {
        if (xu[i] <= xl[i]) {
          std::string serr="Upper limit, "+dtos(xu[i])+", must be greater "+
            "than lower limit, "+dtos(xl[i])+", in mcarlo_vegas::"+
            "vegas_minteg_err().";
          O2SCL_ERR(serr.c_str(),exc_einval);
        }

        if (xu[i]-xl[i] > GSL_DBL_MAX) {
          O2SCL_ERR2("Range of integration is too large, please rescale ",
                         "in mcarlo_vegas::vegas_minteg_err().",exc_einval);
        }
      }

      if (stage == 0) {
        init_grid(xl,xu,dim);
        if (this->verbose>=1) {
          print_lim(xl,xu,dim);
        }
      }
      
      if (stage<=1) {
        wtd_int_sum=0;
        sum_wgts=0;
        chi_sum=0;
        it_num=1;
        samples=0;
        chisq=0;
      }
      
      if (stage <= 2) {

        unsigned int lbins=bins_max;
        unsigned int lboxes=1;

        if (mode != mode_importance_only) {

          /* shooting for 2 calls/box */
          
          // The original GSL code was:
          // boxes=floor (pow (calls/2.0, 1.0/dim));
          // but floor returns double on my machine, so 
          // we explicitly typecast here
          
          lboxes=((unsigned int)(floor(pow(calls/2.0,1.0/dim))));
          mode=mode_importance;

          if (2*lboxes >=  bins_max) {
            /* if bins/box < 2 */
            int box_per_bin=GSL_MAX(lboxes/bins_max,1);

            if (lboxes/box_per_bin<bins_max) lbins=lboxes/box_per_bin;
            else lbins=bins_max;
            lboxes=box_per_bin*lbins;

            mode=mode_stratified;
          }

        }

        //double tot_boxes=gsl_pow_int((double)boxes,dim);
        double tot_boxes=pow((double)lboxes,(double)dim);
        calls_per_box=((unsigned int)(GSL_MAX(calls/tot_boxes,2)));
        calls=((size_t)( calls_per_box*tot_boxes));
        
        /* total volume of x-space/(avg num of calls/bin) */
        jac=vol*pow((double) lbins, (double) dim)/calls;
         
        boxes=lboxes;

        /* If the number of bins changes from the previous invocation, bins
           are expanded or contracted accordingly, while preserving bin
           density */
        
        if (lbins!=bins) {
          resize_grid(lbins);
          if (this->verbose > 2) print_grid(dim);
        }
        if (this->verbose >= 1) {
          print_head(dim,calls,it_num,bins,boxes);
        }
      }

      it_start=it_num;

      cum_int=0.0;
      cum_sig=0.0;

      for (it=0;it<iterations;it++) {

        double intgrl=0.0, intgrl_sq=0.0;
        double tss=0.0;
        double wgt, var, sig;
        size_t lcalls_per_box=calls_per_box;
        double jacbin=jac;

        it_num=it_start+it;

        reset_grid_values();
        init_box_coord(box);

        do {
          volatile double m=0, q=0;
          double f_sq_sum=0.0;

          for (k=0;k<lcalls_per_box;k++) {
            double fval, bin_vol;
            
            random_point(x,bin,bin_vol,box,xl,xu);
            
            fval=func(dim,x);
            fval*=jacbin*bin_vol;

            /* recurrence for mean and variance (sum of squares) */

            {
              double dt=fval-m;
              m+=dt/(k+1.0);
              q+=dt*dt*(k/(k+1.0));
            }

            if (mode != mode_stratified) {
              double f_sq=fval*fval;
              accumulate_distribution(bin,f_sq);
            }
          }

          intgrl+=m*lcalls_per_box;

          f_sq_sum=q*lcalls_per_box;

          tss+=f_sq_sum;

          if (mode == mode_stratified) {
            accumulate_distribution (bin, f_sq_sum);
          }

        } while (change_box_coord(box));

        /* Compute final results for this iteration   */
        
        var=tss/(lcalls_per_box-1.0);
        
        if (var>0) {
          wgt=1.0/var;
        } else if (sum_wgts>0) {
          wgt=sum_wgts/samples;
        } else {
          wgt=0.0;
        }
        
        intgrl_sq=intgrl*intgrl;

        sig=sqrt(var);

        result=intgrl;
        sigma=sig;
        
        if (wgt > 0.0) {
          double lsum_wgts=sum_wgts;
          double m=(sum_wgts > 0) ? (wtd_int_sum/sum_wgts) : 0;
          double q=intgrl-m;
          
          samples++;
          sum_wgts+=wgt;
          wtd_int_sum+=intgrl*wgt;
          chi_sum+=intgrl_sq*wgt;

          cum_int= wtd_int_sum/sum_wgts;
          cum_sig=sqrt(1/sum_wgts);

          /* The original chisq formula from the Lepage paper is
             
             if ( samples > 1) {
             chisq=(chi_sum-wtd_int_sum*cum_int)/(samples-1.0);
             }
             
             This can suffer from cancellations and return a negative
             value of chi squared. We use the new formula from 
             GSL-1.12 instead
          */
          if (samples==1) {
            chisq=0;
          } else {
            chisq*=(samples-2.0);
            chisq+=(wgt/(1+(wgt/lsum_wgts)))*q*q;
            chisq/=(samples-1.0);
          }

        } else {
          cum_int+=(intgrl-cum_int)/(it+1.0);
          cum_sig=0.0;
        }         
        
        if (this->verbose >= 1) {
          print_res(it_num,intgrl,sig,cum_int,cum_sig,chisq);
          if (it+1 ==  iterations &&  this->verbose > 1) {
            print_grid(dim);
          }
        }

        if (this->verbose > 2) {
          print_dist(dim);
        }

        refine_grid ();

        if (this->verbose > 2) {
          print_grid(dim);
        }

      }

      /* 
         By setting stage to 1 further calls will generate independent
         estimates based on the same grid, although it may be rebinned. 
      */
      stage=1;

      res=cum_int;
      err=cum_sig;

      return GSL_SUCCESS;
    }

    virtual ~mcarlo_vegas() {}
  
    /// Integrate function \c func from x=a to x=b.
    virtual int minteg_err(func_t &func, size_t ndim, const vec_t &a, 
                           const vec_t &b, double &res, double &err) {
      allocate(ndim);
      chisq=0;
      bins=bins_max;
      int ret=vegas_minteg_err(0,func,ndim,a,b,res,err);
      return ret;
    }
    
    /** \brief Integrate function \c func over the hypercube from
        \f$ x_i=a_i \f$ to \f$ x_i=b_i \f$ for
        \f$ 0<i< \f$ ndim-1
    */
    virtual double minteg(func_t &func, size_t ndim, const vec_t &a, 
                          const vec_t &b) {
      double res;
      minteg_err(func,ndim,a,b,res,this->interror);
      return res;
    }
    
    /// Return string denoting type ("mcarlo_vegas")
    virtual const char *type() { return "mcarlo_vegas"; }
      
    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

