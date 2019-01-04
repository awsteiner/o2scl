/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_UNIFORM_GRID_H
#define O2SCL_UNIFORM_GRID_H

/** \file uniform_grid.h
    \brief File defining \ref o2scl::uniform_grid and its children
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>

// Forward definition of the uniform_grid class for HDF I/O
namespace o2scl {
  template<class data_t> class uniform_grid;
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  void hdf_input(hdf_file &hf, o2scl::uniform_grid<double> &t, 
		 std::string name);
  void hdf_output(hdf_file &hf, o2scl::uniform_grid<double> &t, 
		  std::string name);
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#ifdef O2SCL_NEVER_DEFINED
  // Forward definition of the uniform_grid class for HDF I/O
  namespace o2scl {
    template<class data_t> class uniform_grid;
  }
  
  // Forward definition of HDF I/O to extend friendship
  namespace o2scl_hdf { 
    class hdf_file; 
    void hdf_input(hdf_file &hf, o2scl::uniform_grid<double> &t, 
		   std::string name);
    void hdf_output(hdf_file &hf, o2scl::uniform_grid<double> &t, 
		    std::string name);
  }
#endif
  
  /** \brief A class representing a uniform linear or logarithmic grid
      
      \note This class has no public constructors and is to be
      instantiated through its children.

      This class should work for any floating-point type 
      compatible with std::pow() .

      Empty grids are those for which g_n_bins is zero. 

      The first and last bin are always exactly equal to the
      originally specified values of "start" and "end", but
      finite-precision errors may affect the inner grid points.

      \future Implement operator==, etc?

      \comment
      \future Add type() classes? 

      Actually type() classes may not be so useful cons
      \comment
  */
  template<class data_t=double> class uniform_grid {

  public:
  
  friend void o2scl_hdf::hdf_output
  (o2scl_hdf::hdf_file &hf, uniform_grid<double> &ug, std::string name);
  
  friend void o2scl_hdf::hdf_input
  (o2scl_hdf::hdf_file &hf, uniform_grid<double> &ug, std::string name);

  protected:
  
  /// The low-side of the first bin
  data_t g_start;
  /// The high-side of the last bin
  data_t g_end;
  /** \brief The width of each bin
	
      This should be always positive and non-zero for linear grids
      and always greater than 1 for logarithmic grids. 
  */
  data_t g_width;
  /// The number of bins
  size_t g_n_bins;
  /// If true, use a logarithmic scale
  bool g_log;

  /** \brief Construct a grid with specified values

      \note This function is not public because it might create
      grids that are non-sensical. We require users to create grid
      objects using one of the children which don't allow
      non-sensical grids.
  */
  uniform_grid(data_t start, data_t end, data_t width, size_t n_bins, 
	       bool log=false) {
    if (n_bins==0) {
      O2SCL_ERR2("Requested zero bins in ",
		 "uniform_grid::uniform_grid().",exc_einval);
    }
    if (start==end) {
      // This helps ensure (but no guarantee) so that one can 
      // use vector_bsearch() for a uniform_grid object.
      O2SCL_ERR2("Requested grid with start==end in ",
		 "uniform_grid::uniform_grid().",exc_einval);
    }
    if (log) {
      if (start>0.0) {
	if ((width<1.0 && start<end) ||
	    (width>1.0 && start>end)) {
	  O2SCL_ERR2("Incommensurate logarithmic grid in ",
		     "uniform_grid::uniform_grid().",exc_einval);
	}
      } else if (start<0.0) {
	if ((width<1.0 && start>end) ||
	    (width>1.0 && start<end)) {
	  O2SCL_ERR2("Incommensurate logarithmic grid in ",
		     "uniform_grid::uniform_grid().",exc_einval);
	}
      } else if (start==0.0 || end==0.0) {
	O2SCL_ERR2("Requested logarithmic grid with either ",
		   "start or end=0 in uniform_grid::uniform_grid().",
		   exc_einval);
      }
      if (width==1.0) {
	O2SCL_ERR2("Requested width 1 for logrithmic grid in ",
		   "uniform_grid::uniform_grid().",exc_einval);
      }
    } else {
      if ((width<0.0 && start<end) ||
	  (width>0.0 && start>end)) {
	O2SCL_ERR2("Incommensurate linear grid in ",
		   "uniform_grid::uniform_grid().",exc_einval);
      }
      if (width==0.0) {
	O2SCL_ERR2("Requested zero width for linear grid in ",
		   "uniform_grid::uniform_grid().",exc_einval);
      }
    }
    g_start=start;
    g_end=end;
    g_width=width;
    g_n_bins=n_bins;
    g_log=log;
  }

  public:

  /// Default constructor
  uniform_grid() {
    g_n_bins=0;
    g_start=0.0;
    g_end=0.0;
    g_width=0.0;
    g_log=false;
  }

  /** \brief Get the number of bins (regions in between grid points)

      This function returns zero if the grid is "empty".
  */
  size_t get_nbins() const {
    return g_n_bins;
  }

  /** \brief Get the number of points in the grid (always get_nbins()+1)
	
      This function will throw an exception if the grid is empty.
  */
  size_t get_npoints() const {
    if (g_n_bins==0) {
      O2SCL_ERR("Grid not set in uniform_grid::get_npoints().",
		exc_einval);
    }
    return g_n_bins+1;
  }

  /** \brief Return true if the grid is logarithmic
	
      This function will throw an exception if the grid is empty.
  */
  bool is_log() const {
    if (g_n_bins==0) {
      O2SCL_ERR("Grid not set in uniform_grid::get_npoints().",
		exc_einval);
    }
    return g_log;
  }

  /** \brief Get the first grid point
   */
  double get_start() {
    return g_start;
  }

  /** \brief Get the last grid point
   */
  double get_end() {
    return g_end;
  }
  
  /** \brief Get the interval between grid points
   */
  double get_width() {
    return g_width;
  }
  
  /** \brief Fill a vector with the specified grid

      If the vector is not big enough to hold the grid, it is
      automatically resized.

      This function will throw an exception if the grid is empty.
  */
  template<class resize_vec_t> void vector(resize_vec_t &v) const {
      
    if (g_n_bins==0) {
      O2SCL_ERR("Grid not set in uniform_grid::get_npoints().",
		exc_einval);
    }

    if (v.size()<g_n_bins+1) {
      v.resize(g_n_bins+1);
    }

    for(size_t i=0;i<g_n_bins+1;i++) {
      v[i]=(*this)[i];
    }
      
    return;
  }

  /** \brief Get the grid point with index \c i 
      (\f$ i \in [0,\mathrm{n_{bins}}] \f$)
  */
  const data_t operator[](size_t i) const {

    if (g_n_bins==0) {
      O2SCL_ERR("Grid not set in uniform_grid::get_npoints().",
		exc_einval);
    }

#if !O2SCL_NO_RANGE_CHECK
    if (i>g_n_bins) {
      std::string str=((std::string)"Index ")+o2scl::szttos(i)+
      " out of range "+o2scl::szttos(g_n_bins)+
      " in uniform_grid::operator[].";
      O2SCL_ERR(str.c_str(),exc_eindex);
    }
#endif

    // Linear case
    if (!g_log) {
      if (i==0) return g_start;
      else if (i==g_n_bins) return g_end;
      return g_start+i*g_width;
    } else {
      // Logarithmic case
      if (i==0) return g_start;
      else if (i==g_n_bins) return g_end;
      return g_start*std::pow(g_width,((data_t)i));
    }
  }

  /// Copy constructor
  uniform_grid(const uniform_grid &ug) {
    g_n_bins=ug.g_n_bins;
    if (ug.g_n_bins>0) {
      g_start=ug.g_start;
      g_end=ug.g_end;
      g_width=ug.g_width;
      g_log=ug.g_log;
    } else {
      g_start=0.0;
      g_end=0.0;
      g_width=0.0;
      g_log=false;
    }
  }
    
  /// Copy from = operator
  uniform_grid &operator=(const uniform_grid &ug) {
    g_n_bins=ug.g_n_bins;
    if (ug.g_n_bins>0) {
      g_start=ug.g_start;
      g_end=ug.g_end;
      g_width=ug.g_width;
      g_log=ug.g_log;
    } else {
      g_start=0.0;
      g_end=0.0;
      g_width=0.0;
      g_log=false;
    }
    return *this;
  }

  };

  /** \brief Linear grid with fixed number of bins and fixed endpoint
   */
  template<class data_t=double> class uniform_grid_end : 
    public uniform_grid<data_t> {
  public:

  /** \brief Create a grid with \c n_bins bins starting at 
      \c start and \c end

      The value of \c n_bins must be larger than zero and 
      \c start must not be the same as \c end. 
  */
  uniform_grid_end(data_t start, data_t end, size_t n_bins) : 
  uniform_grid<data_t>(start,end,(end-start)/((data_t)n_bins),
		       n_bins,false) {
  }
  };
  
  /** \brief Linear grid with fixed number of bins and fixed bin size
   */
  template<class data_t=double> class uniform_grid_width : 
    public uniform_grid<data_t> {
  public:

  /** \brief Create a grid with \c n_bins bins starting at 
      \c start with size \c width

      The value of \c n_bins must be larger than zero and 
      \c width must not be zero.
  */
  uniform_grid_width(data_t start, data_t width, size_t n_bins) :
  uniform_grid<data_t>(start,start+n_bins*width,width,n_bins,false) {
  }
  };

  /** \brief Linear grid with fixed endpoint and fixed bin size
   */
  template<class data_t=double> class uniform_grid_end_width : 
    public uniform_grid<data_t> {
  public:

  /** \brief Create a grid with bins of size \c width starting
      at \c start and ending at \c end

      The value of \c n_bins must be larger than zero and 
      \c start must not be the same as \c end. 
  */
  uniform_grid_end_width(data_t start, data_t end, data_t width) :
  uniform_grid<data_t>(start,end,width,
		       ((size_t)((end-start)/width)),false) {
  }
  };
  
  /** \brief Logarithmic grid with fixed number of bins and fixed endpoint
   */
  template<class data_t=double> class uniform_grid_log_end : 
    public uniform_grid<data_t> {
  public:
  
  /** \brief Create a logarithmic grid with \c n_bins bins starting at 
      \c start and \c end

      The value of \c n_bins must be larger than zero and 
      \c start must not be the same as \c end. 
  */
  uniform_grid_log_end(data_t start, data_t end, size_t n_bins) :
  uniform_grid<data_t>(start,end,std::pow(end/start,1.0/((data_t)n_bins)),
		       n_bins,true) {
  }
  };
  
  /** \brief Logarithmic grid with fixed number of bins and fixed bin size
   */
  template<class data_t=double> class uniform_grid_log_width : 
    public uniform_grid<data_t> {
  public:
  
  /** \brief Create a logarithmic grid with \c n_bins bins starting at 
      \c start with size \c width

      The value of \c n_bins must be larger than zero and 
      \c width must be greater than 1.
  */
  uniform_grid_log_width(data_t start, data_t width, size_t n_bins) :
  uniform_grid<data_t>(start,start*std::pow(width,n_bins),
		       width,n_bins,true) {
  }
  };
  
  /** \brief Logarithmic grid with fixed endpoint and fixed bin size
   */
  template<class data_t=double> class uniform_grid_log_end_width : 
    public uniform_grid<data_t> {
  public:
  
  /** \brief Create a logarithmic grid with bins of size \c width starting
      at \c start and ending at \c end

      The value of \c n_bins must be larger than zero and 
      \c start must not be the same as \c end. 
  */
  uniform_grid_log_end_width(data_t start, data_t end, data_t width) :
  uniform_grid<data_t>(start,end,width,
		       ((size_t)(log(end/start)/log(width))),true) {
  }
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
