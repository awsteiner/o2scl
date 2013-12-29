/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2013, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/hist_ev.h>

using namespace std;
using namespace o2scl;
#if O2SCL_HDF_SVAR
using namespace o2scl_hdf;
#endif

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

hist_ev::hist_ev(uniform_grid<double> g, size_t n_blocks,
		 size_t n_per_block) : expval_base(n_blocks,n_per_block) {

  h.set_bin_edges(g);
  hsize=h.size();

  vals.resize(hsize,nblocks);
  nperblock_bins.resize(hsize);
  iblock_bins.resize(hsize);
  i_bins.resize(hsize);
  
  for(size_t ii=0;ii<nperblock_bins.size();ii++) nperblock_bins[ii]=nperblock;
  for(size_t ii=0;ii<iblock_bins.size();ii++) iblock_bins[ii]=0.0;
  for(size_t ii=0;ii<i_bins.size();ii++) i_bins[ii]=0.0;
  for(size_t ii=0;ii<vals.size1();ii++) {
    for(size_t jj=0;jj<vals.size2();jj++) {
      vals(ii,jj)=0.0;
    }
  }
}

void hist_ev::set_grid_blocks(uniform_grid<double> g, size_t n_blocks,
			      size_t n_per_block) {

  if (hsize>0) {
    O2SCL_ERR2("Cannot call set_grid_blocks before free() in ",
	       "hist_ev::set_grid_blocks().",exc_einval);
  }

  expval_base::set_blocks(n_blocks,n_per_block);
  
  h.clear();
  h.set_bin_edges(g);

  hsize=h.size();

  vals.resize(hsize,nblocks);
  nperblock_bins.resize(hsize);
  iblock_bins.resize(hsize);
  i_bins.resize(hsize);
  
  for(size_t ii=0;ii<nperblock_bins.size();ii++) nperblock_bins[ii]=nperblock;
  for(size_t ii=0;ii<iblock_bins.size();ii++) iblock_bins[ii]=0.0;
  for(size_t ii=0;ii<i_bins.size();ii++) i_bins[ii]=0.0;
  for(size_t ii=0;ii<vals.size1();ii++) {
    for(size_t jj=0;jj<vals.size2();jj++) {
      vals(ii,jj)=0.0;
    }
  }

  return;
}

void hist_ev::add(double x, double val) {

  size_t ibin=h.get_bin_index(x);

  // If all the blocks in the current histogram bin are full
  if (iblock_bins[ibin]==nblocks) {

    // Double up the data
    for(size_t j=0;j<nblocks/2;j++) {
      vals(ibin,j)=(vals(ibin,2*j)+vals(ibin,2*j+1))/2.0;
    }
    // If the number of blocks is even
    if (nblocks%2==0) {
      // Just leave the current bin as is and clear out last half of 'vals'
      i_bins[ibin]=0;
      iblock_bins[ibin]=nblocks/2;
      for(size_t j=nblocks/2;j<nblocks;j++) {
        vals(ibin,j)=0.0;
      }
    } else {
      // Take the odd block from vals and move it to the current bin
      h[ibin]=vals(ibin,nblocks-1);
      i_bins[ibin]=nperblock_bins[ibin];
      iblock_bins[ibin]=nblocks/2;
      for(size_t j=nblocks/2;j<nblocks;j++) {
        vals(ibin,j)=0.0;
      }
    }
    // Double nperblock
    nperblock_bins[ibin]*=2;

  }

  h[ibin]+=(val-h[ibin])/(i+1);
  i++;

  // If the block is full
  if (i==nperblock_bins[ibin]) {
    
    // Store in vals and clear out the current bin
    vals(ibin,iblock)=h[ibin];
    iblock_bins[ibin]++;
    h[ibin]=0.0;
    i_bins[ibin]=0;

  }
    
  return;
}
    
void hist_ev::current_avg_stats
(ubvector &reps, ubvector &avg, ubvector &std_dev, 
 ubvector &avg_err, ubvector_int &m_block, ubvector_int &m_per_block) {
  
  // Reallocate all of the vectors
  reps.resize(hsize);
  avg.resize(hsize);
  std_dev.resize(hsize);
  avg_err.resize(hsize);
  m_block.resize(hsize);
  m_per_block.resize(hsize);
      
  for(size_t j=0;j<hsize;j++) {

    // Obtain the representative coordinate from the histogram class
    reps[j]=h.get_rep_i(j);

    // No full blocks yet
    
    if (iblock_bins[j]==0) {

      if (i_bins[j]==0) {

	avg[j]=0.0;
	std_dev[j]=0.0;
	avg_err[j]=0.0;
	m_block[j]=0;
	m_per_block[j]=0;

      } else {
      
	avg[j]=h[j];
	std_dev[j]=0.0;
	avg_err[j]=0.0;
	m_block[j]=1;
	m_per_block[j]=i_bins[j];

      }
      
    } else if (iblock_bins[j]==1) {
      
      // Only one full block

      avg[j]=vals(j,0);
      std_dev[j]=0.0;
      avg_err[j]=0.0;
      m_block[j]=1;
      m_per_block[j]=nperblock_bins[j];
      
    } else {
      
      typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
      ubmatrix_row ur(vals,j);
      avg[j]=vector_mean(iblock_bins[j],ur);
      std_dev[j]=vector_stddev(iblock_bins[j],ur);
      avg_err[j]=std_dev[j]/sqrt(((double)(iblock_bins[j])));
      m_block[j]=iblock_bins[j];
      m_per_block[j]=nperblock_bins[j];
      
    }
    
  }
  
  return;
}

void hist_ev::current_avg(ubvector &reps, ubvector &avg, ubvector &std_dev, 
			  ubvector &avg_err) {
  ubvector_int m_block, m_per_block;
  current_avg_stats(reps,avg,std_dev,avg_err,m_block,m_per_block);
  return;
}

hist_2d_ev::hist_2d_ev(uniform_grid<double> hxg, uniform_grid<double> hyg, 
		       size_t n_blocks, size_t n_per_block) : 
  expval_base(n_blocks,n_per_block) {
  
  h.set_bin_edges(hxg,hyg);
  hsize_x=h.size_x();
  hsize_y=h.size_y();
  size_t dim[3]={hsize_x,hsize_y,n_blocks};
  vals.allocate(3,dim);
  nperblock_bins.resize(hsize_x,hsize_y);
  iblock_bins.resize(hsize_x,hsize_y);
  i_bins.resize(hsize_x,hsize_y);
      
  for(size_t ii=0;ii<nperblock_bins.size1();ii++) {
    for(size_t jj=0;jj<nperblock_bins.size1();jj++) {
      nperblock_bins(ii,jj)=n_per_block;
    }
  }

  for(size_t ii=0;ii<iblock_bins.size1();ii++) {
    for(size_t jj=0;jj<iblock_bins.size1();jj++) {
      iblock_bins(ii,jj)=0;
    }
  }

  for(size_t ii=0;ii<i_bins.size1();ii++) {
    for(size_t jj=0;jj<i_bins.size1();jj++) {
      i_bins(ii,jj)=0;
    }
  }

  // Set all values in vals to zero
  size_t tot=vals.total_size();
  for(size_t ii=0;ii<tot;ii++) {
    vals.unpack_indices(ii,dim);
    vals.set(dim,0.0);
  }
}

void hist_2d_ev::add(double x, double y, double val) {

  size_t ibin_x, ibin_y;
  h.get_bin_indices(x,y,ibin_x,ibin_y);
      
  // If all the blocks in the current bin are full
  if (iblock_bins(ibin_x,ibin_y)==nblocks) {

    // Double up the data
    for(size_t j=0;j<nblocks/2;j++) {
      vals.get(ibin_x,ibin_y,j)=(vals.get(ibin_x,ibin_y,2*j)+
				 vals.get(ibin_x,ibin_y,2*j+1))/2.0;
    }
    // If the number of blocks is even
    if (nblocks%2==0) {
      // Just leave current as is and clear out last half of 'vals'
      i_bins(ibin_x,ibin_y)=0;
      iblock_bins(ibin_x,ibin_y)=nblocks/2;
      for(size_t j=nblocks/2;j<nblocks;j++) {
        vals.get(ibin_x,ibin_y,j)=0.0;
      }
    } else {
      // Take the odd block from vals and move it to current
      h.set_wgt_i(ibin_x,ibin_y,vals.get(ibin_x,ibin_y,nblocks-1));
      i_bins(ibin_x,ibin_y)=nperblock_bins(ibin_x,ibin_y);
      iblock_bins(ibin_x,ibin_y)=nblocks/2;
      for(size_t j=nblocks/2;j<nblocks;j++) {
        vals.get(ibin_x,ibin_y,j)=0.0;
      }
    }
    // Double nperblock
    nperblock_bins(ibin_x,ibin_y)*=2;

  }

  h.get_wgt_i(ibin_x,ibin_y)+=(val-h.get_wgt_i(ibin_x,ibin_y))/(i+1);
  i++;

  // If the block is full
  if (i==nperblock_bins(ibin_x,ibin_y)) {
    
    // Store in vals and clear out current
    vals.get(ibin_x,ibin_y,iblock)=h.get_wgt_i(ibin_x,ibin_y);
    iblock_bins(ibin_x,ibin_y)++;
    h.set_wgt_i(ibin_x,ibin_y,0.0);
    i_bins(ibin_x,ibin_y)=0;

  }

  return;
}

void hist_2d_ev::current_avg_stats
(ubvector &rep_x, ubvector &rep_y, ubmatrix &avg, ubmatrix &std_dev, 
 ubmatrix &avg_err, ubmatrix_int &m_block, ubmatrix_int &m_per_block) {
  
  // Reallocate all of the vectors
  rep_x.resize(hsize_x);
  rep_y.resize(hsize_y);
  avg.resize(hsize_x,hsize_y);
  std_dev.resize(hsize_x,hsize_y);
  avg_err.resize(hsize_x,hsize_y);
  m_block.resize(hsize_x,hsize_y);
  m_per_block.resize(hsize_x,hsize_y);

  // Get representatives from histogram class

  for(size_t j=0;j<hsize_x;j++) {
    rep_x[j]=h.get_x_rep_i(j);
  }
  for(size_t j=0;j<hsize_y;j++) {
    rep_y[j]=h.get_y_rep_i(j);
  }

  for(size_t j=0;j<hsize_x;j++) {
    for(size_t k=0;k<hsize_y;k++) {

      // No full blocks yet

      if (iblock_bins(j,k)==0) {

	if (i_bins(j,k)==0) {
	  
	  avg(j,k)=0.0;
	  std_dev(j,k)=0.0;
	  avg_err(j,k)=0.0;
	  m_block(j,k)=0;
	  m_per_block(j,k)=0;
	  
	} else {
	  
	  avg(j,k)=h.get_wgt_i(j,k);
	  std_dev(j,k)=0.0;
	  avg_err(j,k)=0.0;
	  m_block(j,k)=1;
	  m_per_block(j,k)=i_bins(j,k);

	}

      } else if (iblock_bins(j,k)==1) {

	// Only one full block
	
	avg(j,k)=vals.get(j,k,0);
	std_dev(j,k)=0.0;
	avg_err(j,k)=0.0;
	m_block(j,k)=1;
	m_per_block(j,k)=nperblock_bins(j,k);
	
      } else {

	// Generic case

	// Create a vector view which fixes all but the last index
	size_t index[3]={j,k,0};
	typedef boost::numeric::ublas::vector<double> ubvector;
	typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
	ubvector_slice vec=vals.vector_slice(2,index);

	// Compute the stats from 'vec'
	avg(j,k)=vector_mean(iblock_bins(j,k),vec);
	std_dev(j,k)=vector_stddev(iblock_bins(j,k),vec);
	avg_err(j,k)=std_dev(j,k)/sqrt(((double)(iblock_bins(j,k))));
	m_block(j,k)=iblock_bins(j,k);
	m_per_block(j,k)=nperblock_bins(j,k);

      }
    }
  }
  
  return;
}

#if O2SCL_HDF_SVAR

void hist_ev::hdf_output(hdf_file &hf, std::string hdf_name) {
  
  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","hist_ev");

  hf.set_szt("hsize",hsize);
  hf.setd_mat("vals",vals);
  hf.set_szt_vec("iblock_bins",iblock_bins);
  hf.set_szt_vec("nperblock_bins",nperblock_bins);
  hf.set_szt_vec("i_bins",i_bins);
  //h.hdf_output(hf,"hist");

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void hist_ev::hdf_input(hdf_file &hf, std::string hdf_name) {
  
  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_group_by_type(hf,"hist_ev",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type hist_ev found in ",
		     "o2scl_hdf::hdf_input().",exc_efailed);
    }
  }

  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Check typename
  string type;
  hf.gets_fixed("o2scl_type",type);
  if (type!="hist_ev") {
    O2SCL_ERR2("Typename in HDF group does not match ",
		   "class in hdf_input().",exc_einval);
  }

  hf.get_szt("hsize",hsize);
  hf.getd_mat("vals",vals);
  hf.get_szt_vec("iblock_bins",iblock_bins);
  hf.get_szt_vec("nperblock_bins",nperblock_bins);
  hf.get_szt_vec("i_bins",i_bins);
  //h.hdf_input(hf,"hist");

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

#ifdef O2SCL_NEVER_DEFINED

void hist_2d_ev::hdf_output(hdf_file &hf, std::string hdf_name) {
  
  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Add typename
  hf.sets_fixed("o2scl_type","hist_2d_ev");

  hf.set_szt("hsize_x",hsize_x);
  hf.set_szt("hsize_y",hsize_y);
  hf.setd_ten("vals",vals);
  hf.set_szt_mat("iblock_bins",iblock_bins);
  hf.set_szt_mat("nperblock_bins",nperblock_bins);
  hf.set_szt_mat("i_bins",i_bins);
  //h.hdf_output(hf,"hist");

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

void hist_2d_ev::hdf_input(hdf_file &hf, std::string hdf_name) {
  
  // If no name specified, find name of first group of specified type
  if (hdf_name.length()==0) {
    hf.find_group_by_type(hf,"hist_2d_ev",hdf_name);
    if (hdf_name.length()==0) {
      O2SCL_ERR2("No object of type hist_2d_ev found in ",
		     "o2scl_hdf::hdf_input().",exc_efailed);
    }
  }

  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(hdf_name);
  hf.set_current_id(group);

  // Check typename
  string type;
  hf.gets_fixed("o2scl_type",type);
  if (type!="hist_2d_ev") {
    O2SCL_ERR2("Typename in HDF group does not match ",
		   "class in hdf_input().",exc_einval);
  }

  hf.get_szt("hsize_x",hsize_x);
  hf.get_szt("hsize_y",hsize_y);
  hf.getd_ten("vals",vals);
  hf.get_szt_mat("iblock_bins",iblock_bins);
  hf.get_szt_mat("nperblock_bins",nperblock_bins);
  hf.get_szt_mat("i_bins",i_bins);
  //h.hdf_input(hf,"hist");

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return;
}

#endif

#endif
