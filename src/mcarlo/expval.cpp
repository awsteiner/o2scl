/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2011-2023, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/expval.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

expval_base::expval_base(size_t n_blocks, size_t n_per_block) {
  if (n_blocks==0 || n_per_block==0) {
    O2SCL_ERR2("Value of zero for n_blocks or n_per_block in ",
	       "expval_base::expval_base().",exc_einval);
  }
  nblocks=n_blocks;
  nperblock=n_per_block;
  iblock=0;
  i=0;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

expval_base::~expval_base() {
}

expval_base::expval_base(const expval_base &ev) {
  nblocks=ev.nblocks;
  nperblock=ev.nperblock;
  iblock=ev.iblock;
  i=ev.i;
  name=ev.name;
  short_name=ev.short_name;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

expval_base &expval_base::operator=(const expval_base &ev) {
  if (&ev!=this) {
    nblocks=ev.nblocks;
    nperblock=ev.nperblock;
    iblock=ev.iblock;
    i=ev.i;
    name=ev.name;
    short_name=ev.short_name;
  }
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  return *this;
}

void expval_base::set_blocks(size_t n_blocks, size_t n_per_block) {
  if (n_blocks==0 || n_per_block==0) {
    O2SCL_ERR2("Value of zero for n_blocks or n_per_block in ",
	       "expval_base::set_blocks().",exc_einval);
  }
  if (i!=0 || iblock!=0) {
    O2SCL_ERR2("Cannot reblock before free() in ",
	       "expval_base::set_blocks().",exc_einval);
  }
  nblocks=n_blocks;
  nperblock=n_per_block;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
  return;
}

void expval_base::get_blocks(size_t &n_blocks, size_t &n_per_block) const {
  n_blocks=nblocks;
  n_per_block=nperblock;
  return;
}

void expval_base::free() {
  i=0;
  iblock=0;
  return;
}

void expval_base::get_block_indices(size_t &i_block, 
				   size_t &i_curr_block) const {
  i_block=iblock;
  i_curr_block=i;
  return;
}

bool expval_base::finished() const {
  if (iblock>=nblocks) return true;
  return false;
}

double expval_base::progress() const {
  return ((double)(iblock*nperblock+i))/((double)nblocks)/
    ((double)(nperblock));
}

expval_scalar::expval_scalar(size_t n_blocks, size_t n_per_block) : 
  expval_base(n_blocks,n_per_block) {
  vals.resize(n_blocks);
  for(size_t ii=0;ii<n_blocks;ii++) vals[ii]=0.0;
  current=0.0;
#if !O2SCL_NO_RANGE_CHECK
  is_valid();
#endif
}

expval_scalar::~expval_scalar() {
  free();
}

expval_scalar::expval_scalar(const expval_scalar &ev) : expval_base() {
  nblocks=ev.nblocks;
  nperblock=ev.nperblock;
  iblock=ev.iblock;
  i=ev.i;
  name=ev.name;
  short_name=ev.short_name;
  current=ev.current;
  vals=ev.vals;
}

expval_scalar &expval_scalar::operator=(const expval_scalar &ev) {
  if (&ev!=this) {
    nblocks=ev.nblocks;
    nperblock=ev.nperblock;
    iblock=ev.iblock;
    i=ev.i;
    name=ev.name;
    short_name=ev.short_name;
    current=ev.current;
    vals=ev.vals;
  }
  return *this;
}

void expval_scalar::set_blocks(size_t n_blocks, size_t n_per_block) {
  // The parent calls the error handler if reset() has not been
  // called and then calls free()
  expval_base::set_blocks(n_blocks,n_per_block);

  // Allocate vectors
  vals.resize(nblocks);
  for(size_t ii=0;ii<vals.size();ii++) vals[ii]=0.0;
  return;
}

void expval_scalar::free() {
  expval_base::free();
  vals.clear();
  return;
}

void expval_scalar::add(double val) {

  // If all blocks are full
  if (iblock==nblocks) {

    if (current!=0.0 || i!=0) {
      O2SCL_ERR2("Current or 'i' nonzero with full blocks in ",
		 "expval_scalar::add()",exc_esanity);
    }
    
    // Double up the data
    for(size_t j=0;j<nblocks/2;j++) {
      vals[j]=(vals[2*j]+vals[2*j+1])/2.0;
    }
    // If the number of blocks is even
    if (nblocks%2==0) {
      // Just readjust iblock
      iblock=nblocks/2;
    } else {
      // Take the odd block from vals and move it to current
      current=vals[nblocks-1];
      i=nperblock;
      iblock=nblocks/2;
    }
    // Clear out last half of 'vals'
    for(size_t j=nblocks/2;j<nblocks;j++) {
      vals[j]=0.0;
    }
    // Double nperblock
    nperblock*=2;

  }

  // Keep track of the rolling average and increment the index
  current+=(val-current)/(i+1);
  i++;

  // If the block is full
  if (i==nperblock) {
    
    // Store in vals and clear out current
    vals[iblock]=current;
    current=0.0;
    iblock++;
    i=0;
    
  }

  return;
}

void expval_scalar::current_avg_stats(double &avg, double &std_dev, 
				      double &avg_err, size_t &m_block,
				      size_t &m_per_block) const {
  
  // Only one block that is partially full
  if (iblock==0) {

    if (i==0) {
      O2SCL_ERR("No data in expval_scalar::current_avg_stats().",
		exc_efailed);
    } else {
      m_block=1;
      m_per_block=i;
      avg=current;
      std_dev=0.0;
      avg_err=0.0;
    }

  } else if (iblock==1) {

    // We're blocking, but have only one full block so far
    m_block=1;
    m_per_block=nperblock;
    avg=vals[0];
    std_dev=0.0;
    avg_err=0.0;

  } else if (iblock<=nblocks) {
    
    // Report the average and std. dev.
    // for all the blocks which have been finished
    m_block=iblock;
    m_per_block=nperblock;
    avg=vector_mean(iblock,vals);
    std_dev=vector_stddev(iblock,vals);
    avg_err=std_dev/sqrt(((double)iblock));

  }
      
  return;
}

void expval_scalar::current_avg(double &avg, double &std_dev, 
			    double &avg_err) const {
  size_t m_block, m_per_block;
  return current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);
} 

void expval_scalar::reblock_avg_stats(size_t new_blocks, double &avg, 
				  double &std_dev, double &avg_err,
				  size_t &m_per_block) const {

  if (new_blocks==0) {
    O2SCL_ERR2("Requested zero blocks in ",
	       "expval_scalar::reblock_avg_stats().",exc_einval);
  }
  
  ubvector dat(new_blocks);
  for(size_t ii=0;ii<dat.size();ii++) dat[ii]=0.0;
  
  // The ratio of the old to new block size
  size_t fact=iblock/new_blocks;
  if (fact==0) {
    O2SCL_ERR2("Not enough data for reblocking ",
	       "in expval_scalar::reblock_avg_stats().",exc_einval);
  }
  size_t iblock2=0;
  // Compute the sum
  for(size_t k=0;k<new_blocks;k++) {
    for(size_t j=0;j<fact;j++) {
      dat[k]+=vals[iblock2];
      iblock2++;
    }
    // Divide to get averages
    dat[k]/=((double)fact);
  }
  // Compute average
  avg=vector_mean(new_blocks,dat);
  // Compute std. dev. and avg. err. if available
  if (new_blocks>1) {
    std_dev=vector_stddev(new_blocks,dat);
    avg_err=std_dev/sqrt(((double)new_blocks));
  } else {
    std_dev=0.0;
    avg_err=0.0;
  }
  // Compute m_per_block
  m_per_block=fact*nperblock;
  
  return;
}

void expval_scalar::reblock_avg(size_t new_blocks, double &avg, 
			    double &std_dev, double &avg_err) const {
  size_t m_per_block;
  return reblock_avg_stats(new_blocks,avg,std_dev,avg_err,m_per_block);
}

const ubvector &expval_scalar::get_data() const {
  return vals;
}

const ubmatrix &expval_vector::get_data() const {
  return vals;
}

const tensor3<> &expval_matrix::get_data() const {
  return vals;
}

const double &expval_scalar::operator[](size_t i_block) const {
  return vals[i_block];
}

double &expval_scalar::operator[](size_t i_block) {
  return vals[i_block];
}

expval_vector::expval_vector() : expval_base(1,1) {
  nvec=0;
}

expval_vector::expval_vector(size_t n, size_t n_blocks, size_t n_per_block) : 
  expval_base(n_blocks,n_per_block) {
  nvec=n;
  if (nvec>0) {
    current.resize(nvec);
    for(size_t ii=0;ii<current.size();ii++) current[ii]=0.0;
    vals.resize(nvec,this->nblocks);
    for(size_t ii=0;ii<vals.size1();ii++) {
      for(size_t jj=0;jj<vals.size2();jj++) {
	vals(ii,jj)=0.0;
      }
    }
  }
}

expval_vector::~expval_vector() {
  if (nvec>0) {
    free();
  }
}

expval_vector::expval_vector(const expval_vector &ev) : expval_base() {
  nblocks=ev.nblocks;
  nperblock=ev.nperblock;
  iblock=ev.iblock;
  i=ev.i;
  name=ev.name;
  short_name=ev.short_name;
  nvec=ev.nvec;
  current=ev.current;
  vals=ev.vals;
}

expval_vector &expval_vector::operator=(const expval_vector &ev) {
  if (&ev!=this) {
    nblocks=ev.nblocks;
    nperblock=ev.nperblock;
    iblock=ev.iblock;
    i=ev.i;
    name=ev.name;
    short_name=ev.short_name;
    nvec=ev.nvec;
    current=ev.current;
    vals=ev.vals;
  }
  return *this;
}

void expval_vector::set_blocks_vec(size_t n, size_t n_blocks, size_t n_per_block) {

  // The parent calls the error handler if reset() has not been
  // called and then calls free()
  expval_base::set_blocks(n_blocks,n_per_block);

  if (n>0) {
    free();
    nvec=n;
    current.resize(nvec);
    for(size_t ii=0;ii<current.size();ii++) current[ii]=0.0;
    vals.resize(nvec,this->nblocks);
    for(size_t ii=0;ii<vals.size1();ii++) {
      for(size_t jj=0;jj<vals.size2();jj++) {
	vals(ii,jj)=0.0;
      }
    }
  }

  return;
}

void expval_vector::free() {
  if (nvec>0) {
    expval_base::free();
    vals.clear();
    current.clear();
  }
  nvec=0;
  return;
}

expval_matrix::expval_matrix() : expval_base(1,1) {
  nr=0;
  nc=0;
}

expval_matrix::expval_matrix
(size_t nrows, size_t ncols, size_t n_blocks, 
 size_t n_per_block) : expval_base(n_blocks,n_per_block) {

  if ((nrows==0 && ncols>0) || (ncols==0 && nrows>0)) {
    O2SCL_ERR2("Tried to use one finite index and one zero index in ",
	       "expval_matrix::expval_matrix().",exc_einval);
  }
  if (nrows>0 && ncols>0) {
    nr=nrows;
    nc=ncols;
    current.resize(nr,nc);
    for(size_t ii=0;ii<current.size1();ii++) {
      for(size_t jj=0;jj<current.size2();jj++) {
	current(ii,jj)=0.0;
      }
    }
    size_t dim[3]={nr,nc,this->nblocks};
    vals.resize(3,dim);
    // Set all values in vals to zero
    size_t tot=vals.total_size();
    for(size_t ii=0;ii<tot;ii++) {
      vals.unpack_index(ii,dim);
      vals.set(dim,0.0);
    }
  }
}

expval_matrix::~expval_matrix() {
  free();
}

expval_matrix::expval_matrix(const expval_matrix &ev) : expval_base() {
  nblocks=ev.nblocks;
  nperblock=ev.nperblock;
  iblock=ev.iblock;
  i=ev.i;
  name=ev.name;
  short_name=ev.short_name;
  nr=ev.nr;
  nc=ev.nc;
  current=ev.current;
  vals=ev.vals;
}

expval_matrix &expval_matrix::operator=(const expval_matrix &ev) {
  if (&ev!=this) {
    nblocks=ev.nblocks;
    nperblock=ev.nperblock;
    iblock=ev.iblock;
    i=ev.i;
    name=ev.name;
    short_name=ev.short_name;
    nr=ev.nr;
    nc=ev.nc;
    current=ev.current;
    vals=ev.vals;
  }
  return *this;
}

void expval_matrix::set_blocks_mat(size_t nrows, size_t ncols, size_t n_blocks, 
			   size_t n_per_block) {

  // The parent calls the error handler if reset() has not been
  // called and then calls free()
  expval_base::set_blocks(n_blocks,n_per_block);

  nr=nrows;
  nc=ncols;
  current.resize(nr,nc);
  for(size_t ii=0;ii<current.size1();ii++) {
    for(size_t jj=0;jj<current.size2();jj++) {
      current(ii,jj)=0.0;
    }
  }
  size_t dim[3]={nr,nc,this->nblocks};
  vals.resize(3,dim);
  // Set all values in vals to zero
  size_t tot=vals.total_size();
  for(size_t ii=0;ii<tot;ii++) {
    vals.unpack_index(ii,dim);
    vals.set(dim,0.0);
  }
  return;
}

void expval_matrix::free() {
  if (nr>0 || nc>0) {
    expval_base::free();
    std::vector<size_t> tmp;
    vals.resize(0,tmp);
    current.clear();
  }
  nr=0;
  nc=0;
  return;
}

void expval_scalar::is_valid() const {
  if (vals.size()!=nblocks) {
    O2SCL_ERR2("Vector size not equal to nblocks in ",
	       "expval_scalar::is_valid().",exc_efailed);
  }
  expval_base::is_valid();
}

