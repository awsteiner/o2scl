/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include "acolm.h"

#include <o2scl/cloud_file.h>
#include <o2scl/vector_derint.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int o2scl_acol_pdma_get_base(void *vp, int &ndim, int &n,
			     double *&low, double *&high) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="prob_dens_mdim_amr") {
    return 1;
  }
  prob_dens_mdim_amr<> &pdma=amp->pdma_obj;
  ndim=pdma.n_dim;
  n=pdma.mesh.size();
  low=&(pdma.low[0]);
  high=&(pdma.high[0]);
  return 0;
}

int o2scl_acol_pdma_get_cube(void *vp, int ix, 
			     double *&low, double *&high,
			     double &frac_vol, double &weight) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  if (amp->type!="prob_dens_mdim_amr") {
    return 1;
  }
  prob_dens_mdim_amr<> &pdma=amp->pdma_obj;
  low=&(pdma.mesh[ix].low[0]);
  high=&(pdma.mesh[ix].high[0]);
  frac_vol=pdma.mesh[ix].frac_vol;
  weight=pdma.mesh[ix].weight;
  return 0;
}

int o2scl_acol_get_cli_options(void *vp, int &n, int *&sizes,
			       char *&chlist) {
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  n=amp->table_obj.get_nlines();
  amp->ctemp.clear();
  amp->itemp.clear();
  vector<string> options=amp->cl->get_option_list();
  for(size_t i=0;i<options.size();i++) {
    amp->itemp.push_back(options[i].length());
    for(size_t j=0;j<options[i].length();j++) {
      amp->ctemp.push_back(options[i][j]);
    }
  }
  n=options.size();
  sizes=&(amp->itemp[0]);
  chlist=&(amp->ctemp[0]);
  return 0;
}

int o2scl_acol_get_cli_options_type(void *vp, char *type,
				    int &n, int *&sizes,
				    char *&chlist) {

  string temp_type=type;
  
  o2scl_acol::acol_manager *amp=(o2scl_acol::acol_manager *)vp;
  n=amp->table_obj.get_nlines();
  amp->ctemp.clear();
  amp->itemp.clear();

  string cur_type=amp->type;
  
  amp->command_del(cur_type);

  if (temp_type.length()>0) {
    amp->command_add(temp_type);
  }
  
  vector<string> options=amp->cl->get_option_list();
  for(size_t i=0;i<options.size();i++) {
    amp->itemp.push_back(options[i].length());
    for(size_t j=0;j<options[i].length();j++) {
      amp->ctemp.push_back(options[i][j]);
    }
  }  
  n=options.size();
  sizes=&(amp->itemp[0]);
  chlist=&(amp->ctemp[0]);

  if (temp_type.length()>0) {
    amp->command_del(temp_type);
  }
  amp->command_add(cur_type);
  
  return 0;
}

