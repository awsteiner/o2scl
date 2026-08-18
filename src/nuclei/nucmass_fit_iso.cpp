/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2026, Andrew W. Steiner

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
#include <o2scl/nucmass_fit_iso.h>
#include <o2scl/string_conv.h>

using namespace std;
using namespace o2scl;

nucmass_fit_iso::nucmass_fit_iso() {
  x=4;
  min_data_per_param=3;
  n_threads=1;
  verbose=0;
}

size_t nucmass_fit_iso::fit(nucmass_fit_base &proto, nucmass &src,
                             int minZ, int maxZ) {

  // Set the number of threads, updating n_threads to the actual
  // value OpenMP grants us (as in o2scl::anneal_para and
  // o2scl::diff_evo_para)
#ifdef O2SCL_SET_OPENMP
  omp_set_num_threads(n_threads);
#pragma omp parallel
  {
    n_threads=omp_get_num_threads();
  }
#else
  n_threads=1;
#endif

  // Each thread gets its own copy of mf, so that its dist and
  // def_mmin (which fit() below mutates) aren't shared across
  // chains fit concurrently. The user's own settings on mf (e.g.
  // mf.def_mmin.ntrial) are preserved by copying mf itself, not a
  // freshly-constructed nucmass_fit.
  std::vector<nucmass_fit> vmf(n_threads,mf);

  int nZ=maxZ-minZ+1;
  std::vector<std::shared_ptr<nucmass_fit_base> > results(nZ);

#ifdef O2SCL_SET_OPENMP
#pragma omp parallel for default(shared)
#endif
  for(int iz=0;iz<nZ;iz++) {

    int Z=minZ+iz;
    int zlo=Z-x;
    int zhi=Z+x;

#ifdef O2SCL_SET_OPENMP
    size_t ithread=omp_get_thread_num();
#else
    size_t ithread=0;
#endif
    
    nucmass_fit &local_mf=vmf[ithread];

    std::shared_ptr<nucmass_fit_base> f(proto.clone());

    std::string expr=((std::string)"Z>="+o2scl::itos(zlo)+
                       " && Z<="+o2scl::itos(zhi));
    nucdist_set(local_mf.dist,src,expr);

    if (verbose>0) {
      std::cout << "nucmass_fit_iso::fit(), Z=" << Z 
                << ", num. params=" << f->nfit
                << " num. nuclei=" << local_mf.dist.size() << " "
                << f->type() << std::endl;
    }
    
    if (local_mf.dist.size()<((size_t)min_data_per_param)*f->nfit) {
      O2SCL_ERR((((std::string)"Cannot fit ")+
                 o2scl::szttos(local_mf.dist.size())+" nuclei with "+
                 o2scl::szttos(f->nfit)+" parameters.").c_str(),
                o2scl::exc_einval);
      
      // Not enough calibration data for this chain; leave results[iz]
      // unset so that any previously-fit chain at this Z (if any) is
      // left untouched below.
      
      //continue;
    }
    
    double res;
    local_mf.fit(*f,res);

    if (local_mf.def_mmin.last_ntrial>=local_mf.def_mmin.ntrial) {
      // Did not converge within the trial budget; skip this chain
      // rather than store a possibly bad fit.
      continue;
    }

    results[iz]=f;
  }

  // Merge the per-Z results into chains serially, in the same order
  // the original sequential loop would have, so a chain which
  // wasn't (re-)fit this call leaves any previously-fit chain at
  // that Z untouched, exactly as before parallelization.
  size_t n_fit=0;
  for(int iz=0;iz<nZ;iz++) {
    if (results[iz]) {
      chains[minZ+iz]=results[iz];
      n_fit++;
    }
  }

  n=chains.size();

  return n_fit;
}

size_t nucmass_fit_iso::fit_interp_pass(nucmass &src) {

#ifdef O2SCL_SET_OPENMP
  omp_set_num_threads(n_threads);
#pragma omp parallel
  {
    n_threads=omp_get_num_threads();
  }
#else
  n_threads=1;
#endif

  // std::map iterators aren't randomly splittable across an OpenMP
  // "for" loop, so collect the (Z, chain) pairs into a vector first
  // -- this is a cheap serial pass over at most a few hundred
  // entries.
  std::vector<std::pair<int,std::shared_ptr<nucmass_fit_base> > >
    entries(chains.begin(),chains.end());
  size_t n_chains=entries.size();
  std::vector<int> results(n_chains);

#ifdef O2SCL_SET_OPENMP
#pragma omp parallel for default(shared)
#endif
  for(size_t j=0;j<n_chains;j++) {

    int Z=entries[j].first;
    int zlo=Z-x;
    int zhi=Z+x;

    // A separate, thread-local dist for every chain -- there is no
    // shared mutable state here between chains fit concurrently,
    // since each chain's fit_interp() (e.g. \ref
    // nucmass_two_interp::fit_interp()) builds its own local \ref
    // nucmass_fit internally rather than reusing \ref mf.
    std::vector<nucleus> dist;
    std::string expr=((std::string)"Z>="+o2scl::itos(zlo)+
                       " && Z<="+o2scl::itos(zhi));
    nucdist_set(dist,src,expr);

    if (verbose>0) {
      std::cout << "nucmass_fit_iso::fit_interp(), Z=" << Z
                << std::endl;
    }
    
    // The default nucmass_fit_base::fit_interp() is a no-op which
    // returns 1; only chains (e.g. nucmass_two_interp) which
    // override it to actually train a secondary stage return 0.
    results[j]=entries[j].second->fit_interp(dist);
  }

  size_t n_done=0;
  for(size_t j=0;j<n_chains;j++) {
    if (results[j]==0) n_done++;
  }

  return n_done;
}

bool nucmass_fit_iso::is_included(int Z, int N) {
  std::map<int,std::shared_ptr<nucmass_fit_base> >::iterator it=
    chains.find(Z);
  if (it==chains.end()) return false;
  return it->second->is_included(Z,N);
}

double nucmass_fit_iso::mass_excess(int Z, int N) {
  std::map<int,std::shared_ptr<nucmass_fit_base> >::iterator it=
    chains.find(Z);
  if (it==chains.end()) {
    O2SCL_ERR((((std::string)"No fitted chain for the requested Z, ")+
               o2scl::itos(Z)+", in nucmass_fit_iso::"+
               "mass_excess().").c_str(),exc_einval);
  }
  return it->second->mass_excess(Z,N);
}

int nucmass_fit_iso::hdf_output(o2scl_hdf::hdf_file &hf,
                                 std::string name) const {

  if (!hf.has_write_access()) {
    O2SCL_ERR2("File not opened with write access in ",
               "nucmass_fit_iso::hdf_output().",exc_efailed);
  }

  // Start group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Add typename and metadata
  hf.sets_fixed("o2scl_type","nucmass_fit_iso");
  hf.seti("x",x);
  hf.seti("min_data_per_param",min_data_per_param);

  // Sorted list of fitted Z values (std::map already iterates in
  // key order, since Z is the key)
  std::vector<int> z_list;
  for(std::map<int,std::shared_ptr<nucmass_fit_base> >::const_iterator
        it=chains.begin();it!=chains.end();it++) {
    z_list.push_back(it->first);
  }
  hf.seti_vec("z_list",z_list);

  // Each chain writes itself, polymorphically, into its own nested
  // subgroup. We're already positioned inside this object's group,
  // so each chain_<Z> subgroup nests correctly underneath it.
  for(std::map<int,std::shared_ptr<nucmass_fit_base> >::const_iterator
        it=chains.begin();it!=chains.end();it++) {
    std::string cname=((std::string)"chain_")+o2scl::itos(it->first);
    it->second->hdf_output(hf,cname);
  }

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return 0;
}

int nucmass_fit_iso::hdf_input_n(o2scl_hdf::hdf_file &hf,
                                  nucmass_fit_base &proto,
                                  std::string &name) {

  // If no name specified, find name of first group of specified type
  if (name.length()==0) {
    hf.find_object_by_type("nucmass_fit_iso",name);
    if (name.length()==0) {
      O2SCL_ERR2("No object of type nucmass_fit_iso found in ",
                 "nucmass_fit_iso::hdf_input_n().",exc_efailed);
    }
  }

  // Discard any chains already present
  clear();

  // Open main group
  hid_t top=hf.get_current_id();
  hid_t group=hf.open_group(name);
  hf.set_current_id(group);

  // Check typename
  std::string type2;
  hf.gets_fixed("o2scl_type",type2);
  if (type2!="nucmass_fit_iso") {
    O2SCL_ERR2("Typename in HDF group does not match class in ",
               "nucmass_fit_iso::hdf_input_n().",exc_einval);
  }

  // Load metadata
  hf.geti("x",x);
  hf.geti("min_data_per_param",min_data_per_param);

  std::vector<int> z_list;
  hf.geti_vec("z_list",z_list);

  // Reconstruct each chain: clone proto, then let the clone load
  // its own nested subgroup, whatever form that data takes.
  for(size_t i=0;i<z_list.size();i++) {
    std::shared_ptr<nucmass_fit_base> f(proto.clone());
    std::string cname=((std::string)"chain_")+o2scl::itos(z_list[i]);
    f->hdf_input(hf,cname);
    chains[z_list[i]]=f;
  }

  n=chains.size();

  // Close group
  hf.close_group(group);

  // Return location to previous value
  hf.set_current_id(top);

  return 0;
}

void nucmass_fit_iso::hdf_input(o2scl_hdf::hdf_file &hf,
                                 nucmass_fit_base &proto,
                                 std::string name) {
  hdf_input_n(hf,proto,name);
  return;
}
