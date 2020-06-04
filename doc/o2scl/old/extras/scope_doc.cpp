// How not to provide subobjects to O2scl classes

void set_integrator(multi_inte<int> &mi) {
  gsl_inte_qag<size_t> one_dim_it[2];
  mi.set_oned_inte(one_dim_it[0],0);
  mi.set_oned_inte(one_dim_it[1],1);
}

void function() {
  composite_inte<int> cit;
  set_integrator(cit);
  cit.minteg_err(func,2,a,b,pa,res,err);
}
