gsl_vector_const_view_array (const double * base, size_t n)
{
  _gsl_vector_const_view view = {{0, 0, 0, 0, 0}};
  
  if (n == 0)
    {
      do { gsl_error ("vector length n must be positive integer", "view_source.c", 28, GSL_EINVAL) ; return view ; } while (0);
      
    }
  {
    gsl_vector v = {0, 0, 0, 0, 0};
    
    v.data = (double *)base ;
    v.size = n;
    v.stride = 1;
    v.block = 0;
    v.owner = 0;
    ((_gsl_vector_view *)&view)->vector = v;
    
    return view;
  }
}

