Todo List
=========

:ref:`O2scl <o2scl>`

.. todo:: 

   Global todo list:

   - examples:

     * add ex_string, ex_tensor, and ex_eos_had_rmf to the documentation
   
   - integ/multip

     * replace the old inte_exp_sinh and inte_tanh_sinh integrators with
       integ_double_exp_boost_multip
     * fix the docs for subdivisions in inte_adapt_cern
     * check to see if polylog.h classes need new integrators (done, created
       new polylog_multip)
     * make a most consistent interface for the multiprecision integrators
       and code reuse via, e.g. a parent class?
     * create multiprecision versions of the GSL integrators using the
       coefficients in the boost headers

   - part

     * improve and calibrate fermion_rel_ld and fermion_rel_cdf25
     * implement fermion_rel_ld in eos_lepton

   - etc

     * implement more code reuse in funct_multip_transform
     * move calc_utf8 and funct_strings classes down in the header
       file hierarchy so we can include polylogs and other related
       functions in calc_utf8

   - acol

     * document the cyl_bessel functions better
     * implement find for more types and fix the problem (what is it?) for
       size_t[]
     * sync read generic and output
     * more options for acol kde to-table, including the options to
       change the upper and lower x limits

.. todolist::
