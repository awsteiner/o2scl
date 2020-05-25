Nuclear structure in the Hartree approximation
==============================================

See class \ref o2scl::nucleus_rmf .

\section ex_nucleus_rmf_sect Nucleus in the Hartree approximation example

This example uses the NL3 EOS as implemented in \ref o2scl::eos_had_rmf and
computes the structure of \f$ ^{208}\mathrm{Pb} \f$ using \ref
o2scl::nucleus_rmf. The results are stored in \ref o2scl::table_units objects
and output to HDF files. It also computes the energy of
three unoccupied neutron and proton levels. 
    
\dontinclude ex_nucleus_rmf.cpp
\skip Example:
\until End of example

\image html ex_nuc_prof.png "Nucleon density distributions in Pb-208"
\comment
\image latex ex_nuc_prof.eps "Nucleon density distributions in Pb-208" width=9cm
\endcomment

Typical output:
\verbinclude ex_nucleus_rmf.scr
