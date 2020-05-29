Finite-temperature Equation of State Tables
===========================================

There are several classes designed to provide a consistent interface
to several EOS tables intended for core-collapse supernovae and
neutron-star mergers. The abstract base class is :ref:`eos_sn_base
<eos_sn_base>`. The child classes correspond to different EOS table
formats:

- :ref:`eos_sn_ls <eos_sn_ls>` - The Lattimer-Swesty EOS tables from Jim's
  webpage ([Lattimer91]_)
- :ref:`eos_sn_stos <eos_sn_stos>` - The H. Shen et al. EOS tables
  ([Shen98]_ and [Shen98b]_)
- :ref:`eos_sn_sht <eos_sn_sht>` - The G. Shen et al. EOS tables 
  ([Shen10a]_, [Shen10b]_, [Shen11]_)
- :ref:`eos_sn_hfsl <eos_sn_hfsl>` - The M. Hempel et al. EOS tables
  ([Hempel10]_ and [Hempel12]_)
- :ref:`eos_sn_oo <eos_sn_oo>` - The Lattimer-Swesty and H. Shen et al. tables
  reformatted by O'Connor and Ott ([OConnor10]_)

The :ref:`O2scl_eos <o2scle>` distribution does not contain the tables
themselves, as they are quite large and most are freely available.
:ref:`O2scl_eos <o2scle>` includes code which parses these tables and
puts them in \ref o2scl::tensor_grid3 objects for analysis by the
user.

The EOSs are stored in a set of :ref:`tensor_grid3
<o2scl:tensor_grid3>` objects on grids with baryon density in
:math:`\mathrm{fm}^{-3}`, electron fraction (unitless) and temperature
in :math:`\mathrm{MeV}`. The choice of baryon density is preferred to
that of 'rest mass density' (commonly denoted :math:`\rho`) because
the rest mass density has no unambiguous definition. This is
especially the case if dense matter contains deconfined quark matter.
On the other hand, baryon number is conserved (to a very good
approximation).

The electron fraction is also a natural variable to use when
tabulating, because charge is also conserved and the electron has
one unit of charge but no baryon number. The electron fraction is
also a sensible choice because of its central role in the relevant
weak interaction physics. This is also the choice made in 
[Lattimer91]_. Some tables are tabulated for constant "proton
fraction", which in this case includes all protons whether or not
they are inside nuclei. Because all baryonic matter in these EOS
tables is made up of neutrons and protons and because there are no
muons, one can assume that the electron fraction is equal to the
proton fraction. The total proton fraction is always larger than
the number fraction of free protons.

Not all tabulated EOSs contain all of the data, in which case the
associated :ref:`tensor_grid3 <tensor_grid3>` object may be empty. For
example, EOSs which do not contain the leptonic contributions do not
provide :cpp:var:`eos_sn_base::E`, :cpp:var:`eos_sn_base::F`,
:cpp:var:`eos_sn_base::S`, and :cpp:var:`o2scl::eos_sn_base::P`. In
these case, the grid is set for these objects but the data is set to
zero. To compute these from the data after loading the EOS table, use
:cpp:func:`o2scl::eos_sn_base::compute_eg()`.

Also, some EOS tables tabulate the 'mass fraction' of the 
various particles, but this is a slight misnomer. What is
actually tabulated is 'baryon number fraction', i.e. the
fraction of baryons which are in particles of type :math:`i`.
These fractions :math:`X_i` are defined by

.. math::
   
   X_i = A_i n_i n_B^{-1} \, ,

where :math:`A_i` is the number of baryons in particle :math:`i`
and :math:`n_i` is the number of particles per unit volume.
In the case of the representative heavy nucleus, the 
baryon number fraction is :math:`X_h = A n_h n_B^{-1}` where
:math:`A` is the baryon number of the representative heavy
nucleus in :cpp:var:`o2scl::eos_sn_base::A`.

The functions named ``load()`` in the children classes load
the entire EOS into memory. Memory allocation is automatically
performed, but not deallocated until ``free()`` or the destructor is
called.

After loading, you can interpolate the EOS by using 
:ref:`o2scl:tensor_grid3::interp_linear()` directly. For example,
the following returns the mass number at an arbitrary
baryon density, electron fraction, and temperature
assuming the table is stored in ``skm.dat``::

  ls_eos ls;
  ls.load("skm.dat");
  double nb=0.01, Ye=0.2, T=10.0;
  cout << ls.A.interp_linear(nb,Ye,T) << endl;

This function performs linear interpolation, however, some of the
grids are logarithmic, so linear interpolation on a logarithmic grid
leads to power-laws in between grid points. Note also that some grids
are not purely linear or purely logarithmic, but a mixture between the
two.

All of these classes are experimental.
