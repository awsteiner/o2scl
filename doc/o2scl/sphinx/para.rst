:ref:`O2scl <o2scl>`

Parallel Programming with O2scl
===============================

Thread safety
-------------

Generally, O\ :sub:`2`\ scl objects are thread-safe in the same way that
classes like ``std::vector<double>`` are thread-safe:
reads are safe and writes are unsafe. It may be useful to make
objects ``const`` to ensure that one is reading data in a
thread-safe way. O\ :sub:`2`\ scl is designed to ensure const methods are
thread-safe, unless noted. (A few classes contain mutable internal
members which mean that const methods are not thread-safe, and
this is noted in the class documentation.)

OpenMP and O2scl
----------------

OpenMP support may be enabled during installation by
``--enable-openmp``. All OpenMP functionality is headers only, but
enabling OpenMP support during installation allows O\ :sub:`2`\ scl to
test the multithreaded behavior of :ref:`anneal_para <anneal_para>`,
:ref:`mcmc_para_base <mcmc_para_base>`, and :ref:`mcmc_para_table
<mcmc_para_table>`.

On some systems, code similar to the following may be required to
ensure that the error handler is valid on each OpenMP thread::
  
  int main(int argc, char *argv[]) {
    cout.setf(ios::scientific);
    MPI_Init(&argc,&argv);
    // Create a new error handler for this thread
    o2scl::err_hnd_cpp ee;
    o2scl::err_hnd=&ee;
    // Do stuff here
    MPI_Finalize();
    return 0;
  }

You can test to see if OpenMP support was enabled during installation
in the \ref o2scl::o2scl_settings object of type \ref
o2scl::lib_settings_class or with ``acol -v``.

MPI and O2scl
-------------

Currently, all MPI calls are in header classes, except for some
support for parallel HDF5 which is currently being developed and not
yet enabled. MPI functions are used in :ref:`anneal_para <anneal_para>`,
:ref:`mcmc_para_base <mcmc_para_base>`, and :ref:`mcmc_para_table
<mcmc_para_table>`.
    
You can test to see if MPI support was enabled during installation in
the \ref o2scl::o2scl_settings object of type \ref
o2scl::lib_settings_class or with ``acol -v``.

