:ref:`O2scl <o2scl>`

Parallel Programming with O2scl
===============================

Thread safety
-------------

Generally, \o2 objects are thread-safe in the same way that
classes like <tt>std::vector&lt;double&gt;</tt> are thread-safe:
reads are safe and writes are unsafe. It may be useful to make
objects \c const to ensure that one is reading data in a
thread-safe way. \o2 is designed to ensure const methods are
thread-safe, unless noted. (A few classes contain mutable internal
members which mean that const methods are not thread-safe, and
this is noted in the class documentation.)

OpenMP and O2scl
----------------

OpenMP support may be enabled during installation by
<tt>--enable-openmp</tt>. All OpenMP functionality is headers
only, but enabling OpenMP support during installation allows \o2
to test the multithreaded behavior of \ref o2scl::anneal_para,
\ref o2scl::mcmc_para_base, and \ref o2scl::mcmc_para_table.

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
o2scl::lib_settings_class or with <tt>acol -v</tt>.

MPI and O2scl
-------------

Currently, all MPI calls are in header classes, except for some
support for parallel HDF5 which is currently being developed and not
yet enabled. MPI functions are used in \ref o2scl::anneal_para, \ref
o2scl::mcmc_para_base, and \ref o2scl::mcmc_para_table .
    
You can test to see if MPI support was enabled during installation in
the \ref o2scl::o2scl_settings object of type \ref
o2scl::lib_settings_class or with <tt>acol -v</tt>.

