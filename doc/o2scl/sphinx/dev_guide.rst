Developer guidelines
====================
    
This set of guidelines is intended for O2scl developers.

Release procedure
-----------------

- Make sure version numbers are updated in

  * ``configure.ac``
  * ``doc/o2scl/doxyfile``
  * ``doc/o2scl/part/doxyfile``
  * ``doc/o2scl/eos/doxyfile``
  * ``doc/o2scl/sphinx/conf.py``
  * ``doc/o2scl/part/sphinx/conf.py``
  * ``doc/o2scl/eos/sphinx/conf.py``
  * ``snap/snapcraft.yaml``

  and the tables in ``doc/o2scl/sphinx/download.rst`` are updated with
  "not yet set" for the new version hashes.
- For most recent commit, make sure the tests, examples, and
  documentation all succeed. Check source installs on ubuntu and
  OS X and a homebrew HEAD install.
- Try a test build on travis-ci.org
- Make sure snaps are working using 
  ``sudo snapcraft -d cleanbuild``
- Update NEWS file with recent changes.
- Make the final commit targeted for release. 
- Check the commit succeeds on travis-ci.org
- Promote the snaps on snapcraft from edge to beta
  at https://snapcraft.io/o2scl/releases .
- Refresh the documentation using ``make o2scl-doc``.
- Create the new distribution and copy to internal svn repo.
- Create github release, tagging the recent commit and uploading
  the distribution.
- Compute hashes with ``md5sum`` and ``openssl dgst -sha256``
  and update dl_page.dox with hashes and github
  release hash. Regenerate documentation with ``make o2scl-doc``.
- Do a ``make install`` on isospin so that ``make utk-sync-doc``
  can copy docs from the post-installation directory
- Copy the new distribution and the new sha256 hash to 
  https://isospin.roam.utk.edu/public_data/o2scl_dists/
- Update homebrew recipe with the new version number and new hash.
- Check installation using homebrew directly.
- Turn on build pushes in travis-ci.org and commit again since
  ``doc/o2scl/sphinx/download.rst`` has changed. 

Procedure for moving to new development version
-----------------------------------------------

- Update to the new development version number and new OLIB numbers
  in ``configure.ac``.
- Update version numbers in doxyfile files and main.dox files
  and in dl_page.dox .
- Update local configure scripts to refer to new version number
  if necessary
    
o2sclpy release procedure
-------------------------

- Update version numbers in ``o2sclpy/__init.py``, 
  ``doc/conf.py``, ``setup.py``, ``snap/snapcraft.yaml`` and ``doc/index.rst``
- Regenerate the o2sclpy documentation using ``make doc``
  and upload it to web using ``make sync-doc``
- Remove old dists in o2sclpy by clearing o2sclpy/dist directory
- Run ``python3 setup.py sdist bdist_wheel``
- Upload a new version of o2sclpy to pypi using
  ``twine upload dist/*``

Coding recommendations and guidelines
-------------------------------------

- Comment code liberally. 
- Avoid goto statements.
- Avoid hard-coded numerical values which are not either:
  (i) well-commented, or (ii) parameters modifiable by a class-user.
- Match the existing coding style when possible.
- Make sure all lines fit in 78 columns.
- Because this is a library (and not a end-user executable),
  it is generally better to call the error handler when an input is 
  incorrectly specified instead of assuming some alternative
  correct value.
- All new classes deserve new documentation and new testing code.
- Documentation and code should be in standard American English.
- All code must be released in GPLv3.
- Global variables and static member data should be avoided.
- When reasonable, put input parameters first and output
  parameters last. 
- When possible, templated vector parameters with size_t arguments
  should appear similar to ``size_t &n, vec_t &v``, and in that 
  order.
- All code should be ANSI-compatible, and, inasmuch as is 
  possible, operating system and platform independent.
- Exception messages should be as informative as possible.
- Exception messages should end with the full class and function name
  from where they are thrown.
- Functions which deallocate memory should never fail and should
  never be required to call the error handler. Also, class
  destructors should never be required to call the error handler.
- Avoid passing pointers. Pass either bare objects (which then
  require copy constructors), shared pointers, or const references.
- For numeric parameters to functions, if not all values
  are permissible, then the error handler should be called
  when a non-permissible value is given by the user and the
  function documentation should clearly explain why.
- Wherever possible, objects should be usable by default.
  Avoid zombie objects which are instantiated but unusable.
  An important exception to this rule is in classes which would
  otherwise have to load data from a file in order to operate
  normally, as file I/O should not occur in constructors
  (so that MPI can instantiate classes without worrying 
  about parallel I/O).
- Objects should thread-safe in the weak sense, that is, 
  different processes should be able to safely access and modify
  different instances of the same class at any time. Functions 
  which read (but not modify) class data should be thread-safe
  in the strong sense, that is, different processes should be
  able to read the same instance of a class at any time.
- Whereever possible, ensure your code compiles without
  warnings using flags analogous to the gcc string::

    -ansi -pedantic -Wno-long-long -Wall -Wno-unused -Wextra 
    -Wconversion -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings

- Avoid 'try' blocks, as a goal is that \o2 should compile
  with -fno-exceptions.
- Functions which return \c void should end with ``return;``.
- All functions which are called by the constructor should be
  documented as doing so

Documentation guidelines
------------------------

- Refer to other classes with \\ref if necessary. Refer
  to function parameters with \\c or embed them in html
  TT (text-type) commands.
- Bibliographic references should be used. When possible,
  include the DOI link which begins with the prefix 
  http://dx.doi.org (not the vendor-specific DOI link). 
- Comment Doxygen documentation with \\comment and \\endcomment.
  (Yes, sometimes comments in comments are useful.)

Git repository
--------------

- Communicate with the lead developer before, during, and after
  any non-trivial development. Communicate your ideas before
  development, so that you don't write many lines of code only to
  find that your pull request will be rejected. Communicate your
  ideas during development to avoid conflicting changes. Communicate
  your ideas after development to ensure they have a chance of being
  implmented. Subversion is not a replacement for real
  communication.
- Pull requests will be integrated into the trunk by the lead
  developer at whatever time they deem appropriate.
- Developer-specific files which are not platform-independent
  should not be added to the repository. Sometimes
  ``.gitignore`` can be used to ignore these files, but this
  should be done sparingly.
    
