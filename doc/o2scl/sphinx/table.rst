Data Tables
===========
    
:ref:`O2scl <o2scl>`

The class :ref:`table <table>` is a container to hold and perform
operations on related columns of data. It supports column operations,
interpolation, column reference by either name or index, binary
searching (in the case of ordered columns), and sorting. It also
supports creating new columns with functions which operate on column
names using the :ref:`calc_utf8 <calc_utf8>` class (see also
:ref:`Algebraic Function Evaluation`). A child class,
:ref:`table_units <table_units>` is similiar except that it
additionally allows one to specify physical units for each column and
convert units using :ref:`convert_units <convert_units>` .

The class :ref:`table3d <table3d>` is a generalization of :ref:`table
<table>` which operates on two-dimensional slices of data rather than
one-dimensional columns.

For higher-dimensional generalizations, the class :ref:`tensor_grid
<tensor_grid>` provides a set of data defined on a grid in an
arbitrary number of dimensions.

Table example
----------------

.. literalinclude:: ../../../examples/ex_table.cpp
   :language: c++		    
   :start-after: sphinx-example-start
