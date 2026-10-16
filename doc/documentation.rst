.. _documentation:

Documentation
=============

The complete **PyPL** reference is reported below.


Huang-Rhys solver
-----------------

.. automodule:: pypl.hr_solver
    :members:
    :undoc-members:
    :show-inheritance:

Lineshape
---------

.. automodule:: pypl.lineshape
    :members:
    :undoc-members:
    :show-inheritance:

Huang-Rhys factors
------------------

.. automodule:: pypl.hr_factors
    :members:
    :undoc-members:
    :show-inheritance:

1D configurational coordinate solver
------------------------------------

.. automodule:: pypl.config_coord_1d_solver
    :members:
    :undoc-members:
    :show-inheritance:

1D configurational coordinate solver: Franck-Condon routes
----------------------------------------------------------

Extends ``config_coord_1d_solver`` without changing it: quadrature of harmonic eigenfunctions, finite-difference
eigenstates of arbitrary (e.g. fitted, anharmonic) potentials, a first-order cubic anharmonic correction, and a
Franck-Condon lineshape with lines at the actual vibrational level energies.

.. automodule:: pypl.config_coord_1d_vof_solver
    :members:
    :undoc-members:
    :show-inheritance:

1D configurational coordinate potentials
----------------------------------------

Fits of a potential energy curve sampled along the configuration coordinate: parabola, polynomial, Morse, cubic and
splines (orders 2-5; natural, not-a-knot or clamped ends; harmonic or smooth C2 tails).

.. automodule:: pypl.config_coord_1d_potential
    :members:
    :undoc-members:
    :show-inheritance:

Intersystem crossing
--------------------

.. automodule:: pypl.isc
    :members:
    :undoc-members:
    :show-inheritance:

Utils
-----

.. automodule:: pypl.utils
    :members:
    :undoc-members:
    :show-inheritance:
