*************************************************************************
RPMDrate - Bimolecular reaction rates via ring polymer molecular dynamics
*************************************************************************

About RPMDrate
==============

RPMDrate is a free, open-source software package for using ring polymer
molecular dynamics simulations to compute properties of chemical systems,
including, but not limited to, chemical reaction rates.

License
=======

RPMDrate is distributed under the terms of the 
`MIT license <http://www.opensource.org/licenses/mit-license>`_::

    Copyright (c) 2013 by Yury V. Suleimanov (ysuleyma@mit.edu, ysuleyma@princeton.edu)
                          Joshua W. Allen (jwallen@mit.edu)
                          William H. Green (whgreen@mit.edu)                         
    
    Permission is hereby granted, free of charge, to any person obtaining a 
    copy of this software and associated documentation files (the "Software"), 
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    and/or sell copies of the Software, and to permit persons to whom the 
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
    DEALINGS IN THE SOFTWARE. 

Installation
============

Dependencies
------------

RPMDrate depends on several other packages in order to provide its full
functional capabilities. The following dependencies are required to compile
and run RPMD:

* `Python <http://www.python.org/>`_ (version 2.5.x or later, including any version of Python 3, is recommended)

* `NumPy <http://numpy.scipy.org/>`_ (version 1.5.0 or later is recommended)

* `FFTW <http://www.fftw.org/>`_ (version 3.3 or later is recommended)

* A standard Fortran 90/95 compiler (e.g. ``gfortran``, ``g95``, ``ifort``, etc.)

RPMDrate uses the ``f2py`` tool (bundled with NumPy) to expose its Fortran code
to Python.

Getting RPMDrate
----------------

The best way to obtain a copy of the repository is to clone it using `git
<http://git-scm.com/>`_::

    $ git clone git://github.com/GreenGroup/RPMDrate.git

This enables you to easy update your local clone with the latest changes. If
you intend to contribute, you should fork the project on GitHub first, and
clone from your fork.

Compiling from Source
---------------------

RPMDrate is installed by invoking the ``setup.py`` script::

    $ python setup.py install

This will compile the Fortran source and Python wrapper code, which may take
some time. Note that you may need administrator privileges to install RPMD.

If you wish to use RPMDrate without installing, simply add the folder containing
this file to your ``PYTHONPATH`` environment variable and compile the source
code in-place using the command ::

    $ python setup.py build_ext --inplace

A Makefile that wraps these commands has been provided. The Makefile also
provides a clean target for deleting all files created during compilation.

List of the program files
---------------------
rpmdrate/_main.f90 - The main Fortran module providing the core RPMD algorithm
rpmdrate/_main.pyf - f2py wrapper script for exposing Fortran subroutines to Python
rpmdrate/_math.f90 - Auxiliary mathematics subroutines
rpmdrate/_surface.f90 - Subroutines defining the reactant and transition state dividing surfaces (including gradients and Hessians)
rpmdrate/__init__.py - Indicates that the directory is a Python package
rpmdrate/blas_lapack.f90 - Subroutines from BLAS and LAPACK (to remove external dependency on these programs)
rpmdrate/constants.py - Definitions of relevant physical constants
rpmdrate/element.py - Definitions of relevant chemical elements and their atomic masses
rpmdrate/input.py - Functionality for parsing RPMDrate input files
rpmdrate/interpolate.py - Classes used for linear, semilogarithmic, and logarithmic interpolation
rpmdrate/main.py - The main Python module used for running RPMD jobs
rpmdrate/quantity.py - Functions for manipulating physical quantities, including unit conversions
rpmdrate/surface.py - Python module that wraps the _surface.f90 Fortran module
rpmdrate/thermostat.py - Classes representing the available thermostats
rpmdrate.py - Main execution script for RPMDrate, used to start RPMD calculations
