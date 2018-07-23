----------------------------------------------
# gen.parRep 
----------------------------------------------

C++ implementation of the Generalized Parallel Replica algorithm

Link to this repository : 
  * https://gitlab.inria.fr/parallel-replica/gen.parRep

Link to the project gitlab page (includes more resources):
  * https://gitlab.inria.fr/parallel-replica

See the following articles:
  * A generalized parallel replica dynamics: Binder, Lelièvre & Simpson, 2015: https://doi.org/10.1016/j.jcp.2015.01.002
  * A new implementation of the Generalized Parallel Replica Dynamics targeting metastable biochemical systems: Hédin & Lelièvre, 2018: https://arxiv.org/abs/1807.02431

----------------------------------------------
## DOCUMENTATION
----------------------------------------------

Please consult the Gitlab wiki for an overview of this software use: 

https://gitlab.inria.fr/parallel-replica/gen.parRep/wikis/home

Doxygen can be used for generating a code reference documentation within the docs directory:
  * doxygen parRep.doxy 

By default an HTML reference is generated (ready to read), and a docs/latex directory is also created, although the pdf version requires to be manually compiled (pdfLaTeX, see the Makefile in ./docs/latex).

The Lua input files, in ./mol, used for starting computations, are self documented, and together
with the bash submission files in the ./run directory should be enough for starting to use the program

----------------------------------------------
## DEPENDANCIES
----------------------------------------------

You need an MPI development framework installed, compatible with the MPI 3.0 or newer standard : tested implementations : 
  * Open MPI 1.10.2 (from Ubuntu 16.04 repositories)
CMake will try to locate it automatically.

You will need to Download and/or Compile the OpenMM library, a required dependancy : 
  * see https://simtk.org/home/openmm
  * and/or https://github.com/pandegroup/openmm
  * tested with version 7.0 to 7.2 (manually compiled and installed)
You may need Nvidia CUDA or AMD OpenCL toolkit for enabling GPU acceleration ; see OpenMM documentation.
You may also need to edit CMakeLists.txt for specifying path to the include and lib directories of OpenMM.
CMake will try to locate it automatically.

You will need a Lua library compatible with Lua API version >= 5.1 :
  * A release of Lua 5.x is usually already installed by default for recent linux versions, try to execute 'lua' and/or 'locate liblua'
  * You can Download and compile the official implementation : http://www.lua.org/download.html
  * The LuaJIT implementation can provide an important speedup : http://luajit.org/download.html
Please download and compile any of the two. See below for hints for setting DCMAKE_PREFIX_PATH in case of a manual install.
CMake will try to locate it automatically.

If a SQLite3 dynamic library is found on your system it will be used for linking; if not there is a provided one in 
the "./external" directory, it will be compiled and statically linked to the program.
See :
  * https://sqlite.org/index.html

The excellent Sol2 (header-only Lua<->c++ inteface) and LuaSQLite3 (a Lua/SQLite3 binding) dependancies are provided
in the './external' directory and are automatically included/compiled if required, and thus should not require extra
configuration.
See : 
  * https://github.com/ThePhD/sol2
  * http://lua.sqlite.org/index.cgi/home

----------------------------------------------
## COMPILE & INSTALL
----------------------------------------------

C and C++ compiler compatible with the C99 and C++14 standards are required.
Tested compilers:
  * gcc/g++ 5.4.0 and 6.2.0 (from Ubuntu 16.04 repositories)
  * clang/clang++ 3.8 (from Ubuntu 16.04 repositories)

Be sure to have CMake installed (http://www.cmake.org/), available on most repositories.
Tested version :
  * CMake 3.5.1 (from Ubuntu 16.04 repositories)

Create a build directory and move to that directory: 
  * mkdir build && cd build

For building a debug or release or an intermediate release with debug information, do:
  * cmake -DCMAKE_BUILD_TYPE=Debug ..
  * cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
  * cmake -DCMAKE_BUILD_TYPE=Release ..         (default, with cpu targeting for the current compilation machine)

Debug builds are possibly slower, usually more memory consuming, but useful when debugging with gdb or valgrind.

If OpenMM or your selected Lua implementation is installed to a non-standard (e.g. not /usr/local) CMake will probably not find them,
in that case add the path where the were installed to CMAKE_PREFIX_PATH, e.g.:
  * cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="$HOME/bin/openmm;$HOME/bin/luajit" ..

Then, once CMake built a Makefile without error, just execute the following in order to build the executable:
  * make

For a verbose make, use: 
  * make VERBOSE=1

For a faster multi-core (on N cores) build, use:
  * make -j N

Please never edit the autogenerated Makefile, edit the CMakeLists.txt instead.

For specifying another compiler on linux, for example clang, or a proprietary one like Intel icc (untested) : 
  * CC=clang CXX=clang++ cmake ..
  * CC=icc CXX=icpc cmake ..
 
In any case, you may need to edit the variables CMAKE_C_FLAGS_* and CMAKE_CXX_FLAGS_*
for setting proper levels of optimization. By default 

When using OpenMPI and forcing another C/C++ compilers version, the following might be required (possibly unsafe !): 
  * OMPI_CC=gcc-5 OMPI_CXX=g++-5 CC=gcc-5 CXX=g++-5 ccmake ..

----------------------------------------------
## LICENSING (all files excepted subdirectory external and its content)
----------------------------------------------
Copyright (c) 2016-2018, Florent Hédin, Tony Lelièvre, and École des Ponts - ParisTech
All rights reserved.

The 3-clause BSD license is applied to this software.

See LICENSE.txt

See in ./external/* for licensing information concerning the embedded libraries provided in './external' 

----------------------------------------------
## NOTES CONCERNING OpenMM
----------------------------------------------

The software should automatically detect the fastest OpenMM Platform available on your computer (i.e. CUDA, OpenCL, ...)
If not it will run with the slow Reference platform : it is most probably because the OpenMM directory with the
'plugins' library is not found.

You may need to export the following OPENMM_PLUGIN_DIR environment variable to solve the problem : 

For example for a custom installation in /home/$USER/bin/openmm
  * export OPENMM_PLUGIN_DIR=$HOME/bin/openmm/lib/plugins

For setting the amount of cpu threads used by each MPI process use the following env. variable :
  * export OPENMM_CPU_THREADS=1

