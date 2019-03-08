# v1.2.0

Mar 08th, 2019

The software was modified in order to allow reproducibility in some cases (the limiting factor is OpenMM which does not
always provides deterministic output even when using the same seeds (!), see http://docs.openmm.org/latest/userguide/library.html#determinism ).

The main executable now has 2 more command line options, '--inp-seeds [fname]' or '--out-seeds [fname]' for respectively loading seeds or writting seeds 
from/to a unique binary file. See rand.hpp/rand.cpp, od the Doxygen doc for more details.

These modifications now allow Continuous Integration (CI) on the infrastructure provided by INRIA : in mol/ci a small test case will be executed at each commit to the repository and compared to reference results.

The two other minor modifications concern the Lua scripts:

* "get_minimised_energy_crdvels" was added to the set of functions that the user can call from the Lua script : it simply combines in one call what 
"get_minimised_energy" and "get_minimised_crdvels" already provided.

* a extra optional parameter is available for the "simulation" parameters block when "simulation.algorithm" is "PARREP_FV" ; this parameter is
"simulation.minAccumulatedObs" : it will enforce that at least minAccumulatedObs observations of an observable have already been accumulated
before the convergence test is performed ; this may be useful if there is a risk of early pseudo-convergence for some of the observables when only a few samples have been accumulated.


**Download sources:**

https://gitlab.inria.fr/parallel-replica/gen.parRep/tags/v1.2.0

or

https://github.com/FHedin/gen.parRep/releases/tag/v1.2.0

----------------------------------------------

# v1.1.0

Nov 16th, 2018

Finished the implementation of a "MD_interface" abstract class for future addition of different Molecular Dynamics engines.
This required important changes to a significant part of the code.
Optimizations to the MPI communications were also performed, therefore the minor version number
of the software was incremented.

**Download sources:**

https://gitlab.inria.fr/parallel-replica/gen.parRep/tags/v1.1.0

or

https://github.com/FHedin/gen.parRep/releases/tag/v1.1.0

----------------------------------------------

# v1.0.1

Aug 16th, 2018

Simply added directory structure to readme file.

**Download sources:**

https://gitlab.inria.fr/parallel-replica/gen.parRep/tags/v1.0.1

or

https://github.com/FHedin/gen.parRep/releases/tag/v1.0.1

----------------------------------------------

# v1.0.0

Aug 13th, 2018

Release of v1.0.0 : several bug fixes like memory leaks addressed, added input files for alanine
dipeptide using transient propagation.

**Download sources:**

https://gitlab.inria.fr/parallel-replica/gen.parRep/tags/v1.0.0

or

https://github.com/FHedin/gen.parRep/releases/tag/v1.0.0

----------------------------------------------

# 1.0.0-rc1

July 6th, 2018

First public release of the software.
It is following the submission of the corresponding article on arXiv ; this is a release candidate 1.0.0-rc1,
the 1.0.0 release should follow within the next weeks.

**Download sources:**

https://gitlab.inria.fr/parallel-replica/gen.parRep/tags/1.0.0-rc1

or

https://github.com/FHedin/gen.parRep/releases/tag/1.0.0-rc1
