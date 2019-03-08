-- Input Lua script providing simulation parameters and functions for the parRep software
-- All lines starting with '--' are comments and thus ignored
-- Alanine dipeptide in gas phase, charmm FF

print("Lua version : ",_VERSION)

--------------------------------------------------------------------------------------------------------------
-- --------------- SIMULATION PARAMETERS ------------------------------
--------------------------------------------------------------------------------------------------------------

-- This is used for defining a maximum allowed running time for the simulation
-- The starting time off the program is saved at initialisation
-- in max_run_time_hours the user stores for example the maximum running time allowed by the queueing system when running on a cluster
--  fractional representation is allowed, e.g. 1.25 would be 1 hour and 15 minutes
max_run_time_hours = 24.0 -- default if not set in this script will be 24 hours
-- in minutes_to_stop_before_max_run_time the user requires the program to stop a few minutes before max_run_time_hours,
--  useful for large systems if the I/O may take some time; should be at least 1 minute, and no fractional value allowed
minutes_to_stop_before_max_run_time = 1 -- default if not set in this script will be 5 minutes

-- This uses OpenMM as MD engine
MD_Engine = "OpenMM"

-- OpenMM platform to use
--  AUTO : let OpenMM find the fastest platform (default)
--  REF  : on cpu, not optimised, not parallellised : slow !
--  CPU  : on cpu, optimised, parallellised with threads
--  OCL  : on cpu or gpu of any brand, or any accelerating device available
--  CUDA : on nvidia gpu only, usually the fastest
-- OMMplatform = "AUTO"
-- OMMplatform = "REF"
OMMplatform = "CPU"
-- OMMplatform = "OCL"
-- OMMplatform = "CUDA"

-- We can define here in an array extra platform specific properties passed to OpenMM
--  see http://docs.openmm.org/latest/userguide/library.html#platform-specific-properties
--  for a list of OpenMM platform specific properties 
-- Internally coerced to a std::map<std::string,std::string> so use a string for indexing (key, before the =)
--  and also be sure to define the values (after the =) also as strings (wrapped within ""), even if it is an integer

-- REF platform has no extra properties
-- ...

-- CPU platform properties : Threads = "1" is equivalent to defining OPENMM_CPU_THREADS=1 in the environnment
OMMplatform_properties = { Threads = "1"}

-- OpenCL platform properties
-- OMMplatform_properties = { 
--   Precision = "mixed", -- or "single" or "double"
--   UseCpuPme = "false", -- or true
--   OpenCLPlatformIndex = "0",
--   DeviceIndex = "0"
-- }

-- CUDA platform properties
-- OMMplatform_properties = {
--   Precision = "mixed", -- or "single" or "double"
--   UseCpuPme = "false", -- or "true"
--   DeviceIndex = "0",
-- --   CudaCompiler = "/path/to/nvcc",
--   UseBlockingSync = "false" -- or "true"
-- }

-- load the integrator parameters from a serialised OpenMM XML file ;
-- no default, error if undefined
--  adapt the path to the file in the following line
integrator = { xml = "./mol/ala2vac/Integrator.xml" }

-- load the OpenMM System from a serialised XML file
-- no default, error if undefined
system = { xml = "./mol/ala2vac/System.xml" }

-- load the OpenMM State from a serialised XML file
-- no default, error if undefined
state = { xml = "./mol/ala2vac/State.xml" }

-- parameters for energy minimisation : tolerance and maximum number of steps (if 0 : no limit, continue until tolerance satsfied)
-- defaults are : minimisation={Tolerance=1e-6,MaxSteps=0}
-- minimisation =
-- {
--   Tolerance = 0.01,
--   MaxSteps = 100
-- }

-- before the parRep algorithm actually starts, perform some steps of equilibration (MPI rank 0 performs dynamics alone)
-- default : 50e3
equilibrationSteps = 0

-- the total number of steps of the simulation ; together with the timestep read from the integrator file
--  it gives to total simulation time
-- no default, error if undefined
numSteps = 100e6

-- Define the type of ParRep simulation to perform, and provide its parameters
--  only one of the two following tables should be defined

-- FV ParRep (Lelièvre et al. algorithm): only "algorithm" and "GRobservables" are mandatory,
--  the others have default (c.f. documentation)
simulation =
{
  -- The parallel replica algorithm based on a Fleming-Viot process as defined by ...
  algorithm = "PARREP_FV",
  -- parameter for the ParRep Fleming-Viot stage
  --  A Fleming-Viot particle process is used for sampling the QSD,
  --  without any a priori defined decorrelation or dephasing time stage.
  --  Convergence is checked with Gelman-Rubin statistics
  -- a frequency (in steps) at which to verify if the system left the current state during FV procedure
  checkFV = 250, -- 0.5 ps
  -- a frequency (in steps) at which to evaluate convergence using the Gelman-Rubin statistics
  checkGR = 10, -- 20 fs
  -- check FV convergence if at least minAccumulatedObs have already been accumulated 
  --  This may be useful if there is a risk of early pseudo-convergence for some of the observables when only a few samples have been accumulated
  --  Therefore the minimum FV convergence time will be minAccumulatedObs*checkGR*dt (here 1000*10*dt = 20000 fs = 20 ps, but it is likely that the convergence time will be larger than that)
  --  Default is 100
  minAccumulatedObs = 1000,
  -- parameter for parallel dynamics stage
  --  checkDynamics is the frequency (in steps) at which to verify if the system left the current state
  checkDynamics = 2500, -- 5 ps
  -- Gelman-Rubin statistics checks convergence of several Observables: we define their name here, then there should be one function matching this name below in this file
  GRobservables = {"getEpot","getEkin","getPhi","getPsi"},
  -- Gelman-Rubin convergence : 0.01 means 1 %
  GRtol = 0.01
}

--------------------------------------------------------------------------------------------------------------
-- --------------- IMPLICIT VARIABLES AND FUNCTIONS ------------------------------
--------------------------------------------------------------------------------------------------------------

-- Together with the previous variables, the following implicit global variables and functions
-- are defined from the c++ code, and are accessible from the code you define below
--
-- ------------------
-- implicit variables (read only) :
-- ------------------
--
-- natoms : number of atoms of the system, read from OMM XML files
--
-- mpi_rank_id   : the id of the MPI process working on this file
-- mpi_num_ranks : the total number of MPI processes running
--
-- epot,ekin,etot : 3 variables containing the values of the potential, kinetic and total energies (kcal/mol),
--  of the system for the current simulation time
--
-- timeStep : the MD timeStep, read from OMM XML files
--
-- referenceTime : the reference clock time, simulation will stop when referenceTime >= (numSteps*timestep),
--                 where timestep is defined within the xml OMM integrator file 
--
-- temperature : T of the system in Kelvin, initial value read from OMM XML files ; use get_temperature() for instantaneous values
--
-- ------------------
-- implicit functions : 
-- ------------------
--
-- exit_from_lua() : call this if it is required to finish the simulation from the lua script : it will terminate MPI properly
--  and flush I/O files ; but it won't perform a DB backup so do it manually before.   
--
-- get_coordinates(n) : function returning a tuple x,y,z containing the coordinates of atom n (in nm)
--  NOTE n is indexed in Lua style starting from 1, internally it will be accessed as (n-1) i.e. C++ style
--
-- get_velocities(n) : function returning a tuple vx,vy,vz containing the velocities of atom n (in nm/picosecond)
--  NOTE n is indexed in Lua style starting from 1, internally it will be accessed as (n-1) i.e. C++ style
--
-- get_all_coordinates() : function returning a table crds containing a safe read-only copy of the coordinates
--  access by using : crds.x[i] crds.y[i] crds.z[i] for respectively getting x,y,z coordinates of atom i 
--  NOTE lua indexing, starting from 1
--
-- get_all_velocities() : same as gel_all_coordinates but for velocities ; returns a table vels such that the access is : vels.x[i], etc.
--
-- get_all_crdvels() : returns a 2-tuple crds,vels containing coordinates and vels : internally it call the 2 above defined functions
--
-- get_pbc() : returns periodic boundary conditions as a table pbc with members pbc.a pbc.b pbc.c, each being a xyz vector (i.e. pbc.a.x pbc.a.y pbc.a.z) 
--
-- set_pbc(pbc) : set the openmm periodic boundary conditions to be used, the atgument is a table pbc a described above
--
-- NOTE possibly slow (especially if a copy to a OCL or CUDA device is required), use rarely for performance
-- set_all_coordinates(crds)  : uses a crds tables (as described above) and use it for setting c++/openMM coordinates
-- set_all_velocities(vels)   : the same with velocities
-- set_all_crdvels(crds,vels) : does with one function only what the 2 above defined do
--
-- get_mass(n) : function returning the mass of atom n in a.m.u.
--  NOTE n is indexed in Lua style starting from 1, internally it will be accessed as (n-1) i.e. C++ style
--
-- get_temperature() : get instantaneous value of the temperature in K
--
-- get_COM() : function returning the center of mass of the system as a tuple x,y,z
--
-- get_COM_idxs(idxs) : function returning the center of mass of a subset of the system as a tuple x,y,z
--  NOTE this time idxs is indexed directly in C++ style
--  for example get_COM_idxs({1,2,3}) to get COM of atoms 1, 2 and 3  (C++ : atoms 0, 1, 2)
--
-- get_minimised_energy(tolerance,maxSteps) : this function returns the minimised energy of the system, using the OpenMM L-BFGS minimiser
--  note that coordinates are not affected, it just returns the minimum epot of the bassin in which dynamics currently evolves
--  it returns a 3-tuple ep,ek,et (potential, kinetic and total energy)
--  the tolerance and maxSteps can be the above defined minimisation.Tolerance and minimisation.MaxSteps
--
-- get_minimised_crdvels(tolerance,maxSteps) : this function returns a 2-tuple (crds,vels) containing
--  a copy of coordinates and velocities after minimisation.
--  crds and vels are both tables with x,y,z members, each of size natoms,  : e.g. crds.x[i] returns the x coordinate of atom i, idem for vels.x[i]
--  note that C++/OpenMM coordinates are not modified if modifying this table : this is a safe read-only copy
--  the tolerance and maxSteps can be the above defined simulation.minimisationTolerance and simulation.minimisationMaxSteps
--  NOTE lua indexing, starting from 1
--
-- get_minimised_energy_crdvels(tolerance,maxSteps) : this returns a 5-tuple (ep,ek,et,crds,vels) 
--  This does the same as the two previous functions but with only one call
--
-- hr_timer() : returns a variable representing a c++ high precision timer : can be used for measuring execution time.
--  do not try to modify it or even read it, it should only be used as argument for the following hr_timediff_* functions.
--
-- hr_timediff_ns(before,after) : returns the time difference in nanoseconds between two hr_timer() 'before' and 'after' : usage:
--
--      local bf = hr_timer()
--      function_to_profile()
--      local af = hr_timer()
--      print('Exec. time of function function_to_profile() is (ns) : ',hr_timediff_ns(bf,af))
--
-- hr_timediff_us() and hr_timediff_ms() : same as above but exec time is returned respectively in microseconds and milliseconds
--

--------------------------------------------------------------------------------------------------------------
-- --------------- USER DEFINED VARIABLES AND FUNCTIONS ------------------------------
--------------------------------------------------------------------------------------------------------------

-- Some of the following VARIABLES and FUNCTIONS are mandatory and called from C++ (if it is the case it is explicitly documented)
-- If not they can be restrited to this file using the local keyword

-- Define here local variables and functions used later within state_init() and check_state_left()

--------------------------------------------------------------------------------------------------------------
-- --------------- FUNCTIONS DEFINING A PARREP STATE ------------------------------
--------------------------------------------------------------------------------------------------------------

-- TWO functions, state_init() and check_state_left(), will be called from c++ code to know if the 
--  dynamics left the current state. You are free to define the state in any way, using variables defined explicitly in this file
--  or implicitly (c++ interface, see above).

-- Define here local variables and functions used later within state_init() and check_state_left()

-- atom index definition of each angle (starting at 1, lua style)
local phi_def = {5,7,9,15}
local psi_def = {7,9,15,17}

local fromState,toState = 'unknown','unknown'

-- calculates a dihedral angle between 2 planes
--  idx contains indices of 4 atoms used for defining the 2 planes
--  if crds==nil then cordinates retrieved using get_coordinates, otherwise they are read from this table crds
local function calcDihe(idx,crds)
  
  -- multiplies a vector by a scalar : c[.] = vec[.] * scalar
  local function mulScalVec(scalar,vec)
    local c={0.0,0.0,0.0}
    c[1] = scalar*vec[1]
    c[2] = scalar*vec[2]
    c[3] = scalar*vec[3]
    return c
  end
  
  -- returns dot product ; expects 2 vectors of length 3
  local function dotProduct(a,b)
    return (a[1]*b[1] + a[2]*b[2] + a[3]*b[3])
  end
  
  -- returns the vector corresponding to the cross product of u and v ; length 3
  local function crossProduct(a,b)
    local c={0.0,0.0,0.0}
    c[1] = a[2]*b[3] - a[3]*b[2]
    c[2] = a[3]*b[1] - a[1]*b[3]
    c[3] = a[1]*b[2] - a[2]*b[1]
    return c
  end
  
  -- returns norm of vector
  local function vecNorm(v)
    return math.sqrt(dotProduct(v,v))
  end
  
  -- see wikipedia : https://en.wikipedia.org/wiki/Dihedral_angle#Calculation_of_a_dihedral_angle
  -- Any plane can  be described by two non-collinear vectors lying in that plane;
  -- taking their cross product yields a normal vector to the plane.
  -- Thus, a dihedral angle can be defined by three vectors, b1, b2 and b3,
  -- forming two pairs of non-collinear vectors.
  
  local x1,y1,z1 = 0.,0.,0.
  local x2,y2,z2 = 0.,0.,0.
  if(crds==nil) then
    x1,y1,z1 = get_coordinates(idx[1])
    x2,y2,z2 = get_coordinates(idx[2])
  else
    x1,y1,z1 = crds.x[idx[1]],crds.y[idx[1]],crds.z[idx[1]]
    x2,y2,z2 = crds.x[idx[2]],crds.y[idx[2]],crds.z[idx[2]]
  end
  
  local b1 = {x2-x1,y2-y1,z2-z1}
  b1 = mulScalVec(-1.0,b1)
  
  if(crds==nil) then
    x1,y1,z1 = get_coordinates(idx[2])
    x2,y2,z2 = get_coordinates(idx[3])
  else
    x1,y1,z1 = crds.x[idx[2]],crds.y[idx[2]],crds.z[idx[2]]
    x2,y2,z2 = crds.x[idx[3]],crds.y[idx[3]],crds.z[idx[3]]
  end
  
  local b2 = {x2-x1,y2-y1,z2-z1}
  
  if(crds==nil) then
    x1,y1,z1 = get_coordinates(idx[3])
    x2,y2,z2 = get_coordinates(idx[4])
  else
    x1,y1,z1 = crds.x[idx[3]],crds.y[idx[3]],crds.z[idx[3]]
    x2,y2,z2 = crds.x[idx[4]],crds.y[idx[4]],crds.z[idx[4]]
  end

  local b3 = {x2-x1,y2-y1,z2-z1}
  
  -- cross-product between b1 and b2
  local cp12 = crossProduct(b1,b2)
  
  -- and between b3 and b2
  local cp32 = crossProduct(b3,b2)
  
  -- cp between the 2 normal vectors
  local cpcp = crossProduct(cp12,cp32)
  
  local y = dotProduct(cpcp, b2)*(1.0/vecNorm(b2))
  local x = dotProduct(cp12, cp32)
  local dihe = math.atan2(y,x) * 180.0/math.pi
  
  return dihe
  
end

-- this function is mandatory and called from C++, program will fail if not defined
--  it should take no arguments
--  it returns nothing
-- Use it if you have global variables used in check_state_left() (or other functions) that you need to initialise
-- It will be called only once after equilibration from C++, before starting the parRep or parRep_FV algorithm
-- It can also be called at any time from this file if required
function state_init()

  -- value of the phi and psi dihedral angles at initialisation
  local phi = calcDihe(phi_def,nil)
  local psi = calcDihe(psi_def,nil)
  
  local domain = nil
  if( (phi > 0.0 and phi < 120.0) and (psi < 0.0 and psi > -150.0) ) then
    domain = 'C_ax'
  else
    domain = 'C_eq'
  end

  fromState = domain
  toState   = domain
  
  print("Initial state is: "..domain.." {"..phi.." "..psi.."} ")
  
end

-- this function is mandatory and called from C++, program will fail if not defined
--  it should take no arguments
--  it should return a boolean : true in case the dynamics left the state, false otherwise
-- You may create as many functions as you want and call them from check_state_left(),
--  but the c++ code will in the end only call check_state_left()
function check_state_left()

  -- value of the phi and psi dihedral angles of ala2 at a given time value
  local phi = calcDihe(phi_def,nil)
  local psi = calcDihe(psi_def,nil)

  local escaped = false
  
  local domain = nil
  -- first check if within a rectangular domain around C_ax
  if( (phi > 0.0 and phi < 120.0) and (psi < 0.0 and psi > -150.0) ) then
    domain = 'C_ax'
  else
    domain = 'C_eq'
  end
  
  toState = domain
  if(fromState == toState) then escaped = false else escaped = true end
  
  return escaped

end

-- this function is mandatory and called from C++, program will fail if not defined
--  it should take no arguments
--  it should return a boolean : true in case the system is currently outside of any of the known states, false otherwise
-- If the configuration space is a partition the system enters a new state as soon as it exits another,
--  and therefore this function can just return false without doing anything
function check_transient_propagation_required()
  -- we have partitioned the (phi,psi) configuration space in two states so this returns false
  return false
end

--------------------------------------------------------------------------------------------------------------
-- --------------- DATABASE OF STATES STORAGE FUNCTIONS ------------------------------
--------------------------------------------------------------------------------------------------------------
--
-- functions and variables use for storing data to a database are defined as members of the table SQLiteDB
--  the following can be defined as they are used from the c++ code : 
--    SQLiteDB.insert_statement_states
--    SQLiteDB.insert_statement_crdvels
--    function SQLiteDB.open()
--    function SQLiteDB.close()
--    function SQLiteDB.insert_state()
--    function SQLiteDB.backup_to_file()
--
-- if not defined within this file there are default empty functions defined, doing nothing

-- only rank 0 manages a database
if(mpi_rank_id==0)
then

  -- parameters for the SQLite3 database used for storing information about parRep states
  -- The database is stored in memory for performance, but a backup is regularly performed to file 'name'
  -- defaults are : database={name="run.db",backupFrequency=500.0}
  database = 
  {
    -- Because the structure of the algorithms, each MPI rank will have its own database and each will be saved to an altered
    --  fileName where the MPI rank id is added; for instance with name='run.db'
    --  files 'run.0.db', 'run.1.db', ... would be generated i.e. id is inserted before extension dot
    name = "ala2vac.db",
    -- frequency, in ps, at which to backup the database to a file
    backupFrequency = 10 -- 10 ps
  }
  
  SQLiteDB.insert_statement_states  = [[ INSERT INTO STATES (PARREP_DONE,REF_TIME,ESC_TIME,STATE_FROM,STATE_TO,TAU)
                                                    VALUES ($lprep,$reft,$esct,$state_from,$state_to,$tau); ]]

  SQLiteDB.insert_statement_crdvels = [[ INSERT INTO CRDVELS (ID,X,Y,Z,VX,VY,VZ)
                                                    VALUES ($id,$x,$y,$z,$vx,$vy,$vz); ]]

  function SQLiteDB.open()

    print('Lua SQLite3 opening an on-disk db named ',database.name)
    
    -- open in memory database and enable foreign keys
    SQLiteDB.db = sqlite3.open(database.name)
    SQLiteDB.db:exec("PRAGMA foreign_keys = ON;")
    
    SQLiteDB.db:exec("BEGIN TRANSACTION;")

    -- create the states table : it contains the escape time for this state
    SQLiteDB.db:exec[[ CREATE TABLE STATES(
                        ID          INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL,
                        PARREP_DONE INTEGER NOT NULL,
                        REF_TIME    REAL,
                        ESC_TIME    REAL,
                        STATE_FROM  TEXT,
                        STATE_TO    TEXT,
                        TAU         REAL); ]]

    -- create the crdvels table : it contains coordinates and velocities of the system for a given 'states' record
    SQLiteDB.db:exec[[ CREATE TABLE CRDVELS(
                        ID INTEGER NOT NULL,
                        X   REAL,
                        Y   REAL,
                        Z   REAL,
                        VX  REAL,
                        VY  REAL,
                        VZ  REAL,
                        FOREIGN KEY(ID) REFERENCES STATES(ID) ); ]]

    SQLiteDB.db:exec("END TRANSACTION;")

  end

  function SQLiteDB.close()

    print('Lua SQLite3 closing db on rank ',mpi_rank_id)
    SQLiteDB.db:close()

  end

  local stateID=0

  -- the first three arguments are always provided from c++ code
  -- extra arguments might be provided by the ... and should be retrieved from Lua's side using the ... special token, converted to a table (see args below)
  function SQLiteDB.insert_state(parRepDone,tauTime,escapeTime)
    
    local ref_time=referenceTime
    local from=fromState
    local to=toState

    SQLiteDB.db:exec("BEGIN TRANSACTION;")
    
    local stmt = SQLiteDB.db:prepare(SQLiteDB.insert_statement_states)
    stmt:bind_names{lprep=parRepDone, reft=ref_time, esct=escapeTime, state_from=from, state_to=to, tau=tauTime}
    stmt:step()
    stmt:finalize()
    
    stateID = stateID+1
    
    for n=1,natoms
    do
      local x,y,z = get_coordinates(n)
      local vx,vy,vz = get_velocities(n)
      
      stmt = SQLiteDB.db:prepare(SQLiteDB.insert_statement_crdvels)
      stmt:bind_names{ id=stateID, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz }
      stmt:step()
      stmt:finalize()
      
    end
    
    SQLiteDB.db:exec("END TRANSACTION;")

  end

end -- if(mpi_rank_id==0)
  
--------------------------------------------------------------------------------------------------------------
-- --------------- GELMAN RUBIN FUNCTIONS ESTIMATING OBSERVABLES ------------------------------
--------------------------------------------------------------------------------------------------------------

-- Define a function for calculating the value of each Observable
-- Those functions should:
-- 1) take no arguments 
-- 2) return a double precision (any Lua numeric value is returned as a double precision)
--
-- Those GR functions are called from C++ if they were listed in simulation.GRobservables
--

-- Definition of the "getEpot" and "getEkin" observables used for Gelman-Rubin statistics
-- Just returns value of the potential and kinetic energy
function getEpot()
  return epot
end

function getEkin()
  return ekin
end

function getEtot()
  return etot
end

-- Definition of the "getPhi" and "getPsi" Observable used for Gelman-Rubin statistics
-- Returns current value of the phi and psi dihedral angles
function getPhi()
  return calcDihe(phi_def,nil)
end

function getPsi()
  return calcDihe(psi_def,nil)
end

--------------------------------------------------------------------------------------------------------------
    
