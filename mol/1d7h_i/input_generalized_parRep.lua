-- All lines starting with '--' are comments and thus ignored
--
-- Input Lua script providing simulation parameters and functions for the parRep software
-- TEST system --> FKBP protein with DMSO ligand (PDB : 1d7h), with implicit solvent model, using Amber14+GAFF+OBC2
-- GOAL --> estimating k_off (vound to unbound) and k_on (unbound to bound) using parRep

print("Lua version : ",_VERSION)

-- show global variables
--  all tables/functions/variables marked as having a default value below should be listed there
-- for key,value in pairs(_G) do
--   print("found Lua global -> " .. key)
-- end

--------------------------------------------------------------------------------------------------------------
-- --------------- SIMULATION PARAMETERS ------------------------------
--------------------------------------------------------------------------------------------------------------

-- This is used for defining a maximum allowed running time for the simulation
-- The starting time off the program is saved at initialisation
-- in max_run_time_hours the user stores for example the maximum running time allowed by the queueing system when running on a cluster
--  fractional representation is allowed, e.g. 1.25 would be 1 hour and 15 minutes
max_run_time_hours = 6.0 -- default if not set in this script will be 24 hours
-- in minutes_to_stop_before_max_run_time the user requires the program to stop a few minutes before max_run_time_hours,
--  useful for large systems if the I/O may take some time; should be at least 1 minute, and no fractional value allowed
minutes_to_stop_before_max_run_time = 2 -- default if not set in this script will be 5 minutes

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
integrator = { xml = "./mol/1d7h_i/Integrator.xml" }

-- load the OpenMM System from a serialised XML file
-- no default, error if undefined
system = { xml = "./mol/1d7h_i/System.xml" }

-- load the OpenMM State from a serialised XML file
-- no default, error if undefined
state = { xml = "./mol/1d7h_i/State.xml" }

-- parameters for energy minimisation : tolerance and maximum number of steps (if 0 : no limit, continue until tolerance satsfied)
-- defaults are : minimisation={Tolerance=1e-6,MaxSteps=0}
-- minimisation =
-- {
--   Tolerance = 0.01,
--   MaxSteps = 100
-- }

-- before the parRep algorithm actually starts, perform some steps of equilibration (MPI rank 0 performs dynamics alone)
-- default : 50e3
-- can be small or 0 if system provided in XML files was already equilibrated
equilibrationSteps = 0

-- the total number of steps of the simulation ;
-- together with the timestep read from the integrator file, it gives an idea of the total simulation time before stopping
-- no default, error if undefined
numSteps = 50e6

-- Define the type of ParRep simulation to perform, and provide its parameters

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
  checkFV = 2000, -- 2 ps
  -- a frequency (in steps) at which to evaluate convergence using the Gelman-Rubin statistics
  checkGR = 100,  -- 100 fs
  -- parameter for parallel dynamics stage
  --  checkDynamics is the frequency (in steps) at which to verify if the system left the current state
  checkDynamics = 2000, -- 2 ps
  -- Gelman-Rubin statistics checks convergence of one or several Observables: we define their name here, then there should be one function matching this name below in this file
  GRobservables = {"getCOMsDist","rmsVelDMSO","get_dist_1","get_dist_2"},
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
--  NOTE this time idxs is indexed in C++ style
--  for example get_COM_idxs({0,1,2}) to get COM of atoms 1, 2 and 3  (C++ : atoms 0, 1, 2)
--
-- get_minimised_energy(tolerance,maxSteps) : this function returns the minimised energy of the system, using the OpenMM L-BFGS minimiser
--  note that coordinates are not affected, it just returns the minimum epot of the bassin in which dynamics currently evolves
--  it returns a 3-tuple ep,ek,et (potential, kinetic and total energy)
--  the tolerance and maxSteps can be the above defined simulation.minimisationTolerance and simulation.minimisationMaxSteps
--
-- get_minimised_crdvels(tolerance,maxSteps) : this function returns a 2-tuple (crds,vels) containing
--  a copy of coordinates and velocities after minimisation.
--  crds and vels are both tables with x,y,z members, each of size natoms,  : e.g. crds.x[i] returns the x coordinate of atom i, idem for vels.x[i]
--  note that C++/OpenMM coordinates are not modified if modifying this table : this is a safe read-only copy
--  the tolerance and maxSteps can be the above defined simulation.minimisationTolerance and simulation.minimisationMaxSteps
--  NOTE lua indexing, starting from 1
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
-- If not they can be restricted to this file using the local keyword

-- Define here local variables and functions used later within state_init() and check_state_left()

-- from Caflisch 2011 article
-- the following are DMSO not-H atoms
local dmso = {1663,1664,1665,1666}

-- From Caflisch 2016 article
-- the first distance we track to know if there is an unbinding is the distance between the amide hydrogen of residue 56 (c++ idx 869) and the oxygen of DMSO (c++ idx 1664)
local index_dist_1 = {869,1664}
-- the second distance is between the COM of the ring of residue 59 (c++ idxs 930,931,933,935,937,939) and the sulfur of the DMSO (c++ idx 1663)
local index_dist_2 = {930,931,933,935,937,939,1663}

-- there are 2 possible states for DMSO : bound ('b') or unbound ('u')
--  intial state will be determined in state_init()
local currState = nil

-- definition of binding <-> unbinding thresholds, following article from
--  A. Caflisch, JCTC 2016
-- a bound state becomes unbound if the distance is GREATER THAN (angstroems):
local b_to_u_d1 = 12.0 
local b_to_u_d2 = 12.0 

-- for the moment we only want to estimate b -> u transitions so after observing a state u
--  we will go back to b : for this we need a backup of coordinates, velocities and PBCs on Lua's side
local b_backup = { crds={}, vels={}, pbc={} }

-- useful for selecting a slice of an array
function table.slice(tbl, first, last, step)
  local sliced = {}

  for i = first or 1, last or #tbl, step or 1 do
    sliced[#sliced+1] = tbl[i]
  end

  return sliced
end

--------------------------------------------------------------------------------------------------------------
-- --------------- FUNCTIONS DEFINING A PARREP STATE ------------------------------
--------------------------------------------------------------------------------------------------------------

-- TWO functions, state_init() and check_state_left(), will be called from c++ code to know if the 
--  dynamics left the current state. You are free to define the state in any way, using variables defined explicitly in this file
--  or implicitly (c++ interface, see above).


-- this function is mandatory and called from C++, program will fail if not defined
--  it should take no arguments
--  it returns nothing
-- Use it if you have global variables used in check_state_left() (or other functions) that you need to initialise
-- It will be called only once after equilibration from C++, before starting the parRep or parRep_FV algorithm
-- It can also be called at any time from this file if required
function state_init()

  -- access coordinates of the 2 atoms defining group index_dist_1
  local d1_h_x,d1_h_y,d1_h_z = get_coordinates(index_dist_1[1]+1)
  local d1_o_x,d1_o_y,d1_o_z = get_coordinates(index_dist_1[2]+1)

  -- access the com of the 6-ring and the S atom defining group index_dist_2
  local d2_ring_x,d2_ring_y,d2_ring_z = get_COM_idxs(table.slice(index_dist_2,1,6))
  local d2_s_x,d2_s_y,d2_s_z = get_coordinates(index_dist_2[7]+1)

  local d1 = 10.0*math.sqrt( (d1_h_x-d1_o_x)^2    + (d1_h_y-d1_o_y)^2    + (d1_h_z-d1_o_z)^2 )
  local d2 = 10.0*math.sqrt( (d2_ring_x-d2_s_x)^2 + (d2_ring_y-d2_s_y)^2 + (d2_ring_z-d2_s_z)^2 )
  
  print('Initial distances d1,d2 are : ',d1,' and ',d2)
  
  -- we require BOTH d1 and d2 to be strictly larger than b_to_u for having an unbound state
  if(d1>b_to_u_d1 and d2>b_to_u_d2) then
    currState = 'u'
    print('Initial DMSO state appears to be : UNBOUND')
  else
    currState = 'b'
    print('Initial DMSO state appears to be : BOUND')
  end
  
  -- backup
  if currState=='b' then
    print('Performing a Lua backup of bound pbc, coordinates and velocities')
    b_backup.pbc = get_pbc()
    b_backup.crds,b_backup.vels = get_all_crdvels()
  else
    print('Found an unbound state ; switching back to a bound one for next parrep iteration')
    set_pbc(b_backup.pbc)
    set_all_crdvels(b_backup.crds,b_backup.vels)
    currState = 'b'
  end

end

-- You may create as many functions as you want and call them from check_state_left(),
--  but the c++ code will in the end only call check_state_left()

-- this function is mandatory and called from C++, program will fail if not defined
--  it should take no arguments
--  it should return a boolean : true in case the dynamics left the state, false otherwise
function check_state_left()

  -- access coordinates of the 2 atoms defining group index_dist_1
  local d1_h_x,d1_h_y,d1_h_z = get_coordinates(index_dist_1[1]+1)
  local d1_o_x,d1_o_y,d1_o_z = get_coordinates(index_dist_1[2]+1)

  -- access the com of the 6-ring and the S atom defining group index_dist_2
  local d2_ring_x,d2_ring_y,d2_ring_z = get_COM_idxs(table.slice(index_dist_2,1,6))
  local d2_s_x,d2_s_y,d2_s_z = get_coordinates(index_dist_2[7]+1)

  local d1 = 10.0*math.sqrt( (d1_h_x-d1_o_x)^2    + (d1_h_y-d1_o_y)^2    + (d1_h_z-d1_o_z)^2 )
  local d2 = 10.0*math.sqrt( (d2_ring_x-d2_s_x)^2 + (d2_ring_y-d2_s_y)^2 + (d2_ring_z-d2_s_z)^2 )
  
  local newState = currState

  -- we require BOTH d1 and d2 to be strictly larger than b_to_u for having an unbound state
  if(d1>b_to_u_d1 and d2>b_to_u_d2) then
    currState = 'u'
  else
    currState = 'b'
  end
  
  print('Distances d1,d2 are : ',d1,' and ',d2)

  print('DMSO state appears to be :',currState)
  
  if (newState == currState) then
    return false
  else
    return true
  end

end

-- this function is mandatory and called from C++, program will fail if not defined
--  it should take no arguments
--  it should return a boolean : true in case the system is currently outside of any of the known states, false otherwise
-- If the configuration space is a partition the system enters a new state as soon as it exits another, and therefore this function can just return false
--  without doing anything
function check_transient_propagation_required()
  -- we have partitioned the configuration space in two states so this returns false
  return false
end

--------------------------------------------------------------------------------------------------------------
-- --------------- DATABASE OF STATES STORAGE FUNCTIONS ------------------------------
--------------------------------------------------------------------------------------------------------------

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
    name = "1d7h.db",
    -- frequency, in ps, at which to backup the database to a file
    backupFrequency = 100.0
  }

  SQLiteDB.insert_statement_states  = [[ INSERT INTO STATES (PARREP_DONE,REF_TIME,ESC_TIME,STATE_FROM,STATE_TO,TAU,EXIT_ORDER)
                                                    VALUES ($lprep,$reft,$esct,$state_from,$state_to,$tau,$exit_order); ]]

  SQLiteDB.insert_statement_crdvels = [[ INSERT INTO CRDVELS (ID,X,Y,Z,VX,VY,VZ)
                                                    VALUES ($id,$x,$y,$z,$vx,$vy,$vz); ]]

  function SQLiteDB.open()

    print('Lua SQLite3 opening an on disk db named ',database.name,' on rank ',mpi_rank_id,' of ',mpi_num_ranks)
    
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
                        TAU         REAL,
                        EXIT_ORDER  INTEGER); ]]

    -- create the crdvels table : it contains coordinates and velocities of the system for a given 'states' record
    SQLiteDB.db:exec[[ CREATE TABLE CRDVELS(
                        ID INTEGER NOT NULL,
                        X   REAL,
                        Y   REAL,
                        Z   REAL,
                        VX  REAL,
                        VY  REAL,
                        VZ  REAL,
                        FOREIGN KEY(ID) REFERENCES STATES(ID)
                                          ); ]]

    SQLiteDB.db:exec("END TRANSACTION;")

  end

  function SQLiteDB.close()

    print('Lua SQLite3 closing in-memory db on rank ',mpi_rank_id,' of ',mpi_num_ranks)
    SQLiteDB.db:close()

  end

  local stateID=0

  -- the first three arguments are always provided from c++ code
  -- extra arguments might be provided by the ... and should be retrieved from Lua's side using the hidden metatable arg
  function SQLiteDB.insert_state(parRepDone,tauTime,escapeTime,...)
    
    local ref_time = referenceTime
    local from,to  = nil,nil
    
    if currState == 'u' then
      from = 'b'
      to   = 'u'
    else
      from = 'u'
      to   = 'b'
    end
    
    -- default exit_order if not provided by c++ is 1
    local order=1
    
    -- extra arguments stored in ... always go as a pair key,value the key always being a string 
    local args = {...}
    if #args>0 then
      
      for i=1,#args,2
      do
        
        local key,value = args[i],args[i+1]
        
        if key=="exit_order" then
          order=value
  --       elseif key=="check_states" then
  --         state_init()
  --         -- ...
        elseif key=="ignore_ref_time" then
          -- in some cases the exit time might be inserted before the exact reference_time is known: in that case the variable is set to NaN
          ref_time = 0/0 -- depending on lua's version it should produce nan or -nan
        end
        
      end
      
    end
    
    SQLiteDB.db:exec("BEGIN TRANSACTION;")
    
    local stmt = SQLiteDB.db:prepare(SQLiteDB.insert_statement_states)
    stmt:bind_names{lprep=parRepDone, reft=ref_time, esct=escapeTime, state_from=from, state_to=to, tau=tauTime, exit_order=order}
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

  end -- SQLiteDB.insert_state

end -- (mpi_rank_id==0)

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

-- Returns instantaneous value of the distance between the COMs defined above, i.e. the criterion defining the state
function getCOMsDist()

  local d_cx,d_cy,d_cz = get_COM_idxs(dmso)
  local cx,cy,cz = get_COM()

  local dist = math.sqrt( (d_cx-cx)^2 + (d_cy-cy)^2 + (d_cz-cz)^2 )

  return dist

end

-- returns the root mean square velocity of DMSO molecule
function rmsVelDMSO()

  local vel = {}

  for k,v in ipairs(dmso)
  do
    -- get_velocities(n) expects n as Lua index but in DMSO they are c++ style (starting at 0), so we add 1
    local vx,vy,vz = get_velocities(v+1)
    vel[k] = vx*vx + vy*vy + vz*vz
  end

  -- then calculate and return the RMS velocity
  local rms = 0
  for i=1,#vel
  do
    rms = rms + vel[i]
  end

  return (1/#vel)*math.sqrt(rms)

end

function get_dist_1()

  -- access coordinates of the 2 atoms defining group index_dist_1
  local d1_h_x,d1_h_y,d1_h_z = get_coordinates(index_dist_1[1]+1)
  local d1_o_x,d1_o_y,d1_o_z = get_coordinates(index_dist_1[2]+1)

  local d1 = math.sqrt( (d1_h_x-d1_o_x)^2 + (d1_h_y-d1_o_y)^2 + (d1_h_z-d1_o_z)^2 )
  
  return d1
  
end

function get_dist_2()

  -- access the com of the 6-ring and the S atom defining group index_dist_2
  local d2_ring_x,d2_ring_y,d2_ring_z = get_COM_idxs(table.slice(index_dist_2,1,6))
  local d2_s_x,d2_s_y,d2_s_z = get_coordinates(index_dist_2[7]+1)

  local d2 = math.sqrt( (d2_ring_x-d2_s_x)^2 + (d2_ring_y-d2_s_y)^2 + (d2_ring_z-d2_s_z)^2 )
  
  return d2
  
end

--------------------------------------------------------------------------------------------------------------

