#! /usr/bin/env python2.7

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

psf = CharmmPsfFile('A.psf')
pdb = PDBFile('A.pdb')
params = CharmmParameterSet('top_all27_prot_lipid.rtf', 'par_all27_prot_lipid.prm')

# System Configuration
nonbondedMethod = CutoffNonPeriodic
nonbondedCutoff = 1.4*nanometers
constraints = HBonds
constraintTolerance = 0.00001

system = psf.createSystem(params,nonbondedMethod=nonbondedMethod,nonbondedCutoff=nonbondedCutoff,constraints=constraints)

# Integration Options
dt = 0.002*picoseconds
temperature = 300.00*kelvin
friction = 2/picosecond
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)

# do minimization, perform equilibration, then save the state of the simulation
equilibrationSteps = 25000
platform = Platform.getPlatformByName('CPU')
platformProperties = {'CpuThreads':'1'}

simulation = Simulation(psf.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(pdb.positions)

# Minimize and Equilibrate
print('Performing energy minimization...')
simulation.minimizeEnergy()
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

simulation.currentStep = 0
simulation.context.setTime(0.0)

state = simulation.context.getState(getPositions=True, getVelocities=True,
		                 getForces=True, getEnergy=True, getParameters=True)

PDBFile.writeFile(simulation.topology,state.getPositions(),open('after_equil_gas.pdb','w'))

#serialize the state to xml for load in parRep later
print('Serializing state : it contains coordinates, velocities, etc. for parRep...')
f = open('State.xml','w')
f.write(XmlSerializer.serialize(state))
f.close()

#serialize the system to xml for load in parRep later
print('Serializing system : it contains all definition of forces for parRep...')
f = open('System.xml','w')
f.write(XmlSerializer.serialize(system))
f.close()

#serialize the integrator to xml for load in parRep later
print('Serializing integrator for parRep...')
f = open('Integrator.xml','w')
f.write(XmlSerializer.serialize(integrator))
f.close()

