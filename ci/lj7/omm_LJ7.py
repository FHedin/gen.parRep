from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import numpy as np
from math import *
from sys import stdout

platform = Platform.getPlatformByName('Reference')

# periodic conditions
# periodic vectors
a = Vec3(20.0,0.0,0.0)
b = Vec3(0.0,20.0,0.0)
c = Vec3(0.0,0.0,20.0)

center = Vec3(10.,10.,10.)

print("a = ",a)
print("b = ",b)
print("c = ",c)
print("center = ",center)

natoms = 7
mass = 39.948*amu
temperature = 15.0*kelvin

integrator = LangevinIntegrator(temperature, 10/picosecond, 0.001*picosecond)
integrator.setConstraintTolerance(1e-8)

integrator.setRandomNumberSeed(123456789)

nb = NonbondedForce()
nb.setNonbondedMethod(NonbondedForce.NoCutoff)
#nb.setUseSwitchingFunction(True)
#nb.setSwitchingDistance(1.0*nanometer)
#nb.setCutoffDistance(1.2*nanometer)
nb.setUseDispersionCorrection(False)

# see http://www.sklogwiki.org/SklogWiki/index.php/Argon for argon LJ parameters
k_B = 0.0083144621*(kilojoule/(mole*kelvin))
lj_epsilon = 119.8*kelvin * k_B # in kJ/mol
lj_sigma = 0.3405*nanometer  # in nm

# All seven particles will initially be placed on a 2D plane
# Each particle is by construction at distance r_min of another atom
# r_min = lj_sigma * 2^(1./6.) is the distance at which the LJ potential is minimum
# See on the right for the numbering of particles
#
#     *   *               3   2
#   *   *   *           4   0   1
#     *   *               5   6
#
r_min = 2.**(1./6.) * lj_sigma._value
xpos  = r_min*cos(pi/3.)
ypos  = r_min*sin(pi/3.)

iniPos = np.zeros((natoms,3))
iniPos[0,:] = (center[0],center[1],center[2])
iniPos[1,:] = (center[0]+r_min,center[1],center[2])
iniPos[2,:] = (center[0]+xpos,center[1]+ypos,center[2])
iniPos[3,:] = (center[0]-xpos,center[1]+ypos,center[2])
iniPos[4,:] = (center[0]-r_min,center[1],center[2])
iniPos[5,:] = (center[0]-xpos,center[1]-ypos,center[2])
iniPos[6,:] = (center[0]+xpos,center[1]-ypos,center[2])

system = System()
for i in range(natoms):
  system.addParticle(mass)
  nb.addParticle(0.0,lj_sigma,lj_epsilon)
  
system.setDefaultPeriodicBoxVectors(a,b,c)
system.addForce(nb)
system.addForce(CMMotionRemover(100))

# generate basic topology object to use with simulation class (and to be able to write pdb files)
top = Topology()
top.setPeriodicBoxVectors([a,b,c])
chain = top.addChain("A")
for i in range(natoms):
  top.addResidue("AR",chain)
for r in top.residues():
  top.addAtom("AR",Element.getBySymbol("Ar"),r)


simulation = Simulation(top,system,integrator,platform)#,platformProperties)

simulation.context.setPositions(iniPos)

state = simulation.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
with open("start.pdb", 'w') as f:
  PDBFile.writeFile(top,state.getPositions(),f)
  
print('Minimizing energy...')
simulation.minimizeEnergy(1e-5)

simulation.context.setVelocitiesToTemperature(temperature)

state = simulation.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
with open("minim.pdb", 'w') as f:
  PDBFile.writeFile(top,state.getPositions(),f)
  
nsteps   = 1e6 
# save each 1 ps
saveFreq = 1000

simulation.reporters.append(StateDataReporter(stdout, saveFreq, step=True, 
    time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, progress=True,
    remainingTime=True, speed=True,
    totalSteps=nsteps, separator='\t'))

simulation.reporters.append(DCDReporter(file='sim.dcd', reportInterval=saveFreq, enforcePeriodicBox=True))

print("Running dynamics...")
simulation.step(nsteps)

state = simulation.context.getState(getPositions=True,getVelocities=True,enforcePeriodicBox=True)
with open("1ns.pdb", 'w') as f:
  PDBFile.writeFile(top,state.getPositions(),f)
  
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

