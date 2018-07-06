from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

prmtop = AmberPrmtopFile('1d7h_i.prmtop')

platform = Platform.getPlatformByName('CPU')

simulation = Simulation(topology=prmtop.topology,system="System.xml", integrator="Integrator.xml", platform=platform, state="State.xml")

eqsteps = 500000
svfreq = 500
simulation.reporters.append(StateDataReporter(stdout, svfreq, step=True, time=True, potentialEnergy=True,
                                              kineticEnergy=True, totalEnergy=True, 
                                              temperature=True, progress=True, remainingTime=True, speed=True,
                                              totalSteps=eqsteps, separator='\t'))
simulation.reporters.append(DCDReporter('run.dcd', svfreq))
simulation.step(eqsteps)

