from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

print("Creating System...")

prmtop = AmberPrmtopFile('1d7h_i.prmtop')
inpcrd = AmberInpcrdFile('1d7h_i.inpcrd')

# System Configuration
nonbondedMethod = CutoffNonPeriodic
#nonbondedMethod = NoCutoff
nonbondedCutoff = 1.6*nanometers
constraints = HBonds

implicitSolvent=OBC2

# those salt concentration or kappa features caused bug on openMM 7.1.1
#  use only with master git branch (future 7.2) commit >= feb79f770145bb12c99bbe058e29d077a381465d (29 th june 2017)
#  disabled for now because too slow on cpus
#saltConc=0.154*molar
#saltkappa = 1.2651*(1/nanometer)
temperature = 310.00*kelvin

system = prmtop.createSystem(
			    nonbondedMethod=nonbondedMethod,
			    nonbondedCutoff=nonbondedCutoff,
			    constraints=constraints,
			    implicitSolvent=implicitSolvent #,
			    #implicitSolventSaltConc=saltConc,
			    #implicitSolventKappa=saltkappa,
			    #temperature=temperature
			    )

print("Creating CustomCentroidBondForce...")

# we put a constraint to keep DMSO ligand around its initial distance from the pocket when equilibrating, using a CustomCentroidBondForce
force   = CustomCentroidBondForce(2,'50*max(0,distance(g1,g2)-0.36)^2')
protein = [373,375,377,380,381,388,390,392,393,693,695,697,700,711,854,856,866,868,870,872,885,886,918,920,922,925,930,937,939,940,1505,1507,1512]
dmso    = [1663,1664,1665,1666]
# make groups from atom lists
g_protein = force.addGroup(protein)
g_dmso    = force.addGroup(dmso)
#register a virtual CentroidBond between the two above defined groups
force.addBond([g_protein,g_dmso])
force_index = system.addForce(force)

print("Creating Langevin integrator...")
dt = 2*femtoseconds
friction = 2/picosecond
integrator = LangevinIntegrator(temperature, friction, dt)

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

print("Creating Platform and Simulation ...")

platform = Platform.getPlatformByName('CPU')
platformProperties = {'CpuThreads':'4'}

#platform = Platform.getPlatformByName('OpenCL')
#platformProperties = {'OpenCLPrecision':'mixed'}

#platform = Platform.getPlatformByName('CUDA')
#platformProperties = {'CudaPrecision': 'mixed'}

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(inpcrd.positions)

print('Minimizing...')
simulation.minimizeEnergy()

simulation.context.setVelocitiesToTemperature(temperature)

simulation.saveState('minim.state')
simulation.saveCheckpoint('minim.chk')

print('Equilibrating NVT for 1.0 ns with an harmonic constraint on DMSO...')

eqsteps = 500000
svfreq = 500
simulation.reporters.append(StateDataReporter(stdout, svfreq, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
                                              temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=eqsteps, separator='\t'))
simulation.reporters.append(DCDReporter('nvt.dcd', svfreq))
simulation.step(eqsteps)

simulation.saveState('nvt.state')
simulation.saveCheckpoint('nvt.chk')

simulation.currentStep = 0
simulation.context.setTime(0.0)

state = simulation.context.getState(getPositions=True, getVelocities=True,
                                    getForces=True, getEnergy=True,
                                    getParameters=True, enforcePeriodicBox=True)

#serialize the state to xml for load in parRep later
print('Serializing state : it contains coordinates, velocities, etc. for parRep...')
f = open('State.xml','w')
f.write(XmlSerializer.serialize(state))
f.close()

