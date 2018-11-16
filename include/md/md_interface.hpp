/**
 * \file md_interface.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#ifndef MD_INTERFACE_H
#define MD_INTERFACE_H

#include <vector>
#include <memory>
#include <chrono>
#include <map>
#include <string>

#include "global.hpp"

/**
 * @brief List of currently supported MD engines ; 
 *        For each of them there should be an interface derived from the above defined "MD_interface" class
 * 
 * @todo TODO use it everywhere
 * 
 */
enum MD_ENGINES
{
  OPENMM,
  UNKNOWN_ENGINE
};

/**
 * @brief List of supported distance units
 * 
 * @todo TODO use it everywhere
 * 
 */
enum MD_DISTANCE_UNIT
{
  NANOMETER,
  ANGSTROEM,
  UNKNOWN_DIST_UNIT
};

/**
 * @brief List of supported energy units
 * 
 * @todo TODO use it everywhere
 * 
 */
enum MD_ENERGY_UNIT
{
  KJ_PER_MOL,
  KCAL_PER_MOL,
  UNKNOWN_ENERGY_UNIT
};

/**
 * @brief This abstract class represents an interface to a MD engine.
 * 
 * The rest of the code interacts with a class derived from this one but never directly with the MD code,
 * in order to use implementations using different MD engines.
 * 
 */
class MD_interface
{

public:
  
  /**
  * @brief Constructor to be called by derived class representing an MD engine
  * 
  * @param _dat Common simulation parameters
  * @param _engine_type The type of the engine, see the above defined MD_ENGINES enum encoding it
  * @param _engine_description A string description of the MD engine, at least a human readable name plus possible a version number at the end
  * @param _distance_unit The unit used by the engine for distances, see the above defined MD_DISTANCE_UNIT enum encoding it. Default is nanometer (nm)
  * @param _energy_unit The unit used by the engine for energies, see the above defined MD_ENERGY_UNIT enum encoding it. Default is kilojoules per mol (kJ/mol)
  * @param _engine_supports_groups_splitting A boolean encoding whether the MD engine supports retrieving simulation data for only a subset of the system (a "group"). Default is false
  */
  MD_interface(DATA& _dat,
               MD_ENGINES _engine_type,
               const std::string& _engine_description,
               MD_DISTANCE_UNIT _distance_unit = MD_DISTANCE_UNIT::NANOMETER,
               MD_ENERGY_UNIT _energy_unit = MD_ENERGY_UNIT::KJ_PER_MOL,
               bool _engine_supports_groups_splitting = false
  );
  
  /**
  * @brief Destructor of the abstract class MD_interface does nothing at the moment
  * 
  */
  virtual ~MD_interface();

  /**
   * @brief Let the MD simulation run for a given number of steps
   * 
   * @param numSteps Desired number of steps
   */
  virtual void doNsteps(uint32_t numSteps) = 0;
  
  /**
   * @brief Let the MD simulation run until at least t milliseconds of CPU time elapsed
   * 
   * This function performs numerous simulation steps until at least t milliseconds
   * of CPU time was spend in the dynamics. The exact number of simulation steps it took
   * fo reaching t milliseconds of simulaton is returned.
   * 
   * @param t Desired simulation time in milliseconds
   * 
   * @return The exact number of simulation steps performed during the time interval t
   */
  virtual uint32_t runForPhysicalTime(std::chrono::milliseconds t) = 0;
  
  /**
   * @brief Perform energy minimisation and get energies
   * 
   * @param tol Minimisation tolerance
   * @param maxsteps Maximum number of minimisation steps
   * @param ener Energy structure (by reference)
   */
  virtual void minimise(double tol, uint32_t maxsteps, ENERGIES& ener) = 0;

  /**
   * @brief Same as minimise but coordinates are copied before minimisation and restored once done
   * 
   * @param tol Minimisation tolerance
   * @param maxsteps Maximum number of minimisation steps
   * @param ener Energy structure (by reference)
   */
  virtual void minimiseWithCopy(double tol, uint32_t maxsteps, ENERGIES& ener) = 0;
  
  /**
   * @brief Same as #minimiseWithCopy but coordinates and velocities are also returned to provided vectors
   * 
   * @param tol Minimisation tolerance
   * @param maxsteps Maximum number of minimisation steps
   * @param ener Energy structure (by reference)
   * @param pos  Vector of coordinates (by reference)
   * @param vels Vector of velocities (by reference)
   */
  virtual void minimiseWithCopy(double tol, uint32_t maxsteps, ENERGIES& ener,
                                std::vector<XYZ>& pos,
                                std::vector<XYZ>& vels) = 0;

  /**
   * @brief Set the simulation clock time to a given value provided by user
   * 
   * @param t Desired time in ps
   */
  virtual void setSimClockTime(double t) = 0;
  
  /**
   * @brief Update coordinates and velocities : from main code to MD engine
   * 
   * @param at  Atomic system
   * @param dat Simulation parameters
   */
  virtual void setCrdsVels(ATOM at[]) = 0;
  
  /**
   * @brief Generate a new boltzmann distribution for velocities, for a given temperature
   * 
   * @param dat Simulation parameters
   */
  virtual void randomiseVelocities() = 0;
  
  /**
   * @brief copy simulation's data back to main programm from MD engine
   * If one of the pointer is set to nullptr, copy is ignored for the corresponding element.
   * 
   * This function implementation may support the splitting of the simulation system in different groups:
   *  if so the user may want to retrieve data only for a given group (for example the energy of a group of atoms).
   * For this the last variable groups can be used for indicating which group to retrieve energy from.
   * 
   * @param timeInPs pointer where to store simulation time, or nullptr
   * @param energies pointer where to store energies, or nullptr
   * @param currentTemperature pointer where to store temperature, or nullptr
   * @param atoms Atomic system coordinates, or nullptr
   * @param groups Which group to retrieve data for; by convention -1 means all groups. Groups splitting may be unsupported
   *               depending on the implementation
   */
  virtual void getSimData(double* timeInPs, ENERGIES* energies,
                          double* currentTemperature, ATOM atoms[],
                          int32_t groups = -1) = 0;
  
  /**
   * @brief Copy back to atoms.params LJ parameters, charge, mass from the MD code
   * 
   * @param atoms Atomic system
   */
  virtual void getParticlesParams(ATOM atoms[]) = 0;
  
  /**
   * @brief Returns current simulation time
   * 
   * @return The simulation time in ps
   */
  virtual double getTime() = 0;
  
  /**
   * @brief Returns the target temperature of the integrator
   * 
   * @return The temperature in Kelvin
   */
  virtual double getTemperature() = 0;
  
  /**
   * @brief Returns whether the md engine supports the possibility to split the system into several groups
   *        This allows retrieving simulation data using getSimData(...) for only a part of the system.
   * 
   * @return true if the feature is available, false otherwise.
   */
  bool engine_supports_groups_splitting() const;
  
protected:

  // general simulation parameters
  DATA& dat;
  
  // current simulation clock time
  double time;
  
  // temperature of the system
  double T;
  
  // pressure of the system
  double P;
  
  /**
   * @brief The type of the MD engine, encoded as an enum
   */
  const MD_ENGINES engine_type;
  
  /**
   * @brief A string description (at least the name, plus also a possible version requirement) of the available MD engines
   */
  const std::string engine_decription;
  
  /**
   * @brief The distance unit used by the MD engine ; usually either Angstroems or nanometers
   */
  const MD_DISTANCE_UNIT distance_unit;
  
  /**
   * @brief The energy unit used by the MD engine ; usually either kJ/mol or kcal/mol
   */
  const MD_ENERGY_UNIT energy_unit;
  
  /**
   * @brief set this to true (in the derived class) if the MD engine supports retrieving simulation data for only a subset of the system (a "group")
   */
  const bool supports_groups_splitting;

};



#endif // MD_INTERFACE_H
