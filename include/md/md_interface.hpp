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

#include "global.hpp"

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
  
  MD_interface(DATA& _dat) : dat(_dat) {}
  
  virtual ~MD_interface(){};

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
   * @brief copy state's data back to main programm from MD engine
   * If one of the pointer is set to nullptr, copy is ignored for the corresponding element
   * 
   * @param timeInPs pointer where to store simulation time, or nullptr
   * @param energies pointer where to store energies, or nullptr
   * @param currentTemperature pointer where to store temperature, or nullptr
   * @param atoms Atomic system coordinates, or nullptr
   * @param dat Simulation parameters
   */
  virtual void getState(double* timeInPs, ENERGIES* energies,
                        double* currentTemperature, ATOM atoms[]) = 0;
  
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
   * @return The temperature time in Kelvin
   */
  virtual double getTemperature() = 0;
  
protected:

  // general simulation parameters
  DATA& dat;
  
  // current simulation clock time
  double time;
  
  // temperature of the system
  double T;

};



#endif // MD_INTERFACE_H
