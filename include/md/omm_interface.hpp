/**
 * \file omm_interface.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef OMM_INTERFACE_H
#define OMM_INTERFACE_H

#include <memory>

#include "md_interface.hpp"

#include "OpenMM.h"

class luaInterface;

/**
 * @brief Enumerated type representing supported OpenMM platforms
 * 
 */
enum PLATFORMS
{
  AUTO,  ///< AUTO : let OpenMM find the fastest platform
  REF,   ///< REF  : on cpu, not optimised, nor parallellised
  CPU,   ///< CPU  : on cpu, simd optimised, parallellised with threads
  OCL,   ///< OCL  : on cpu or gpu or any accelerating device available
  CUDA   ///< CUDA : on nvidia gpu only probably the fastest
};

/**
 * @brief Enumerated type representing supported OpenMM integrators
 * 
 */
enum INTEGRATORS
{
  LANGEVIN,     ///< code will use a Langevin integrator from OpenMM
  BROWNIAN      ///< code will use a Brownian integrator (i.e. overdamped Langevin) from OpenMM
};

typedef std::map<std::string,std::string> PLATFORM_PROPERTIES;

/**
 * @brief This class is a custom wrapper arround the OpenMM C++ library
 * 
 * The rest of the code interacts with this class but never directly with the OpenMM code, in order to make possible
 * a future implementation using different MD engines.
 * 
 */
class OMM_interface final : public MD_interface
{

public:
  
  /**
  * @brief The constructor ; it can be called directly, but it is easier to use the static
  *        method initialise_fromSerialization(...)
  *        taking care of XML files unserializing
  * 
  * @param _dat Common simulation parameters
  * @param omm_name_version A string representation of the enine version (e.g. "Openmm 7.3")
  * @param syst An OpenMM system
  * @param integr An OpenMM integrator, should be Langevin or Brownian (or custom version of those)
  * @param state An OpenMM state with at least positions
  * @param platformDesired The platform to use for performing the computations
  * @param platformProperties An associative map containing key-value string pairs of specific per platform parameters (to be forwarded to OpenMM)
  * 
  */
  OMM_interface(DATA& _dat,
                const std::string& omm_name_version,
                OpenMM::System& syst,
                OpenMM::Integrator& integr,
                OpenMM::State& state,
                PLATFORMS platformDesired,
                const PLATFORM_PROPERTIES& platformProperties
  );

  ~OMM_interface();
  
  virtual void doNsteps(uint32_t numSteps) override;
  
  virtual uint32_t runForPhysicalTime(std::chrono::milliseconds t) override;
  
  virtual void minimise(double tol, uint32_t maxsteps, ENERGIES& ener) override;

  virtual void minimiseWithCopy(double tol, uint32_t maxsteps, ENERGIES& ener) override;
  
  virtual void minimiseWithCopy(double tol, uint32_t maxsteps, ENERGIES& ener,
                                std::vector<XYZ>& pos, std::vector<XYZ>& vels) override;

  virtual void setSimClockTime(double t) override;
  
  virtual void setCrdsVels(ATOM at[]) override;
  
  virtual void randomiseVelocities() override;
  
  virtual void getSimData(double* timeInPs, ENERGIES* energies,
                          double* currentTemperature, ATOM atoms[], int32_t groups=-1) override;
  
  virtual void getParticlesParams(ATOM atoms[]) override;
  
  virtual double getTime() override;
  
  virtual double getTemperature() override;

  /**
   * @brief Call this for serializing this object for later doing a simulation restart
   * 
   * This is  based on context.createCheckpoint
   * 
   */
  void backupOMMobject();
  
  /**
   * @brief Call this for unserializing this object for later doing a simulation restart
   * 
   * This is  based on context.loadCheckpoint
   * 
   */
  void restoreOMMobject();
  
  ////////////////////////////////////////////
  // STATIC METHODS HERE
  ////////////////////////////////////////////

  /**
   * @brief do initialisation of OpenMM code but using previously serialized states
   * 
   * @param dat Simulation parameters
   * @param luaItf A unique_ptr to the lua interface from which parsed parameters are retrieved
   * 
   * @return OMM_interface* A pointer to an initialised OpenMM interface
   */
  static std::unique_ptr<OMM_interface> initialise_fromSerialization(DATA& dat, std::unique_ptr<luaInterface>& luaItf);
  
  static const std::map<PLATFORMS,std::string> omm_platforms_names;
  
private:
  
  /**
   * constant used when calculating temperature from kinetic energy
   * T = (2./(3.*K_Boltzmann))*Ekin
   * So EkinToK = 2./(3.*K_Boltzmann) with K_Boltzmann in units of kJ/(mol.K)
   */
  static const double EkinToK;
  
  std::unique_ptr<OpenMM::System>         system  = nullptr;
  std::unique_ptr<OpenMM::Context>        context = nullptr;
  std::unique_ptr<OpenMM::Integrator>     integrator  = nullptr;
  
  OpenMM::LangevinIntegrator* lint = nullptr;
  OpenMM::BrownianIntegrator* bint = nullptr;
  
  INTEGRATORS integType;
  
  // periodic vectors
  std::array<OpenMM::Vec3,3> pbc;
  
  std::string platformName;
  
  std::vector<OpenMM::Vec3> pos;
  std::vector<OpenMM::Vec3> vel;
  
  /// number of degreees of freedom
  int32_t dofs; 
  
  /*
   * This is used when calling runForPhysicalTime
   * The first time it is called we estimate how much of CPU physical time will cost one integration step,
   * and this is stored in cost_for_1step_dyna.
   * Then the future calls to runForPhysicalTime() will perform the required steps until the desired amount of CPU time is reached
   */
  bool first_call_of_runForPhysicalTime = true;
  std::chrono::nanoseconds cost_for_1step_dyna;
  
  /**
   * @brief Register a platform for use
   * 
   * @param platformDesired The desired platform
   * @param  properties An associative map containing key-value string pairs of specific per platform parameters (to be forwarded to OpenMM)
   * @param state A state to be used for creating the platform
   */
  void addPlatform(PLATFORMS platformDesired, const PLATFORM_PROPERTIES& properties, OpenMM::State& state);

};

#endif // OMM_INTERFACE_H
