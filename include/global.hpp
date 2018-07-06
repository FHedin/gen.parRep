/**
 * \file global.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#ifndef GLOBAL_HPP_INCLUDED
#define GLOBAL_HPP_INCLUDED

#include <array>
#include <chrono>
#include <cstdint>

//define where is the null file : should work with any unix-like OS
#ifdef __unix__
#define NULLFILE "/dev/null"
#else //for MS windows (untested, not useful ?)
#define NULLFILE "nul"
#endif

/**
 * @brief This can be used for storing either coordinates or velocities, it simply contains an array of 3 doubles.
 * 
 * It should be flexible enough for providing compatibility with various MD codes
 */
class XYZ final
{

public:

  XYZ()
  {
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  
  XYZ(double x, double y, double z)
  {
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  
  double operator[](size_t i) const
  {
    return xyz[i];
  }
  
  double& operator[](size_t i)
  {
    return xyz[i];
  }
  
  /*
   * The following are by ref redefined operator, altering the object values
   */
  XYZ& operator/=(double rhs)
  {
    xyz[0] /= rhs;
    xyz[1] /= rhs;
    xyz[2] /= rhs;
    return *this;
  }
  
  XYZ& operator*=(double rhs)
  {
    xyz[0] *= rhs;
    xyz[1] *= rhs;
    xyz[2] *= rhs;
    return *this;
  }
  
  XYZ& operator+=(const XYZ& rhs)
  {
    xyz[0] += rhs[0];
    xyz[1] += rhs[1];
    xyz[2] += rhs[2];
    return *this;
  }
  
  XYZ& operator-=(const XYZ& rhs)
  {
    xyz[0] -= rhs[0];
    xyz[1] -= rhs[1];
    xyz[2] -= rhs[2];
    return *this;
  }
  
  /*
   * The following return a new XYZ object
   */
  XYZ operator*(double rhs) const
  {
    XYZ ret(xyz[0]*rhs,xyz[1]*rhs,xyz[2]*rhs);
    return ret;
  }
  
  double x() const {return xyz[0];}
  double y() const {return xyz[1];}
  double z() const {return xyz[2];}
  
  double& x() {return xyz[0];}
  double& y() {return xyz[1];}
  double& z() {return xyz[2];}
  
  
private:
  
  std::array<double,3> xyz;

};

/**
 * \brief We also store a Periodic Boundary vector (PBC) in a similar fashion, but without redefinition of 
 *  common operators beacause we just need storage, no advanced data manipulation
 */
class PBC final
{
  
public:
  
  PBC()
  {
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  
  PBC(double x, double y, double z)
  {
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  
  double operator[](size_t i) const
  {
    return xyz[i];
  }
  
  double& operator[](size_t i)
  {
    return xyz[i];
  }
  
  double x() const {return xyz[0];}
  double y() const {return xyz[1];}
  double z() const {return xyz[2];}
  
  double& x() {return xyz[0];}
  double& y() {return xyz[1];}
  double& z() {return xyz[2];}
  
private:
  
  std::array<double,3> xyz;
};

/**
 * @brief This structure holds useful variables used across the simulations,
 * it is almost always passed by pointer/reference from one function to another one.
 */
struct DATA final
{

  uint32_t natom ;                ///< Number of atoms
  uint64_t nsteps ;               ///< Number of steps as a 64 bits integer to allow long simulations
  
  double T ;                      ///< Temperature : in Kelvin
  double timestep;                ///< Timestep for Langevin/Brownian integrator : in ps
  
  /*
   * this uses c++11 high resolution clock capabilities in order to provide accurate timers
   */
  std::chrono::time_point<std::chrono::steady_clock> start_time; ///< stores the starting time of the simulation, used for conditional stop of the simulation based on max_run_time

  std::chrono::duration<double,std::chrono::hours::period>   max_run_time_hours; ///< the maximum physical running time in hours, user defined (with default) in input file, fractional representation allowed
  
  /*
   * a duration in minutes indicating when to stop, before the max_run_time_hours: for large systems it gives time 
   * for compressing files, or deleting some of them, etc.
   * user defined (with default) in input file, should be at least one minute
   */
  std::chrono::minutes minutes_to_stop_before_max_run_time; ///< a duration in minutes indicating when to stop, before the max_run_time_hours
  
  /// a b c periodic vectors
  std::array<PBC,3> pbc;
  
  PBC& a() {return pbc[0];}
  PBC& b() {return pbc[1];}
  PBC& c() {return pbc[2];}
};


/**
 * @brief A structure representing an atom : contains coordinates (x,y,z), velocities (v,vy,vz) and mass
 */
struct ATOM final
{
  XYZ crds;  ///< OMM triplet of coordinates
  XYZ vels;  ///< OMM triplet of velocities
  double mass;        ///< atomic mass
  
  double& x(){return crds[0];} ///< accessor by reference for x coordinate
  double& y(){return crds[1];} ///< accessor by reference for y coordinate
  double& z(){return crds[2];} ///< accessor by reference for z coordinate
  
  double x() const {return crds[0];} ///< returns constant value for x coordinate
  double y() const {return crds[1];} ///< returns constant value for y coordinate
  double z() const {return crds[2];} ///< returns constant value for z coordinate

  double& vx(){return vels[0];} ///< accessor by reference for vx velocity
  double& vy(){return vels[1];} ///< accessor by reference for vy velocity
  double& vz(){return vels[2];} ///< accessor by reference for vz velocity
  
  double vx() const {return vels[0];} ///< returns constant value for vx velocity
  double vy() const {return vels[1];} ///< returns constant value for vy velocity
  double vz() const {return vels[2];} ///< returns constant value for vz velocity

};

/**
 * @brief A structure holding energy terms at a curent time of the simulation
 * The energy can be accessed either by term of as an array
 */
struct ENERGIES final
{
  ENERGIES()
  {
    ene[0] = ene[1] = ene[2] = 0.;
  }
  
  double& epot(){return ene[0];} ///< returns reference to potential energy
  double& ekin(){return ene[1];} ///< returns reference to kinetic energy
  double& etot(){return ene[2];} ///< returns reference to total energy
  
  double epot() const {return ene[0];}   ///< returns constant value of potential energy
  double ekin() const {return ene[1];}   ///< returns constant value of kinetic energy
  double etot() const {return ene[2];}   ///< returns constant value of total energy

  std::array<double,3> ene;
  
};

#endif // GLOBAL_HPP_INCLUDED
