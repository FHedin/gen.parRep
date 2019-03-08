/**
 * \file observable.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef OBSERVABLE_HPP_INCLUDED
#define OBSERVABLE_HPP_INCLUDED

#include <vector>
#include <string>

/**
 * @brief  Contains a vector of double precision observations and basic methods for analysing them
 * 
 * The is used together with the class GelmanRubin for tracking convergence to the QSD.
 * This can probably handle anything define within the lua script so we probably won't change anymore this file.
 * 
 * internally it just keeps a vector<double> and the associated name ; the access methods provide ways
 * of adding and removing obsrvations, together with analysis features (mean, variance, etc.)
 * 
 */
class Observable  final
{
public:

  /**
  * @brief Empty constructor; required for some C++ STL containers 
  * 
  */
  Observable(){}
  
  /**
  * @brief Named constructor
  * 
  * @param _name A name to assign to this observable; by default the name of the Lua function taken from the input file
  */
  Observable(const std::string& _name);

  
  /**
  * @brief Add one double precision observation to this observable
  * 
  * @param new_obs The observation to add
  */
  void add(const double new_obs);
  
  /**
  * @brief Add multiple double precision observations (stored in a vector) to this observable
  * 
  * @param new_obs The observations to add
  */
  void add(const std::vector<double>& new_obs);
  
  /**
  * @brief Add multiple double precision observations (stored in standard array) to this observable
  * 
  * @param new_obs The observations to add
  * @param length  The number of observations to add
  */
  void add(const double new_obs[], const size_t length);
  
  /**
  * @brief Add multiple double precision observations (via vector iterators from-to) to this observable
  * 
  * @param from An iterator pointing to the first element to add
  * @param to An iterator pointing to the last element to add
  */
  void add(const std::vector<double>::const_iterator& from,
           const std::vector<double>::const_iterator& to);
  
  /**
  * @brief Returns the mean over all the available observations
  * 
  * @return The mean
  */
  double get_mean() const;
  
  /**
  * @brief Returns the mean of an external vector of observations
  * @note  This is a static method
  * @return The mean of obs
  */
  static double get_mean(const std::vector<double>& obs);
  
  /**
  * @brief Returns the variance of an external vector of observations
  * 
  * @return The variance of observations
  */
  double get_var() const;
  
  /**
  * @brief Returns the standard deviation of an external vector of observations
  * 
  * @return The standard deviation of observations
  */
  double get_sd() const;
  
  /**
  * @brief Returns the square of the mean of an external vector of observations
  * 
  * @return The square of the mean of observations
  */
  double get_mean2() const;
  
  /**
  * @brief Returns the number (size) of observations
  * 
  * @return The number of registered observations
  */
  size_t get_numObs() const;
  
  /**
  * @brief Returns the name of the observable as defined during initialization
  * 
  * @return The name of the observable
  */
  std::string get_name() const;
  
  /**
  * @brief Returns a const reference to the vector of observations
  * 
  * @return const ref to observations
  */
  const std::vector<double>& get_observations() const {return observations;}
  
  /**
  * @brief Replaces the content of the vector of observations by the content of another one 
  * 
  * @param extern_observations New values to put in observations
  */
  void replace_observations(const std::vector<double>& extern_observations);
  
  /**
  * @brief Clears the vector of observations
  */
  void reset_observations();

private:

  std::string name;                   ///< The name of this observable
  std::vector<double> observations;   ///< The vector for storing the observations
  
};


#endif /* OBSERVABLE_HPP_INCLUDED */
