/**
 * \file GelmanRubin.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef GELMANRUBIN_HPP_INCLUDED
#define GELMANRUBIN_HPP_INCLUDED

#include <vector>
#include <map>
#include <string>
#include <memory>

#include <cstdint>

#include "observable.hpp"

/**
 * @brief The Gelman and Rubin diagnostic is used to check the convergence of multiple chains run in parallel.
 * It compares the within-chain variance to the between-chain variance
 * 
 * This class represents such diagnostic, it contains methods for performing the required analyses
 * 
 * Diagnostic is performed on one or more user-defined Observable.
 */
class GelmanRubinAnalysis final
{

public:

  /**
  * @brief Constructor for a Gelman Rubin Analysis object; requires at least the number of chains (replica) as first argument.
  * 
  * @param _num_chains The number of G-R chains (same as the number of replicas usually).
  * @param _no_convergence_if_less_than Optional argument : if set there won't be convergence check until at least _no_convergence_if_less_than observations have been accumulated for each observable : this may hel avoiding pseudo-convergence.
  */
  GelmanRubinAnalysis(const uint32_t _num_chains, const uint32_t _no_convergence_if_less_than=0);
  
  
  /**
  * @brief After initialization, adds one observable and a corresponding tolerance level.
  * 
  * @param obs An Observable
  * @param tol A tolerance level for the Observable
  */
  void registerObservable(const Observable& obs, const double tol);
  
  /**
  * @brief After initialization, adds several observable and a corresponding tolerance level for each.
  * 
  * @param obs A vector of Observable to register
  * @param tol A vector of tolerance level, one for each Observable
  */
  void registerObservables(const std::vector<Observable>& obs,
                           const std::vector<double>& tol);
  

  /**
   * @brief Add one observation for one observable, from one chain
   * 
   * @param obs A pair (name,obsrvation) to add
   * @param chain_id The id of the corresponding chain
   */
  void addObservation(const std::pair<std::string,double>& obs,
                      uint32_t chain_id);
  
  /**
   * @brief Add one observation for each observable, from one chain
   * 
   * @param observations A map of (names,observations) to add
   * @param chain_id The id of the corresponding chain
   */
  void addObservations(const std::map<std::string,double>& observations,
                       const uint32_t chain_id);
  
  /**
   * @brief Add N observations for each observable, from one chain
   * 
   * @param n_observations A map of (names,vector<observations>) to add
   * @param chain_id The id of the corresponding chain
   */
  void addNObservations(const std::map<std::string,std::vector<double>>& n_observations,
                        const uint32_t chain_id);
  
  /**
   * @brief Add N observations of one observable, from one chain
   * 
   * @param obName An Observable name
   * @param observations A vector of observations
   * @param chain_id The id of the corresponding chain
   */
  void addNObservations(const std::string& obName,
                        const std::vector<double>& observations,
                        const uint32_t chain_id);
  
  /**
   * @brief Add N observations of one observable (C style array), from one chain
   * 
   * @param obName An Observable name
   * @param observations An array of observations
   * @param length Number of items in array observations
   * @param chain_id The id of the corresponding chain
   */
  void addNObservations(const std::string& obName,
                        const double observations[],
                        const size_t length,
                        const uint32_t chain_id);
  /**
   * @brief Add N observations of one observable using iterators, from one chain.
   * 
   * @param obName An Observable name
   * @param from Iterator for a vector of observations, front
   * @param to Iterator for a vector of observations, back
   * @param chain_id The id of the corresponding chain
   */
  void addNObservations(const std::string& obName,
                        const std::vector<double>::const_iterator& from,
                        const std::vector<double>::const_iterator& to,
                        const uint32_t chain_id
                       );
  
  /**
  * @brief Update G-R statistics : applies formula to check if each observable has converged
  */
  void updateStatistics();
  
  /**
  * @brief Describe G-R statistics for a given chain
  * 
  * @param chain_id The id of the corresponding chain
  */
  void describeChain(const uint32_t chain_id);
  
  /**
  * @brief Shorter version of describeChain, but describing all chains
  */
  void describe();
  
  /**
  * @brief Check if a given observable has converged
  * 
  * @param obName The observable name
  * @return Whether the observable has converged or not 
  */
  bool check_convergence(const std::string& obName);
  
  /**
  * @brief Check if several observables have converged
  * 
  * @param obList A vector of observable names to check
  * @return A vector of bool to know if the observables have converged or not
  */
  std::vector<bool> check_convergence(const std::vector<std::string>& obList);
  
  /**
  * @brief  Check if all observables have converged
  * 
  * @return A vector of bool to know if the observables have converged or not
  */
  std::vector<bool> check_convergence_all();

  /**
  * @brief Copy observables of ne chain (i.e. replica) to another one
  * 
  * @param from_id Candidate chain from which to take the values
  * @param to_id Target chain to which to copy the values
  */
  void chain_copy(const uint32_t from_id, const uint32_t to_id);
  
  /**
  * @brief Resets all chains by cleaning all the associated obervables
  * 
  */
  void reset_all_chains();
  
private:
  
  /**
  * @brief Returns the G-R ratio for a given observable
  * 
  * @param obName The observable name
  * 
  * @return The corresponding G-R ratio
  */
  double get_ratio(std::string& obName);

  std::unique_ptr<std::map<std::string,Observable>[]> observablesList; ///< Array of pointers: one for each chain
  std::vector<std::string> obsTypes;                  ///< The type (name) of all the observables
  std::map<std::string,double> tolerance;             ///< Allows one tolerance treshold per observable type
  std::map<std::string,std::vector<double>> ratio;    ///< Allows one ratio vector per observable type
  
  uint32_t num_chains = 0;                          ///< The total number of chains
  uint32_t min_num_observations_before_check = 0;   ///< There won't be convergence check until at least this amount ob observations have been accumulated for each obseravble
  
};



#endif /* GELMANRUBIN_HPP_INCLUDED */
