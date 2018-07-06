/**
 * \file observable.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <numeric>
#include <algorithm>
#include <cmath>

#include "observable.hpp"
#include "mpi_utils.hpp"

using namespace std;

Observable::Observable(const string& _name){
  this->name = _name;
}

void Observable::add(const double new_obs)
{
  observations.push_back(new_obs);
}

void Observable::add(const vector<double>& new_obs)
{
  observations.insert(observations.end(),new_obs.cbegin(),new_obs.cend());
}

void Observable::add(const double new_obs[], const size_t length)
{
  for(size_t i=0; i<length; i++)
    this->add(new_obs[i]);
}

void Observable::add(const vector<double>::const_iterator& from,
                     const vector<double>::const_iterator& to)
{
  observations.insert(observations.end(),from,to);
}

double Observable::get_mean() const
{
  // accumulate using a STL function 
  const double sum = accumulate(observations.begin(), observations.end(), 0.0);
  const double mean = sum / (double) observations.size();
  return mean;
}

// static method
double Observable::get_mean(const vector<double>& obs)
{
  // accumulate using a STL function 
  const double sum = accumulate(obs.begin(), obs.end(), 0.0);
  const double mean = sum / (double) obs.size();
  return mean;
}

double Observable::get_var() const
{
  const double mean = get_mean();
  
  const size_t num_obs = observations.size();
  
  vector<double> diff(observations.size());
  // given the mean we apply a lambda function to the set of all observations in order to create a vector of deviation from the mean
  transform(observations.begin(), observations.end(), diff.begin(), [mean](double x) { return x - mean; });
  
  // calculate the sum of the square of diff members
  //   inner product == sum of products
  const double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  
  // divide by number of observations to get variance
  const double variance = sq_sum/((double)num_obs-1.0);
  
  return variance;
}

double Observable::get_sd() const
{
  // standard deviation is just square root of the variance
  const double variance = get_var();
  const double sd = sqrt(variance);
  return sd;
}

// get square of the mean
double Observable::get_mean2() const
{
  double m2 = get_mean();
  m2 *= m2;
  return m2;
}

size_t Observable::get_numObs() const
{
  return observations.size();
}

string Observable::get_name() const
{
  return name;
}

void Observable::replace_observations(const vector<double>& extern_observations)
{
  if(observations.size() == extern_observations.size())
  {
    reset_observations();
    
    //data copy
    observations = vector<double>(extern_observations);
  }
  else
  {
    fprintf(stderr,"Error in %s line %d ; see error log file for details\n",__FILE__,__LINE__);
    MPI_CUSTOM_ABORT_MACRO();
  }
}

void Observable::reset_observations()
{
  observations.clear();
}
