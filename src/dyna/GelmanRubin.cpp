/**
 * \file GelmanRubin.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <cstdio>
#include <cmath>

#include <random>
#include <functional>
#include <algorithm>

#include "GelmanRubin.hpp"
#include "logger.hpp"

using namespace std;

GelmanRubinAnalysis::GelmanRubinAnalysis(const uint32_t _num_chains, const uint32_t _discard)
{
  num_chains = _num_chains;
  discard_first = _discard;
  
  // allocate one map<string,Observable> per chain
  observablesList = unique_ptr<std::map<std::string,Observable>[]>(new map<string,Observable>[num_chains]);
  
}

// register an observable ; for all chains
void GelmanRubinAnalysis::registerObservable(const Observable& obs, const double tol)
{
  for(uint32_t n=0; n<num_chains; n++)
  {
    map<string,Observable>& chainObsList = observablesList[n];
  
    chainObsList[obs.get_name()] = obs;
  }
  
  tolerance[obs.get_name()] = tol;
  obsTypes.push_back(obs.get_name());
}

// register all observables ; for all chains
void GelmanRubinAnalysis::registerObservables(const vector<Observable>& obs, const vector<double>& tol)
{
  for(uint32_t n=0; n<num_chains; n++)
  {
    map<string,Observable>& chainObsList = observablesList[n];
    
    for(uint32_t o=0; o<obs.size(); o++)
    {
      const Observable& ob = obs[o];
      chainObsList[ob.get_name()] = ob;
      tolerance[ob.get_name()] = tol[o];
      obsTypes.push_back(ob.get_name());
    }
    
  }
}

// add 1 observation for one given chain
void GelmanRubinAnalysis::addObservation(const pair<string,double>& obs, const uint32_t chain_id)
{
  map<string,Observable>& chainObsList = observablesList[chain_id];
  const string& name = obs.first;
  const double& ob   = obs.second;
  chainObsList[name].add(ob);
}

void GelmanRubinAnalysis::addObservations(const map<string,double>& observations, const uint32_t chain_id)
{
  map<string,Observable>& chainObsList = observablesList[chain_id];
  
  for(const pair<string,double>& iter : observations)
  {
    const string& name = iter.first;
    const double& obs = iter.second;
    chainObsList[name].add(obs);
  }
}

void GelmanRubinAnalysis::addNObservations(const std::map<std::string,std::vector<double>>& n_observations,
                                           const uint32_t chain_id)
{
  map<string,Observable>& chainObsList = observablesList[chain_id];
  
  for(const pair<string,vector<double>>& iter : n_observations)
  {
    const string& name = iter.first;
    const vector<double>& obs = iter.second;
    chainObsList[name].add(obs);
  }
  
}

void GelmanRubinAnalysis::addNObservations(const std::string& obName, const std::vector<double>& observations,
                                           const uint32_t chain_id)
{
  map<string,Observable>& chainObsList = observablesList[chain_id];
  chainObsList[obName].add(observations);
}

void GelmanRubinAnalysis::addNObservations(const std::string& obName, const double observations[],
                                           const size_t length, const uint32_t chain_id)
{
  map<string,Observable>& chainObsList = observablesList[chain_id];
  chainObsList[obName].add(observations,length);
}

void GelmanRubinAnalysis::addNObservations(const std::string& obName,
                                           const std::vector<double>::const_iterator& from,
                                           const std::vector<double>::const_iterator& to,
                                           const uint32_t chain_id)
{
  map<string,Observable>& chainObsList = observablesList[chain_id];
  chainObsList[obName].add(from,to);
}

void GelmanRubinAnalysis::updateStatistics()
{
  // preallocation trick
  const Observable& fob = (*observablesList[0].cbegin()).second;
  vector<double> sqobs(fob.get_numObs());
  
  // now for each observable type update the convergence ratio
  for(string& obtype : obsTypes)
  {
    vector<double>& ratio_vector = ratio[obtype];

    double sumo2kt(0.0);
    double sumokt(0.0);
    double sumokt2(0.0);
    for(uint32_t n=0; n<num_chains; n++)
    {
      Observable& ob = observablesList[n][obtype];
      transform(ob.get_observations().begin(),ob.get_observations().end(),sqobs.begin(),[](double x){return x*x;});
      sumo2kt += Observable::get_mean(sqobs);
      sumokt  += ob.get_mean();
      sumokt2 += ob.get_mean2();
    }
    
    sumo2kt *= 1.0/num_chains;
    
    sumokt *= 1.0/num_chains;
    sumokt *= sumokt;
    
    sumokt2 *= 1.0/num_chains;
    
    const double num = sumo2kt - sumokt;
    const double den = sumo2kt - sumokt2;
    ratio_vector.push_back(num/den);

  } // loop over observables

}

void GelmanRubinAnalysis::describeChain(const uint32_t chain_id)
{
  LOG_PRINT(LOG_DEBUG,"GelmanRubinAnalysis describing chain %d (of a total of %d) :\n",chain_id,num_chains);
  
  map<string,Observable>& chainObsList = observablesList[chain_id];
  
  for(const pair<string,Observable>& pair : chainObsList)
  {
    const string& name = pair.first;
    const Observable& ob = pair.second;
    const double tol = tolerance[ob.get_name()];
    LOG_PRINT(LOG_DEBUG,"\tObservable '%s' :\n",name.c_str());
    LOG_PRINT(LOG_DEBUG,"\t\tObservations\t%lu\n",ob.get_numObs());
    LOG_PRINT(LOG_DEBUG,"\t\tTolerance\t%lf\n",tol);
    LOG_PRINT(LOG_DEBUG,"\t\tMean\t%lf\n",ob.get_mean());
    LOG_PRINT(LOG_DEBUG,"\t\tVariance\t%lf\n",ob.get_var());
    LOG_PRINT(LOG_DEBUG,"\t\tStddev.\t%lf\n",ob.get_sd());
  }
  
}

void GelmanRubinAnalysis::describe()
{
  LOG_PRINT(LOG_DEBUG,"GelmanRubinAnalysis desciption :\n");
  
  for(string& name : obsTypes)
  {
    LOG_PRINT(LOG_INFO,"\tObservable '%s' :\n",name.c_str());
    LOG_PRINT(LOG_INFO,"\t\tGR ratio\t%lf\n",get_ratio(name));
    LOG_PRINT(LOG_INFO,"\t\tConverged ?\t%s\n",(check_convergence(name))?"TRUE":"FALSE");
  }
  
}

double GelmanRubinAnalysis::get_ratio(string& obName)
{
  double lratio = ratio[obName].back();
  return lratio;
}

bool GelmanRubinAnalysis::check_convergence(const string& obName)
{
  
  vector<double>& rvec = ratio[obName];
  
  bool check;
  if(rvec.size()<discard_first)
    check = false;
  else
    check = (fabs(1.0-rvec.back()) < tolerance[obName])? true : false ;
  
  return check;
}

vector<bool> GelmanRubinAnalysis::check_convergence(const vector<string>& obList)
{
  vector<bool> check;
  
  for(const string& st : obList)
    check.push_back(check_convergence(st));
  
  return check;
}

vector<bool> GelmanRubinAnalysis::check_convergence_all()
{
  vector<bool> check;
  
  for(const string& st : obsTypes)
    check.push_back(check_convergence(st));
  
  return check;
}

void GelmanRubinAnalysis::chain_copy(const uint32_t from_id, const uint32_t to_id)
{
  map<string,Observable>& from_map = observablesList[from_id];
  map<string,Observable>& to_map   = observablesList[to_id];
  
  for(const string& st : obsTypes)
  {
    Observable& ob_from = from_map[st];
    Observable& ob_to   = to_map[st];
    ob_to.replace_observations(ob_from.get_observations());
  }
}

void GelmanRubinAnalysis::reset_all_chains()
{
  
  for(uint32_t n=0; n<num_chains; n++)
  {
    map<string,Observable>& chainMap = observablesList[n];
    
    for(const string& st : obsTypes)
    {
      Observable& ob = chainMap[st];
      ob.reset_observations();
      LOG_PRINT(LOG_DEBUG,"For chain %u, reset observable '%s'\n",n,st.c_str());
    }
    
  }
  
  for(const string& st : obsTypes)
    ratio[st].clear();
}
