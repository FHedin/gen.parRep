/**
 * \file md_interface.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#include "md_interface.hpp"

using namespace std;

MD_interface::MD_interface(DATA& _dat,
                           MD_ENGINES _engine_type,
                           const string& _engine_description,
                           MD_DISTANCE_UNIT _distance_unit,
                           MD_ENERGY_UNIT _energy_unit,
                           bool _engine_supports_groups_splitting
                          ) :
                          dat(_dat),
                          engine_type(_engine_type),
                          engine_decription(_engine_description),
                          distance_unit(_distance_unit),
                          energy_unit(_energy_unit),
                          supports_groups_splitting(_engine_supports_groups_splitting)
{
}

MD_interface::~MD_interface()
{  
}

bool MD_interface::engine_supports_groups_splitting() const
{
  return supports_groups_splitting;
}
