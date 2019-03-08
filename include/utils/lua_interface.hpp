/**
 * \file lua_interface.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef LUA_INTERFACE_HPP_INCLUDED
#define LUA_INTERFACE_HPP_INCLUDED

#define SOL_CHECK_ARGUMENTS

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <functional>
#include <tuple>
#include <ostream>

#include "sol.hpp"

#include "global.hpp"
#include "md_interface.hpp"

#ifdef PARREP_DEBUG_BUILD
#define USE_SOL_PROTECTED_FUNCTIONS
#endif

#ifdef USE_SOL_PROTECTED_FUNCTIONS
typedef sol::protected_function SOL_FUNCTION;
#else
typedef sol::function SOL_FUNCTION;
#endif

// typedefs for sqlite db functions
typedef SOL_FUNCTION SQLiteDB_open;
typedef SOL_FUNCTION SQLiteDB_close;
typedef SOL_FUNCTION SQLiteDB_insert;
typedef SOL_FUNCTION SQLiteDB_backup;

/**
 * @brief a function initialising, on lua's side, all the state checking code
 */
typedef SOL_FUNCTION   ParRep_function_state_init;

/**
 * @brief a function returning a boolean describing if the system left a given state or not.
 * Definition of the state is handled by the lua function
 */
typedef SOL_FUNCTION   ParRep_function_check_state;

/** @brief a function returning a boolean describing if the system is within ametastable state
 *  All tests are performed within the lua script
 *  is true if transient propagation should be performed
 *  is false otherwise; if the configurational space has been partitionned this is always false
 */
typedef SOL_FUNCTION   ParRep_function_check_transient;

/**
 * @brief a function returning a lua table corresponding to an opaque (c++ code shouldn't touch it) serialized representation of the states
 */
typedef SOL_FUNCTION   ParRep_function_get_serialized_state;

/**
 * @brief a function returning pushing to lua a table corresponding to an opaque serialized representation of the states
 */
typedef SOL_FUNCTION   ParRep_function_put_serialized_state;

/**
 * @brief A function used for Gelman-Rubin statistics should return a double observable, without taking any parameter
 * indeed, the c++ interface exposes globally some variable to the lua interface that the user can used directly
 * more details given in the lua input script
 */
typedef SOL_FUNCTION   GR_function;

/** 
 * @brief Simply the name of a GR_function
 */
typedef std::string    GR_function_name;

// for defining a map of parsed parameters form the lua script
typedef std::string LuaParameter;
typedef std::string LuaValue;
typedef std::map<LuaParameter,LuaValue> lua_ParVal_map;

typedef std::map<std::string,std::string> omm_platform_properties;

/**
* @brief This class is in charge of parsing the lua input script, and of bridging C++ and Lua codes together.
*/
class luaInterface final
{
public:

  /**
  * @brief Constructs the Lua interface ; 
  * 
  * @param _inputFile The Lua file to be parsed
  * @param _dat Internal simulation data
  * @param _at Atomic structure
  * @param _md The MD interface
  */
  luaInterface(const std::string _inputFile,
               DATA& _dat,
               std::unique_ptr<ATOM[]>& _at,
               std::unique_ptr<MD_interface>& _md
              );
  
  /**
  * @brief This functions reads and parses the Lua input script
  * 
  */
  void parse_lua_file();

  /**
  * @brief Provides read only access to the Gelman-Rubin functions name
  * 
  * @return A const reference to a vector of names
  */
  const std::vector<GR_function_name>& get_gr_names() const {return gr_funcs_names;}
  
  /**
  * @brief Provides read only access to the Gelman-Rubin functions
  * 
  * 
  * @return A const reference to a map of (names,G-R functions)
  */
  const std::map<GR_function_name,GR_function>& get_gr_functions() const {return gr_funcs;}
  
  /**
  * @brief Returns a const reference to the Lua function initializing state definition
  * 
  * @return the reference to the function
  */
  const ParRep_function_state_init&  get_function_state_init()  const {return func_state_init;}
  
  /**
  * @brief Returns a const reference to the Lua function checking the state
  * 
  * @return the reference to the function
  */
  const ParRep_function_check_state& get_function_check_state() const {return func_check_state;}
  
  /**
   * @brief Returns a const reference to the Lua function checking if transient propagation is required
   * 
   * @return the reference to the function
   */
  const ParRep_function_check_state& get_function_check_transient() const {return func_check_transient;}
  
  /**
  * @brief Returns a const reference to the Lua function for obtaining from Lua a serialized state
  * 
  * @return the reference to the function
  */
  const ParRep_function_get_serialized_state&  get_function_get_serialized_state()  const {return func_get_serialized_state;}
  
  /**
  * @brief Returns a const reference to the Lua function for pushing to Lua a serialized state
  * 
  * @return the reference to the function
  */
  const ParRep_function_put_serialized_state&  get_function_put_serialized_state() const {return func_put_serialized_state;}
  
  /**
  * @brief Returns a const reference to the Lua function opening a SQLite DB
  * 
  * @return the reference to the function
  */
  const SQLiteDB_open&    get_db_open()   const {return db_open;}
  
  /**
  * @brief Returns a const reference to the Lua function closing a SQLite DB
  * 
  * @return the reference to the function
  */
  const SQLiteDB_close&   get_db_close()  const {return db_close;}
  
  /**
  * @brief Returns a const reference to the Lua function inserting a record in SQLite DB
  * 
  * @return the reference to the function
  */
  const SQLiteDB_insert&  get_db_insert() const {return db_insert;}
  
  /**
  * @brief Returns a const reference to the Lua function bperforming backup of a SQLite DB
  * 
  * @return the reference to the function
  */
  const SQLiteDB_backup&  get_db_backup() const {return db_backup;}
  
  /**
  * @brief Returns a const reference to the map where parsed Lua objects are stored
  * 
  * @return the reference to the map
  */
  const lua_ParVal_map& get_parsed_parameters_map() const {return pvMap;}
  
  /**
   * @brief Returns a const reference to the map where parsed OpenMM platform specific properties are stored
   * 
   * @return the reference to the map
   */
  const omm_platform_properties& get_omm_platform_properties() const {return omm_props;}

  /**
  * @brief add or update a Lua variable
  */
  template<typename T>
  void set_lua_variable(const std::string& name, T var)
  {
    lua[name] = var;
  }
  
private:

  /**
  * @brief Exposes to Lua a set of C++ lambda functions
  * 
  */
  void register_default_lua_functions();
  
  /**
  * @brief Exposes to Lua a set of C++ variables
  * 
  */
  void register_default_lua_variables();
  
  /**
   * @brief parse OpenMM MD engine parameters
   */
  void parse_omm_params();
  
  /**
  * @brief Get parrep parameters from Lua
  * 
  */
  void parse_parrep_params();
  
  /**
  * @brief Get gen. parrep parameters from Lua
  * 
  */
  void parse_parrepFV_params();
  
  /* it is crucial to have this declared first : as variable below may depend on this sol::state
   * indeed according to c++ standards it will be the last to be destroyed if it was declared first
   * see this issue : https://github.com/ThePhD/sol2/issues/69
   */
  
  sol::state lua;  ///< The core Lua object 

  std::string inputFile;                ///< path to the input file
  DATA& dat;                            ///< simulation parameters
  std::unique_ptr<ATOM[]>& at;          ///< coordinates and velocities
  std::unique_ptr<MD_interface>& md;    ///< MD interface class
//   std::string md_engine_name;
  
  std::vector<GR_function_name> gr_funcs_names;       ///< vector of Gelman-Rubin functions' names
  std::map<GR_function_name,GR_function> gr_funcs;    ///< map of Gelman-Rubin (names,functions)
  
  lua_ParVal_map pvMap;      ///< map of lua parsed parameters
  
  //a set of functions for initialisiing and using the sqlite3 database via the lua interface
  SQLiteDB_open     db_open;      ///< function to Lua SQLite open function
  SQLiteDB_close    db_close;     ///< function to Lua SQLite close function
  SQLiteDB_insert   db_insert;    ///< function to Lua SQLite insert function
  SQLiteDB_backup   db_backup;    ///< function to Lua SQLite backup function
  
  ParRep_function_state_init      func_state_init;      ///< Lua function initialising state definition
  ParRep_function_check_state     func_check_state;     ///< Lua function checking state status
  ParRep_function_check_transient func_check_transient; ///< Lua function checking if transient propagation is required
  
  ParRep_function_get_serialized_state func_get_serialized_state;   ///< Lua function for getting from Lua a serialized state
  ParRep_function_put_serialized_state func_put_serialized_state;   ///< Lua function for putting to Lua a serialized state

  omm_platform_properties omm_props;
  
  int32_t my_id     = -1;
  int32_t num_procs = -1;
};

#endif // LUA_INTERFACE_HPP_INCLUDED
