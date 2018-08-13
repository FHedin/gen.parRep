# Copyright (c) 2016-2018, Florent Hédin, Tony Lelièvre, and École des Ponts - ParisTech
#All rights reserved.
#
#The 3-clause BSD license is applied to this software.
#
#See LICENSE.txt

cmake_minimum_required (VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug | Release (default) | RelWithDebInfo | MinSizeRel"
      FORCE)
endif()

set(TGT "parRep" CACHE STRING "The executable name")

project(${TGT} C CXX)

# Some of the following are required for building the project, if not found you will have to install them manually 
find_package(OpenMM REQUIRED)
find_package(MPI REQUIRED)
find_package(LuaJIT)
if(NOT LUAJIT_FOUND)
  find_package(Lua REQUIRED)
endif()
find_package(sqlite3)

# list of include folders with .h/.hpp files
include_directories(
# found by cmake
${MPI_INCLUDE_PATH}
${OpenMM_INCLUDE_DIR}
${LUA_INCLUDE_DIR}
# project files
include
include/dyna
include/md
include/utils
# external headers (dependancies provided with this project)
external/sol2
)

# if sqlite3 requires to be locally built
if(SQLITE3_FOUND)
  include_directories(${SQLITE3_INCLUDE_DIRS})
endif()

# list all source files
set(
SRCS
# core files
src/logger.cpp
src/main.cpp
src/rand.cpp
# dynamics files
src/dyna/observable.cpp
src/dyna/GelmanRubin.cpp
src/dyna/runSim.cpp
src/dyna/parRep.cpp
src/dyna/parRep_FV.cpp
src/dyna/parRep_FV_multiExits.cpp
# MD engine interfaces files
src/md/omm_interface.cpp
# utilities files
src/utils/mpi_utils.cpp
src/utils/lua_interface.cpp
# external dependancies provided with this project
external/luasqlite3/lsqlite3.c
)

if(NOT SQLITE3_FOUND)
  include_directories(external/sqlite3)
  set(SRCS ${SRCS} external/sqlite3/sqlite3.c)
endif()

add_executable(${TGT} ${SRCS})

set(cWarnings "-Wall -Wextra")
# more warnings for c++
# set(cxxWarnings "-Wall -Wextra -Wformat=2 -Wshadow -Wconversion -Wuseless-cast")
set(cxxWarnings "-Wall -Wextra")

# set(cFeatures   "-fstack-protector-all -fstack-check")
# set(cxxFeatures "-fstack-protector-all -fstack-check")

# enabling -ffast-math may provide extra performance by enabling supplementary vectorization
#set(cFeatures   "-ffast-math -fopt-info-vec")
#set(cxxFeatures "-ffast-math -fopt-info-vec")

# may be required when using clang
#set(cxxFeatures "-stdlib=libstdc++")

# C compiler flags
set(CMAKE_C_FLAGS_DEBUG          "-std=c11 ${cWarnings} ${cFeatures} -O0 -g -DPARREP_DEBUG_BUILD")
set(CMAKE_C_FLAGS_RELEASE        "-std=c11 ${cWarnings} ${cFeatures} -O3    -march=native")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-std=c11 ${cWarnings} ${cFeatures} -O2 -g -march=native")
set(CMAKE_C_FLAGS_MINSIZEREL     "-std=c11 ${cWarnings} ${cFeatures} -Os -s -march=native")

# C++ compiler flags
set(CMAKE_CXX_FLAGS_DEBUG          "-std=c++14 ${cxxWarnings} ${cxxFeatures} -O0 -g -DPARREP_DEBUG_BUILD")
set(CMAKE_CXX_FLAGS_RELEASE        "-std=c++14 ${cxxWarnings} ${cxxFeatures} -O3    -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-std=c++14 ${cxxWarnings} ${cxxFeatures} -O2 -g -march=native")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-std=c++14 ${cxxWarnings} ${cxxFeatures} -Os -s -march=native")

target_link_libraries(${TGT} ${OpenMM_LIBRARY} ${LUA_LIBRARIES} ${MPI_LIBRARIES} dl pthread)

if(SQLITE3_FOUND)
  target_link_libraries (${TGT} ${SQLITE3_LIBRARIES})
endif()

