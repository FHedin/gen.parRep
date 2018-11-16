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
find_package(LuaJIT REQUIRED)
find_package(sqlite3 REQUIRED)

# list of include folders with .h/.hpp files
include_directories(
# found by cmake
${MPI_INCLUDE_PATH}
${OpenMM_INCLUDE_DIR}
${LUA_INCLUDE_DIR}
${SQLITE3_INCLUDE_DIRS}
# project files
include
include/dyna
include/md
include/utils
# external headers (dependancies provided with this project)
external/sol2
)

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
src/md/md_interface.cpp
src/md/omm_interface.cpp
# utilities files
src/utils/mpi_utils.cpp
src/utils/lua_interface.cpp
# external dependancies provided with this project
external/luasqlite3/lsqlite3.c
)

add_executable(${TGT} ${SRCS})

set(cWarnings "-Wall -Wextra")
# more warnings for c++
# set(cxxWarnings "-Wall -Wextra -Wformat=2 -Wshadow -Wconversion -Wuseless-cast")
set(cxxWarnings "-Wall -Wextra -Wno-reorder")

set(cFeatures   "-ffast-math")
# This is required if we try to bind the executable to a binary version of the OpenMM library compiled either with  g++ < 5 or with clang
set(cxxFeatures "-ffast-math -D_GLIBCXX_USE_CXX11_ABI=0")

# C compiler flags
set(CMAKE_C_FLAGS                "-fopenmp -std=c11 ${cWarnings} ${cFeatures}")
set(CMAKE_C_FLAGS_DEBUG          "-O0 -g -DPARREP_DEBUG_BUILD")
set(CMAKE_C_FLAGS_RELEASE        "-O3    -march=native")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g -march=native")
set(CMAKE_C_FLAGS_MINSIZEREL     "-Os -s -march=native")

# C++ compiler flags
set(CMAKE_CXX_FLAGS                "-fopenmp -std=c++14 ${cxxWarnings} ${cxxFeatures}")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g -DPARREP_DEBUG_BUILD")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3    -march=native")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -march=native")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -s -march=native")

target_link_libraries(${TGT} ${OpenMM_LIBRARY} ${LUA_LIBRARIES} ${MPI_LIBRARIES} ${SQLITE3_LIBRARIES})
