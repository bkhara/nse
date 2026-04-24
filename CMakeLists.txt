cmake_minimum_required(VERSION 3.11)

project(fracture LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules")

# Optionally let a user.cmake set MFEM_DIR (or MFEM_{RELEASE,DEBUG}_DIR)
if(EXISTS "${CMAKE_SOURCE_DIR}/user.cmake")
  include("${CMAKE_SOURCE_DIR}/user.cmake")
endif()

# Setup build type (defaults to release)
if (NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  # default build type is Release
  set(CMAKE_BUILD_TYPE "RELEASE")
endif()
if (NOT DEFINED MFEM_DIR)
# Check if MFEM_RELEASE_DIR and MFEM_DEBUG_DIR are set
if (NOT DEFINED MFEM_RELEASE_DIR OR NOT DEFINED MFEM_DEBUG_DIR)
  message(FATAL_ERROR "MFEM_RELEASE_DIR and MFEM_DEBUG_DIR must be set in user.cmake.")
endif()
# Determine the MFEM directory based on the build type
string(TOUPPER "${CMAKE_BUILD_TYPE}" _build_type)
if (_build_type STREQUAL "DEBUG")
  set(MFEM_DIR "${MFEM_DEBUG_DIR}")
elseif (_build_type STREQUAL "RELEASE")
  set(MFEM_DIR "${MFEM_RELEASE_DIR}")
else()
  message(FATAL_ERROR "Unknown build type: ${CMAKE_BUILD_TYPE}. Please set it to either 'Debug' or 'Release'.")
endif()
endif()

find_package(MPI REQUIRED COMPONENTS C CXX)
# Bring in MFEM's CMake config (defines the 'mfem' target)
#include("${MFEM_DIR}/MFEMConfig.cmake")
find_package(MFEM CONFIG REQUIRED)
## Try to find libconfig (C++ version)
find_package(Config++)
if (NOT ${CONFIG++_FOUND})
  message(FATAL_ERROR "Libconfig++ could not be located.")
endif()

# --- Core interface library (no compiled output) ---
add_library(core INTERFACE)
target_include_directories(core INTERFACE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
        ${CONFIG++_INCLUDE_DIR}
)
target_link_libraries(core INTERFACE
        mfem
        ${CONFIG++_LIBRARY}
)

# COMPILE the new source files below
set(MY_SOURCES
        src/frac.cpp
        src/FEMachinery.cpp
        src/TimeLevelFields.cpp
        src/Checkpointing.cpp
)

#foreach (src ${MY_SOURCES})
#  get_filename_component(exe "${src}" NAME_WE)
#
#  add_executable(${exe} "${src}")
#  target_link_libraries(${exe} PRIVATE core)
#endforeach ()

add_executable(frac ${MY_SOURCES})
target_link_libraries(frac PRIVATE core)
