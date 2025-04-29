cmake_minimum_required(VERSION 3.10)
project(planck_cpp)

# Set the minimum C++ standard and use the most recent standard supported by the compiler
set(CMAKE_CXX_STANDARD 20)  # Minimum C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/opt/homebrew/Cellar/boost/1.87.0_1/include")
# set(CMAKE_CXX_EXTENSIONS OFF)  # Disable compiler-specific extensions

# Print the C++ standard being used
message(STATUS "Using C++ Standard: ${CMAKE_CXX_STANDARD}")

# Set optimization flags based on OS
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Wall -pedantic -ftree-vectorize")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin") # macOS
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -Wall -pedantic")
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++${CMAKE_CXX_STANDARD} /O2 /W4 /permissive-")
else()
    message(WARNING "Unknown system; defaulting to generic flags.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Wall -pedantic")
endif()

# if (NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE Debug) # Default to Debug if no build type is specified
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
# endif()

message(STATUS "Building on : ${CMAKE_SYSTEM_NAME} using following compiler flags : ${CMAKE_CXX_FLAGS}")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++${CMAKE_CXX_STANDARD} -O3 -Rpass=loop-vectorize -Wall -pedantic")

# find_package(Boost REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})

# FetchContent setup
include(FetchContent)
set(FETCHCONTENT_QUIET FALSE)

# Fetch Eigen
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0 # You can also use master
)
FetchContent_MakeAvailable(eigen)

# Fetch Boost
FetchContent_Declare(
  Boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-1.88.0
#   GIT_TAG b6928ae5c92e21a04bbe17a558e6e066dbe632f6
)
FetchContent_MakeAvailable(Boost)
# include_directories(${Boost_INCLUDE_DIRS})

# Gather source files
file(GLOB MAIN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
file(GLOB BASE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/base/*.cpp)
file(GLOB LOOKUP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/lookup/*.cpp)
file(GLOB HELPER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/helper/*.cpp)
file(GLOB INTEGRAL_DRIVER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/*.cpp)
file(GLOB INTEGRAL_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/*/*.cpp)
# file(GLOB TWO_ELECTRON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/two_electron/*.cpp)
file(GLOB MATH_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/math/*.cpp)

# Add the executable target
add_executable(planck 
    ${MAIN_SOURCES}
    ${BASE_SOURCES}
    ${MATH_SOURCES}
    ${HELPER_SOURCES}
    ${LOOKUP_SOURCES}
    ${INTEGRAL_SOURCES}
    ${INTEGRAL_DRIVER_SOURCES}
    # ${TWO_ELECTRON_SOURCES}
)

# Get the dependencies
set(MSYM_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym/install)
set(MSYM_BUILD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym/build)
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install)

include(ExternalProject)

ExternalProject_Add(
    libmsym
    PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/src/external/libmsym
    GIT_REPOSITORY https://github.com/mcodev31/libmsym.git
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MSYM_INSTALL_DIR} -DBUILD_SHARED_LIBS:BOOL=OFF -DMSYM_BUILD_EXAMPLES:BOOL=OFF
)

ExternalProject_Get_Property(libmsym install_dir)

# Include directories
target_include_directories(planck PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/base 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/lookup
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/helper
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/huzinaga
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/obarasakia
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals/hermite
    ${CMAKE_CURRENT_SOURCE_DIR}/src/integrals
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${MSYM_INSTALL_DIR}/include  # Include the external library headers
    ${eigen_SOURCE_DIR}
    # ${Boost_INCLUDE_DIR}
)

# set(Boost_USE_STATIC_LIBS ON) 
# set(Boost_USE_MULTITHREADED OFF)  
# set(Boost_USE_STATIC_RUNTIME OFF) 
# set(Boost_NO_SYSTEM_PATHS TRUE)

# Link the external library
target_link_libraries(planck 
                        ${MSYM_INSTALL_DIR}/lib/libmsym.a 
                        Boost::date_time 
                        Boost::property_tree 
                        Boost::foreach 
                        Boost::math
                        Boost::json
                        )

# Add dependencies to ensure the external project is built first
add_dependencies(planck libmsym eigen)

# Install the executable
install(TARGETS planck DESTINATION)

# copy the basis set files
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/basis DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# write the path of basis sets
message(STATUS "writing planck.defaults file")
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/planck.defaults "${CMAKE_INSTALL_PREFIX}/bin/basis\n")

# Build option for the separate package
option(BUILD_UTILITIES "Build the utilities" OFF)

# Conditional inclusion of the separate package
if(BUILD_UTILITIES)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/utilities)
endif()