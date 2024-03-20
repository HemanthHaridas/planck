cmake_minimum_required(VERSION 3.0)
project("Planck")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Option to enable/disable OpenMP support
option(ENABLE_OPENMP "Enable OpenMP support" OFF)

# Check if OpenMP is requested and available
if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED)
    if(OpenMP_CXX_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Check if Eigen is available, otherwise download and build it
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
    include(ExternalProject)
    ExternalProject_Add(
        Eigen
        PREFIX ${CMAKE_BINARY_DIR}/Eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        TIMEOUT 10
        UPDATE_DISCONNECTED 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(Eigen SOURCE_DIR)
    set(EIGEN3_INCLUDE_DIR ${SOURCE_DIR})
endif()

# Check if Boost is available, otherwise download and build it
find_package(Boost QUIET COMPONENTS REQUIRED)
if(NOT Boost_FOUND)
    include(ExternalProject)
    ExternalProject_Add(
        Boost
        PREFIX ${CMAKE_BINARY_DIR}/Boost
        URL https://boostorg.jfrog.io/artifactory/main/release/1.75.0/source/boost_1_75_0.tar.gz
        URL_HASH SHA256=eb153c775fe99de5a7e21e4b57b847b2ed36f7c4c17ed0b8c1b90b6d774b80f7
        TIMEOUT 10
        UPDATE_DISCONNECTED 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ""
    )
    ExternalProject_Get_Property(Boost SOURCE_DIR)
    set(BOOST_INCLUDE_DIR ${SOURCE_DIR})
endif()

# Add the executable target
add_executable(planck src/planck.cpp)

# Include directories
target_include_directories(planck PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/auxiliary
    ${CMAKE_CURRENT_SOURCE_DIR}/src/molecule 
    ${CMAKE_CURRENT_SOURCE_DIR}/src/scf
)

# Add source files from auxiliary directory using file globbing
file(GLOB AuxiliarySources src/auxiliary/*.cpp)

# Add source files from math directory using file globbing
file(GLOB MathSources src/math/*.cpp)

# Add source files from molecule directory using file globbing
file(GLOB MoleculeSources src/molecule/*.cpp)

# Add source files from scf directory using file globbing
file(GLOB ScfSources src/scf/*.cpp)

# Add source files to the target
target_sources(planck PRIVATE
    ${MoleculeSources}
    ${ScfSources}
    ${AuxiliarySources}
    ${MathSources}
)

# Set the output directory for the planck executable
set_target_properties(planck PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# Copy the basis set directory to the bin folder
file(COPY ${CMAKE_SOURCE_DIR}/basis DESTINATION ${CMAKE_BINARY_DIR}/bin)