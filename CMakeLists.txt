cmake_minimum_required(VERSION 3.0)
project("Planck")

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    set(EIGEN3_INCLUDE_DIRS ${SOURCE_DIR})
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
    set(BOOST_INCLUDE_DIRS ${SOURCE_DIR})
endif()

# Add source files
file(GLOB EXECUTABLE_SOURCES "src/*.cpp")

# Create an executable for each source file
foreach(EXECUTABLE_SOURCE ${EXECUTABLE_SOURCES})
    get_filename_component(EXECUTABLE_NAME ${EXECUTABLE_SOURCE} NAME_WE)
    add_executable(${EXECUTABLE_NAME} ${EXECUTABLE_SOURCE})
    # Link libraries
    target_link_libraries(${EXECUTABLE_NAME} Eigen3::Eigen ${Boost_LIBRARIES})
    # Set output directory for binary
    set_target_properties(${EXECUTABLE_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endforeach()
